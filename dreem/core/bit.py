from __future__ import annotations

from abc import ABC, abstractmethod
from functools import cache, cached_property
from itertools import product
from logging import getLogger
import re
from typing import Callable, Iterable

import numpy as np
import pandas as pd

from .sect import index_to_seq, index_to_pos
from .seq import (seq_to_int_array, DNA,
                  MATCH, DELET, INS_5, INS_3, SUB_A, SUB_C, SUB_G, SUB_T)

logger = getLogger(__name__)


class SemiBitCaller(object):
    """ Convert mutation vectors (whose elements are 8-bit values that
    indicate the possible relationships between a read and a reference)
    into bit vectors (whose elements are boolean values that indicate
    whether positions are mutations or matches). """

    ref_bases = "ACGT"
    read_bases = "ACGTDI"
    refreads = list(product(ref_bases, read_bases))
    mut_bits = bytes([SUB_A, SUB_C, SUB_G, SUB_T, DELET, INS_5 | INS_3])
    fmt_plain = "{}{}"
    fmt_fancy = "{} -> {}"
    ptrn_plain = re.compile(f"([{ref_bases.lower()}])([{read_bases.lower()}])")
    ptrn_fancy = re.compile(f"([{ref_bases}]) -> ([{read_bases}])")

    @classmethod
    def as_match(cls, code: str) -> re.Match[str]:
        """
        Return a re.Match object if the code matches either the plain or
        the fancy format. Otherwise, raise ValueError.
        """
        # If code matches ptrn_plain, cls.ptrn_plain.match(code.lower())
        # is truthy, so short-circuit the OR and return the plain match.
        # If code matches ptrn_fancy, cls.ptrn_fancy.match(code.upper())
        # is truthy, so match becomes truthy and is returned.
        if match := (cls.ptrn_plain.match(code.lower())
                     or
                     cls.ptrn_fancy.match(code.upper())):
            return match
        raise ValueError(f"Failed to match code: '{code}'")

    @classmethod
    def as_plain(cls, code: str):
        """
        Convert a ref-read code into plain format, as follows:

        - 2-character lowercase string
        - 1st character is reference base
        - 2nd character is read base, 'd' (deletion), or 'i' (insertion)

        Examples:

        - 'ag' means an A in the reference is a G in the read
        - 'cd' means a C in the reference is deleted in the read
        """
        return cls.fmt_plain.format(*cls.as_match(code).groups()).lower()

    @classmethod
    def as_fancy(cls, code: str):
        """
        Convert a ref-read code into fancy format, as follows:

        - 6-character uppercase string
        - 1st character is reference base
        - 6th character is read base, 'D' (deletion), or 'I' (insertion)
        - 2nd to 5th characters form an arrow: ' -> '

        Examples:

        - 'A -> G' means an A in the reference is a G in the read
        - 'C -> D' means a C in the reference is deleted in the read
        """
        return cls.fmt_fancy.format(*cls.as_match(code).groups()).upper()

    @classmethod
    def compile(cls, codes: Iterable[str]):
        """
        Given one or more codes in plain or fancy format, return a dict
        that maps each reference base to a query byte that will match
        all and only the codes given for that reference base.

        This function is the inverse of `cls.decompile`.
        """
        # Create a dict that maps each reference base to a query byte,
        # which is an integer in the range [0, 256). Initialize to 0.
        queries: dict[str, int] = {ref: 0 for ref in cls.ref_bases}
        # For each code given, get the ref and read bases by converting
        # the code to plain format, then to uppercase, then to a tuple.
        for ref, read in map(str.upper, map(cls.as_plain, codes)):
            # Update the query byte for the reference base. If the read
            # and reference bases are equal, then this code represents
            # a match, so update using the match byte, MATCH.
            # Otherwise, update using the mutation bit that corresponds
            # to the read base (it is at the same index in cls.mut_bytes
            # as the read base is in cls.read_bases). Update by taking
            # the bitwise OR so that all query bytes are accumulated.
            queries[ref] |= (MATCH if read == ref
                             else cls.mut_bits[cls.read_bases.index(read)])
        return queries

    @classmethod
    def decompile(cls, queries: dict[str, int]):
        """
        For each reference base and its one-byte query, yield all codes
        that the query will count.

        This function is the inverse of `cls.compile`.
        """
        # Check each query byte.
        for ref, query in queries.items():
            if ref not in cls.ref_bases:
                raise ValueError(f"Invalid reference base: '{ref}'")
            if query & MATCH:
                # The query byte has its match bit set to 1, so the code
                # in which this ref base matches the read base counts.
                yield cls.as_fancy(f"{ref}{ref}")
            # For each mutation bit, check whether the query has the bit
            # set to 1.
            for mut_bit, read in zip(cls.mut_bits, cls.read_bases, strict=True):
                if query & mut_bit:
                    # If the mutation bit is set to 1 in the query byte,
                    # then the code where the ref base becomes the read
                    # base (or deletion/insertion) counts.
                    yield cls.as_fancy(f"{ref}{read}")

    def __init__(self, *codes: str):
        # Compile the codes into a query.
        self.queries = self.compile(codes)
        logger.debug(f"Instantiated new {self.__class__.__name__}"
                     f"From: {codes}\nTo: {self.queries}")

    @cache
    def seq_query(self, seq: DNA):
        """ Convert the query dictionary into an array with one element
        per position in the sequence. """
        # Cast the sequence from DNA (subclass of bytes) to a NumPy array.
        seqarr = seq_to_int_array(seq)
        if seqarr.ndim != 1:
            raise ValueError(f"seq must have 1 dimension, but got {seqarr.ndim}")
        # Initialize an empty query array: one element per base in seq.
        query = np.zeros_like(seqarr, dtype=np.uint8)
        # Set the elements of the query array for each type of ref base.
        for ref_base, ref_query in self.queries.items():
            # Locate all elements of seq with the given type of ref base
            # and set the corresponding positions in query to the query
            # byte for the ref base.
            query[np.equal(seqarr, ord(ref_base))] = ref_query
        return query

    def call(self, relvecs: pd.DataFrame) -> pd.DataFrame:
        """
        Query the given mutation vectors. Return a boolean DataFrame
        indicating which elements of mutation vectors matched the query.

        Parameters
        ----------
        relvecs: DataFrame
            Relation vectors

        Returns
        -------
        DataFrame
            Boolean DataFrame of the same shape as mut_vectors wherein
            element i,j is True if element i,j of mut_vectors matches
            the query, otherwise False.
        """
        # Get the query array for the sequence of the relation vectors.
        query = self.seq_query(index_to_seq(relvecs.columns))
        # Determine whether each element of the relation vectors counts
        # given the query. A byte counts if and only if it equals or is
        # a bitwise subset of the query byte. This check is implemented
        # by computing the bitwise OR of the query byte and the byte in
        # the relation vector. The bitwise OR will equal the query byte
        # if the relation vector byte is a bitwise subset of the query.
        # If not, it will have bits set to 1 that are not in the query.
        bits = np.equal(query, np.bitwise_or(relvecs, query))
        logger.debug(f"Queried mutation vectors\n{relvecs}\nwith\n{query}\n"
                     f"and returned\n{bits}")
        return bits

    def to_report_format(self):
        """ Return the types of counted mutations as a list. """
        codes = list(self.decompile(self.queries))
        logger.debug(f"Decompiled query for {self.__class__.__name__}"
                     f"From: {self.queries}\nTo: {codes}")
        return codes

    @classmethod
    def from_report_format(cls, mut_codes: Iterable[str]):
        return cls(*list(mut_codes))

    @classmethod
    def from_counts(cls, *,
                    count_ref: bool = False,
                    count_sub: bool = False,
                    count_del: bool = False,
                    count_ins: bool = False,
                    discount: Iterable[str] = ()):
        """
        Return a new SemiBitCaller by specifying which general types of
        relationships are to be counted.

        Parameters
        ----------
        count_ref: bool = False
            Whether to call True all matches between the read and ref.
        count_sub: bool = False
            Whether to call True all substitutions in the read.
        count_del: bool = False
            Whether to call True all deletions in the read.
        count_ins: bool = False
            Whether to call True all insertions in the read.
        discount: Iterable[str] = ()
            Do not count any of these relationships between the read and
            the reference, even if they would be counted according to
            any of the other parameters. Should be an iterable of str in
            either plain or fancy format (except case-insensitive).

        Returns
        -------
        SemiBitCaller
            New SemiBitCaller instance that counts the specified bytes.
        """
        codes: set[str] = set()
        if count_ref:
            # Count all matches between the read and reference.
            codes.update(cls.as_plain(2 * base) for base in cls.ref_bases)
        if count_sub:
            # Count all substitutions in the read.
            codes.update(cls.as_plain(f"{base1}{base2}")
                         for base1, base2 in product(cls.ref_bases, repeat=2)
                         if base1 != base2)
        if count_del:
            # Count all deletions in the read.
            codes.update(cls.as_plain(f"{base}D") for base in cls.ref_bases)
        if count_ins:
            # Count all insertions in the read.
            codes.update(cls.as_plain(f"{base}I") for base in cls.ref_bases)
        # Remove all the codes to be discounted.
        codes -= set(map(cls.as_plain, discount))
        logger.debug(f"Converted counts for {cls.__name__}\n"
                     f"ref: {count_ref}\nsub: {count_sub}\n"
                     f"del: {count_del}\nins: {count_ins}\n"
                     f"dis: {discount}\nTo codes: {sorted(codes)}")
        return cls(*codes)


class BitCaller(object):
    def __init__(self, ref_caller: SemiBitCaller, mut_caller: SemiBitCaller):
        self.ref_caller = ref_caller
        self.mut_caller = mut_caller

    def call(self, mut_vectors: pd.DataFrame):
        """
        Query the given mutation vectors. Return two boolean DataFrames
        indicating, respectively, which elements of the mutation vectors
        are informative and which are mutated.

        Parameters
        ----------
        mut_vectors: DataFrame
            Mutation vectors

        Returns
        -------
        tuple[pd.DataFrame, pd.DataFrame]
            The informative and mutated positions, respectively.
        """
        # Using each SemiBitCaller, determine which elements match the
        # reference sequence and which are mutated.
        refs = self.ref_caller.call(mut_vectors)
        muts = self.mut_caller.call(mut_vectors)
        # Determine which positions are informative, which means that
        # they are unambiguously either matches or mutations. Logically,
        # this condition corresponds to an exclusive or (XOR) operation.
        info: pd.DataFrame = np.logical_xor(refs, muts)
        # For every uninformative element in info (i.e. which is False),
        # mask the element of muts at the same coordinate to False.
        muts: pd.DataFrame = np.logical_and(info, muts)
        return info, muts

    @classmethod
    def from_counts(cls,
                    count_del: bool = False,
                    count_ins: bool = False,
                    discount: Iterable[str] = ()):
        discount = list(discount)
        return cls(ref_caller=SemiBitCaller.from_counts(count_ref=True,
                                                        discount=discount),
                   mut_caller=SemiBitCaller.from_counts(count_sub=True,
                                                        count_del=count_del,
                                                        count_ins=count_ins,
                                                        discount=discount))


class BitVectorBase(ABC):
    def __init__(self, _: BitCaller, __: Iterable[pd.DataFrame]):
        pass

    @property
    @abstractmethod
    def reads(self):
        """ Read names. """
        return pd.Index()

    @property
    @abstractmethod
    def seqpos(self):
        """ Column names. """
        return pd.Index()

    @cached_property
    def positions(self):
        """ Numeric positions. """
        return pd.Index(index_to_pos(self.seqpos))

    @property
    @abstractmethod
    def nvec(self):
        """ Number of bit vectors. """
        return 0

    @property
    def npos(self):
        """ Number of positions in each bit vector. """
        return self.seqpos.size

    @property
    @abstractmethod
    def ninfo_per_pos(self):
        return pd.Series()

    @property
    @abstractmethod
    def ninfo_per_vec(self):
        return pd.Series()

    @property
    @abstractmethod
    def nmuts_per_pos(self):
        return pd.Series()

    @property
    @abstractmethod
    def nmuts_per_vec(self):
        return pd.Series()

    @property
    def finfo_per_pos(self) -> pd.Series:
        """ Fraction of informative bits for each position. """
        return self.ninfo_per_pos / self.nvec

    @property
    def finfo_per_vec(self) -> pd.Series:
        """ Fraction of informative bits for each vector. """
        return self.ninfo_per_vec / self.npos

    @property
    def fmuts_per_pos(self) -> pd.Series:
        """ Fraction of mutated bits for each position. """
        return self.nmuts_per_pos / self.ninfo_per_pos

    @property
    def fmuts_per_vec(self) -> pd.Series:
        """ Fraction of mutated bits for each vector. """
        return self.nmuts_per_vec / self.ninfo_per_vec

    def check_duplicate_reads(self):
        """ If a read name occurs more than once, raise ValueError. """
        if self.reads.has_duplicates:
            dups = self.reads[self.reads.duplicated(keep="first")]
            raise ValueError(f"Duplicate read names: {dups}")


class BitVectorSet(BitVectorBase):
    """ Collects and operates on one or more batches of bit vectors. """

    def __init__(self,
                 bit_caller: BitCaller,
                 mut_vectors: Iterable[pd.DataFrame]):
        super().__init__(bit_caller, mut_vectors)
        # Call the informative and mutated bits in every batch.
        info: list[pd.DataFrame] = list()
        muts: list[pd.DataFrame] = list()
        for batch in mut_vectors:
            if info:
                # Verify that the new columns match the previous ones.
                if not batch.columns.equals(info[0].columns):
                    raise ValueError("New batch positions differ from previous "
                                     f"({batch.columns} ≠ {info[0].columns})")
            info_batch, muts_batch = bit_caller.call(batch)
            info.append(info_batch)
            muts.append(muts_batch)
        if len(info) == 0:
            raise ValueError("No mutation vectors were passed")
        if len(info) == 1:
            # Only one batch of mutation vectors was passed.
            self._info = info[0]
            self._muts = muts[0]
        else:
            # Merge the informative bits and mutated bits into data frames.
            self._info = pd.concat(info, axis=0)
            self._muts = pd.concat(muts, axis=0)
        self.check_duplicate_reads()

    @property
    def info(self):
        return self._info

    @property
    def muts(self):
        return self._muts

    @property
    def reads(self) -> pd.Index:
        """ Read names. """
        return self.info.index

    @property
    def seqpos(self):
        """ Column names. """
        return self.info.columns

    @property
    def nvec(self):
        return self.reads.size

    def _count_bits(self, bits: pd.DataFrame, by_vec: bool):
        if bits is not self._info and bits is not self._muts:
            raise ValueError(f"bits must be info or muts, but got {bits}")
        # Count the number of True bits at each position (if by_vec is
        # False) or in each vector (if by_vec is True).
        return pd.Series(np.count_nonzero(bits, axis=int(by_vec)),
                         self.reads if by_vec else self.seqpos)

    @property
    def ninfo_per_pos(self):
        """ Number of informative bits for each position. """
        return self._count_bits(self.info, by_vec=False)

    @property
    def ninfo_per_vec(self):
        """ Number of informative bits for each vector. """
        return self._count_bits(self.info, by_vec=True)

    @property
    def nmuts_per_pos(self):
        """ Number of mutated bits for each position. """
        return self._count_bits(self.muts, by_vec=False)

    @property
    def nmuts_per_vec(self):
        """ Number of mutated bits for each vector. """
        return self._count_bits(self.muts, by_vec=True)

    def drop_vec(self, to_drop: np.ndarray):
        """ Remove bit vectors in-place. """
        read_names = self.reads[to_drop]
        self._info.drop(index=read_names, inplace=True)
        self._muts.drop(index=read_names, inplace=True)
        return read_names

    def get_unique_muts(self):
        return UniqMutBits(self.muts)


class BitCounter(BitVectorBase):
    def __init__(self,
                 bit_caller: BitCaller,
                 mut_vectors: Iterable[pd.DataFrame], *,
                 filters: Iterable[Callable[[BitVectorSet], np.ndarray]] = ()):
        super().__init__(bit_caller, mut_vectors)
        self._ninfo_per_pos: pd.Series | None = None
        self._nmuts_per_pos: pd.Series | None = None
        self._ninfo_per_vec = list()
        self._nmuts_per_vec = list()
        filters = list(filters)
        self.nvec_filtered = [0] * len(filters)
        self.nvec_given = 0
        for batch in mut_vectors:
            if self._ninfo_per_pos is None:
                # Initialize the numbers of informative and mutated bits
                # to zero for each position.
                self._ninfo_per_pos = pd.Series(np.zeros(batch.shape[1],
                                                         dtype=int),
                                                index=batch.columns)
                self._nmuts_per_pos = pd.Series(np.zeros(batch.shape[1],
                                                         dtype=int),
                                                index=batch.columns)
            else:
                # Verify that the new positions match the previous ones.
                if not batch.columns.equals(self._ninfo_per_pos.index):
                    raise ValueError(f"Inconsistent positions: {batch.columns} "
                                     f"≠ {self._ninfo_per_pos.index})")
            # Create a BitVectorSet from this batch.
            bvec = BitVectorSet(bit_caller, [batch])
            self.nvec_given += bvec.reads.size
            # Filter the bit vectors in this batch.
            for i, filt in enumerate(filters):
                self.nvec_filtered[i] += bvec.drop_vec(filt(bvec)).size
            # Add the count of each to the total per position.
            self._ninfo_per_pos += np.count_nonzero(bvec.info, axis=0)
            self._nmuts_per_pos += np.count_nonzero(bvec.muts, axis=0)
            # Accumulate the count of each per vector.
            self._ninfo_per_vec.append(pd.Series(np.count_nonzero(bvec.info,
                                                                  axis=1),
                                                 index=bvec.reads))
            self._nmuts_per_vec.append(pd.Series(np.count_nonzero(bvec.muts,
                                                                  axis=1),
                                                 index=bvec.reads))
        if self._ninfo_per_pos is None or self._nmuts_per_pos is None:
            raise ValueError("No mutation vectors were passed")
        self.check_duplicate_reads()

    @property
    def ninfo_per_pos(self) -> pd.Series:
        return self._ninfo_per_pos

    @property
    def ninfo_per_vec(self) -> pd.Series:
        return pd.concat(self._ninfo_per_vec, axis=0)

    @property
    def nmuts_per_pos(self) -> pd.Series:
        return self._nmuts_per_pos

    @property
    def nmuts_per_vec(self) -> pd.Series:
        return pd.concat(self._nmuts_per_vec, axis=0)

    @property
    def reads(self):
        """ Read names. """
        return self.ninfo_per_vec.index

    @property
    def nvec(self):
        return sum(batch.size for batch in self._ninfo_per_vec)

    @property
    def read_batches(self):
        """ Read names in each batch. """
        return (batch.index for batch in self._ninfo_per_vec)

    @property
    def seqpos(self):
        """ Column names. """
        return self._ninfo_per_pos.index


class UniqMutBits(object):
    def __init__(self, muts: np.ndarray):
        """ Find the unique bit vectors. """
        if muts.ndim != 2:
            raise ValueError(f"muts needs 2 dimensions, but got {muts.ndim}")
        # Count each unique bit vector.
        uniq, self._inv, self._count = np.unique(muts, axis=0,
                                                 return_inverse=True,
                                                 return_counts=True)
        # For each position, find the indexes of the unique bit vectors
        # with a mutation. Storing only the mutations requires much less
        # memory than storing the entire sparse matrix (uniq) because
        # mutations are relatively rare.
        self._uniq_idxs = tuple(map(np.flatnonzero, uniq.T))

    @property
    def indexes(self):
        """ For each position, the indexes of all unique bit vectors
        that have a mutation at that position. """
        return self._uniq_idxs

    @property
    def counts(self) -> np.ndarray:
        """ Number of times each unique bit vector occurs in the given
        set of possibly redundant bit vectors. """
        return self._count

    @property
    def inverse(self) -> np.ndarray:
        """ Indexes to map the unique bit vectors back to the given set
        of possibly redundant bit vectors. """
        return self._inv

    @property
    def n_uniq(self):
        """ Number of unique bit vectors. """
        return self.counts.size

    @property
    def n_pos(self):
        """ Number of positions in each bit vector. """
        return len(self.indexes)

    def get_full(self):
        """ Full boolean matrix of the unique bit vectors. """
        # Initialize an all-False matrix with one row for each unique
        # bit vector and one column for each position.
        full = np.zeros((self.n_uniq, self.n_pos), dtype=bool)
        # For each position (j), set the mutated elements to True.
        for j, indexes in enumerate(self.indexes):
            full[indexes, j] = True
        return full
