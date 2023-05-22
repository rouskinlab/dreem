from __future__ import annotations

from abc import ABC, abstractmethod
from functools import cache, cached_property
from itertools import product
from logging import getLogger
import re
from typing import Callable, Iterable

import numpy as np
import pandas as pd

from .sect import cols_to_seq, cols_to_pos
from .seq import DNA, MATCH, DELET, INS_5, INS_3, SUB_A, SUB_C, SUB_G, SUB_T

logger = getLogger(__name__)


class SemiBitCaller(object):
    """ Convert mutation vectors (whose elements are 8-bit values that
    indicate the possible relationships between a read and a reference)
    into bit vectors (whose elements are boolean values that indicate
    whether positions are mutations or matches). """

    refs = "ACGT"
    reads = "ACGTDI"
    mut_bytes = bytes([SUB_A, SUB_C, SUB_G, SUB_T, DELET, INS_3])
    pattern_plain = re.compile(f"([{refs.lower()}])([{reads.lower()}])")
    pattern_fancy = re.compile(f"([{refs}]) -> ([{reads}])")

    @classmethod
    def _validate_refread(cls, refread: tuple[int, int]):
        ref, read = refread
        if chr(ref) not in cls.refs:
            raise ValueError(f"Invalid ref base code: {chr(ref)}")
        if chr(read) not in cls.reads:
            raise ValueError(f"Invalid read base code: {chr(read)}")
        return ref, read

    @classmethod
    def format_plain(cls, refread: tuple[int, int]):
        ref, read = cls._validate_refread(refread)
        return f"{chr(ref)}{chr(read)}".lower()

    @classmethod
    def format_fancy(cls, refread: tuple[int, int]):
        ref, read = cls._validate_refread(refread)
        return f"{chr(ref)} -> {chr(read)}"

    @classmethod
    def all_refreads(cls):
        return ((ord(ref), ord(read)) for ref, read in product(cls.refs,
                                                               cls.reads))

    @classmethod
    def refreads_to_queries(cls, refreads: Iterable[tuple[int, int]]
                            ) -> list[tuple[int, int]]:
        queries: dict[int, int] = {ord(ref): 0 for ref in cls.refs}
        for ref, read in refreads:
            if read == ref:
                # Match
                query = MATCH | INS_5
            else:
                try:
                    # Mutation
                    query = cls.mut_bytes[cls.reads.index(chr(read))]
                except ValueError:
                    logger.error(f"Invalid read code: '{chr(read)}'")
                    break
            try:
                # Update the query for the reference base.
                queries[ref] |= query
            except KeyError:
                logger.error(f"Invalid ref code: '{chr(ref)}'")
        return list(queries.items())

    @classmethod
    def query_to_refreads(cls, refs_queries: Iterable[tuple[int, int]]):
        """ For each reference base and its one-byte query, yield all
        tuples of (reference-base, read-base) that match the query. """
        for ref, query in refs_queries:
            if chr(ref) not in cls.refs:
                logger.error(f"Invalid ref code: '{chr(ref)}'")
                continue
            if query & MATCH:
                # Match
                yield ref, ref
            for mut_byte, read in zip(cls.mut_bytes, cls.reads, strict=True):
                if query & mut_byte:
                    # Mutation
                    yield ref, ord(read)

    @classmethod
    def refreads_to_codes(cls, refreads: Iterable[tuple[int, int]],
                          fancy: bool = False):
        """ For each tuple of (reference-base, read-base), yield the
        corresponding mutation code (either plain or fancy). """
        return map(cls.format_fancy if fancy else cls.format_plain, refreads)

    @classmethod
    def codes_to_refreads(cls, mut_codes: Iterable[str], fancy: bool = False):
        """ For each mutation code (either plain or fancy), yield the
        corresponding tuple of (reference-base, read-base). """
        pattern = cls.pattern_fancy if fancy else cls.pattern_plain
        for mut_code in mut_codes:
            if match := pattern.match(mut_code):
                ref, mut = match.groups()
                yield ord(ref.upper()), ord(mut.upper())
            else:
                logger.warning(f"Invalid mutation code: {mut_code}")

    def __init__(self, *codes: str):
        self.queries = self.refreads_to_queries(self.codes_to_refreads(codes))

    @cache
    def _full_query(self, seq: DNA):
        """ Convert the query dictionary into an array with one element
        per position in the sequence. """
        seqarr = np.frombuffer(seq, dtype=np.uint8)
        if seqarr.ndim != 1:
            raise ValueError(f"seq must have 1 dimension, but got {seqarr.ndim}")
        query = np.empty_like(seqarr, dtype=np.uint8)
        for base, base_query in self.queries:
            query[np.equal(seqarr, base)] = base_query
        return query.reshape((1, -1))

    def call(self, mut_vectors: pd.DataFrame) -> pd.DataFrame:
        """
        Query the given mutation vectors. Return a boolean DataFrame
        indicating which elements of mutation vectors matched the query.

        Parameters
        ----------
        mut_vectors: DataFrame
            Mutation vectors

        Returns
        -------
        DataFrame
            Boolean DataFrame of the same shape as mut_vectors wherein
            element i,j is True if element i,j of mut_vectors matches
            the query, otherwise False.
        """
        # Get the reference sequence of the mutation vectors as an array
        # of 8-bit unsigned integers.
        seq = cols_to_seq(mut_vectors.columns.to_list())
        query = self._full_query(seq)
        bits = np.equal(query, np.bitwise_or(mut_vectors, query))
        logger.debug(f"Queried mutation vectors\n{mut_vectors}\nwith\n{query}\n"
                     f"and returned\n{bits}")
        return bits

    def to_report_format(self):
        """ Return the types of counted mutations as a list. """
        return list(self.refreads_to_codes(self.query_to_refreads(self.queries),
                                           fancy=True))

    @classmethod
    def from_report_format(cls, mut_codes: Iterable[str]):
        return cls(*cls.refreads_to_codes(cls.codes_to_refreads(mut_codes,
                                                                fancy=True)))

    @classmethod
    def from_counts(cls, *,
                    count_ref: bool = False,
                    count_sub: bool = False,
                    count_del: bool = False,
                    count_ins: bool = False,
                    discount: Iterable[str] = ()):
        counts: set[str] = set()
        if count_ref:
            # Matches
            counts.update(cls.refreads_to_codes((ord(base), ord(base))
                                                for base in cls.refs))
        if count_sub:
            # Substitutions
            counts.update(cls.refreads_to_codes((ord(base1), ord(base2))
                                                for base1 in cls.refs
                                                for base2 in cls.refs
                                                if base1 != base2))
        if count_del:
            # Deletions
            counts.update(cls.refreads_to_codes((ord(ref), ord("D"))
                                                for ref in cls.refs))
        if count_ins:
            # Deletions
            counts.update(cls.refreads_to_codes((ord(ref), ord("I"))
                                                for ref in cls.refs))
        discounts = set(discount)
        if extras := discounts - set(cls.refreads_to_codes(cls.all_refreads())):
            logger.warning(f"Invalid codes of mutations to discount: {extras}")
        return cls(*(counts - discounts))


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
        return pd.Index(cols_to_pos(self.seqpos))

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
            i, m = bit_caller.call(batch)
            info.append(i)
            muts.append(m)
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


class BitVectorClusters(BitVectorSet):
    def __init__(self,
                 bit_caller: BitCaller,
                 mut_vectors: pd.DataFrame,
                 resps: pd.DataFrame):
        super().__init__(bit_caller, [mut_vectors])
        if not resps.index.equals(self.reads):
            raise ValueError(f"mut_vectors and resps have different read names "
                             f"({self.reads} ≠ {resps.index})")
        self._resps = resps

    @property
    def resps(self):
        return self._resps

    @property
    def clusters(self):
        return self.resps.columns

    @cached_property
    def nvec(self) -> pd.Series:
        """ Count the bit vectors in each cluster. """
        return self.resps.sum(axis=0)

    def _count_bits(self, bits: pd.DataFrame, by_vec: bool):
        if bits is not self._info and bits is not self._muts:
            raise ValueError(f"bits must be info or muts, but got {bits}")
        if by_vec:
            # Count the number of True bits in each vector and weight by
            # the responsibility of each read in each cluster.
            return self.resps * np.count_nonzero(bits, axis=1)
        # Find the vector and position of every True bit.
        vec_idxs, pos_idxs = np.nonzero(bits)
        # Count the bits at each position in each cluster.
        counts = pd.DataFrame(index=self.seqpos,
                              columns=self.clusters,
                              dtype=float)
        for cluster, resps in self.resps.items():
            counts[cluster] = np.bincount(pos_idxs,
                                          weights=resps[vec_idxs],
                                          minlength=self.npos)
        return counts


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
        self.nvec_given = 0
        self.nvec_filtered = list()
        filters = list(filters)
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
            for filt in filters:
                self.nvec_filtered.append(bvec.drop_vec(filt(bvec)).size)
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
        # memory than storing the entire sparse matrix (unique) because
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
