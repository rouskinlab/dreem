"""
DREEM -- Bit Caller Module
"""

from __future__ import annotations

from functools import cache as ftcache, reduce
from itertools import product
from logging import getLogger
import re
from typing import Callable, Iterable

import numpy as np
import pandas as pd

from .bitv import BitBatch
from .sect import index_to_seq
from .seq import (seq_to_int_array, DNA,
                  MATCH, DELET, INS_5, INS_3, SUB_A, SUB_C, SUB_G, SUB_T)

logger = getLogger(__name__)


class SemiBitCaller(object):
    """ Convert relation vectors (whose elements are 8-bit values that
    indicate the possible relationships between a read and a reference)
    into bit vectors (whose elements are boolean values that indicate
    whether positions are mutations or matches). """

    ref_bases = "ACGT"
    read_bases = "ACGTDI"
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
        if match := (cls.ptrn_plain.match(code.lower()) or
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

    def __init__(self, *codes: str, cache_all: bool = False):
        # Compile the codes into a query.
        self.queries = self.compile(codes)
        self._cache_all = cache_all
        self._cache = dict()
        logger.debug(f"Instantiated new {self.__class__.__name__}"
                     f"From: {codes}\nTo: {self.queries}")

    @ftcache
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

    def call(self, relvecs: pd.DataFrame, cache: bool = False) -> pd.DataFrame:
        """
        Query the given mutation vectors. Return a boolean DataFrame
        indicating which elements of mutation vectors matched the query.

        Parameters
        ----------
        relvecs: DataFrame
            Relation vectors
        cache: bool = False
            Whether to cache the result using the indexes and columns of
            `relvecs` as the keys. Assume that if `self.call()` receives
            a `relvecs` whose indexes and columns are identical to those
            of a previous call, then the previous result should be used.

        Returns
        -------
        DataFrame
            Boolean DataFrame of the same shape as mut_vectors wherein
            element i,j is True if element i,j of mut_vectors matches
            the query, otherwise False.
        """
        if cache or self._cache_all:
            # Determine the indexes and columns of relvecs.
            idxs = tuple(relvecs.index)
            cols = tuple(relvecs.columns)
            # Check if a DataFrame with those indexes and columns was
            # already cached.
            cache_key = idxs, cols
            try:
                cached_bits = self._cache[cache_key]
            except KeyError:
                # Those indexes and columns have not yet been cached.
                logger.debug(f"Found no cached bits for {relvecs} in {self}")
            else:
                # Those indexes and columns have been cached: return.
                logger.debug(f"Found cached bits for {relvecs} in {self}")
                return cached_bits
        else:
            # Do not cache the result of this call.
            cache_key = None
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
        if cache_key is not None:
            # Cache the result of this call.
            self._cache[cache_key] = bits
            logger.debug(f"Cached the bits for {relvecs} in {self}")
        return bits

    def to_report_format(self):
        """ Return the types of counted mutations as a list. """
        codes = list(self.decompile(self.queries))
        logger.debug(f"Decompiled query for {self.__class__.__name__}"
                     f"From: {self.queries}\nTo: {codes}")
        return codes

    @classmethod
    def from_report_format(cls, mut_codes: Iterable[str],
                           cache_all: bool = False):
        return cls(*list(mut_codes), cache_all=cache_all)

    @classmethod
    def from_counts(cls, *,
                    count_ref: bool = False,
                    count_sub: bool = False,
                    count_del: bool = False,
                    count_ins: bool = False,
                    discount: Iterable[str] = (),
                    cache_all: bool = False):
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
        cache_all: bool = False
            Whether to cache every invocation of `self.call()`.

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
        return cls(*codes, cache_all=cache_all)

    @classmethod
    def _junction(cls, operation: Callable, *callers: SemiBitCaller,
                  cache_all: bool = False):
        return cls.from_report_format(reduce(operation,
                                             map(set, map(cls.to_report_format,
                                                          callers))),
                                      cache_all=cache_all)

    @classmethod
    def union(cls, *callers: SemiBitCaller, cache_all: bool = False):
        """ Return the union of one or more SemiBitCallers. """
        return cls._junction(set.union, *callers, cache_all=cache_all)

    @classmethod
    def inter(cls, *callers: SemiBitCaller, cache_all: bool = False):
        """ Return the union of one or more SemiBitCallers. """
        return cls._junction(set.intersection, *callers, cache_all=cache_all)


class BitCaller(object):
    def __init__(self, nos_call: SemiBitCaller, yes_call: SemiBitCaller):
        self.nos_call = nos_call
        self.yes_call = yes_call

    def call(self, relvecs: pd.DataFrame,
             filters: dict[str, Callable[[BitBatch], pd.Index]] | None = None):
        # Using each SemiBitCaller, determine which elements match the
        # reference sequence and which are mutated.
        nos = self.nos_call.call(relvecs)
        yes = self.yes_call.call(relvecs)
        # Determine which positions are informative, which means that
        # they are unambiguously either yes or no. The corresponding
        # logical operation is exclusive or (XOR).
        info: pd.DataFrame = np.logical_xor(yes, nos)
        # For every uninformative element in info (i.e. which is False),
        # mask the element of yes at the same coordinate to False.
        yes: pd.DataFrame = np.logical_and(yes, info)
        # Create a BitBatch from the data, and optionally mask reads.
        return BitBatch(info, yes, filters)

    def iter(self, rv_batches: Iterable[pd.DataFrame],
             filters: dict[str, Callable[[BitBatch], pd.Index]] | None = None):
        """ Run `self.call()` on each batch and yield the result. """
        return (self.call(batch, filters) for batch in rv_batches)

    @classmethod
    def from_counts(cls,
                    count_del: bool = False,
                    count_ins: bool = False,
                    discount: Iterable[str] = ()):
        """ Return a new BitCaller by specifying which general types of
        mutations are to be counted, with optional ones to discount. """
        discount = list(discount)
        return cls(nos_call=SemiBitCaller.from_counts(count_ref=True,
                                                      discount=discount),
                   yes_call=SemiBitCaller.from_counts(count_sub=True,
                                                      count_del=count_del,
                                                      count_ins=count_ins,
                                                      discount=discount))
