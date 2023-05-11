from collections import Counter
from functools import cache
from logging import getLogger
from sys import byteorder
from typing import Any, Iterable, Sequence

import numpy as np
import pandas as pd

from .load import VectorLoader
from ..util.report import lookup_key
from ..util.sect import filter_gu, filter_polya, filter_pos, sects_to_pos
from ..util.seq import MATCH, DELET, INS_5, INS_3, SUB_N

logger = getLogger(__name__)

QEQ = "="
QSB = "<"
QSP = ">"
QEB = "{"
QEP = "}"
QMETHOD = QSB, QEB, QEQ, QEP, QSP


def get_mut_info_bits(muts: pd.DataFrame, refs: pd.DataFrame | None = None):
    """ Ensure muts and refs are both boolean and have the same axes,
    determine which bits are informative, mask uninformative bits in
    muts and refs to zero, and return the validated data frames of
    mutations, matches, and informative bits (all with the same shape).
    If refs is not given, initialize it to the logical not of muts. """
    # Ensure each data frame is boolean.
    mut_bits = muts.astype(bool, copy=True)
    if refs is None:
        # If no refs data frame was given, assume that it is the logical
        # not of muts, so all bits are informative.
        info_bits = pd.DataFrame(np.ones_like(mut_bits, dtype=bool),
                                 index=mut_bits.index,
                                 columns=mut_bits.columns)
    else:
        # Ensure the indexes and columns of muts and refs match.
        if not muts.index.equals(refs.index):
            raise ValueError(f"Got different indexes for muts {muts.index} "
                             f"and refs {refs.index}")
        if not muts.columns.equals(refs.columns):
            raise ValueError(f"Got different columns for muts"
                             f"{mut_bits.columns} and refs {refs.columns}")
        # Determine which bits are informative.
        info_bits = mut_bits ^ refs.astype(bool, copy=False)
        # Mask any uninformative bits in mutsb to zero.
        mut_bits &= info_bits
    # Return a boolean data frame indicating the informative bits.
    return mut_bits, info_bits


def mvec_to_bvec(vectors: pd.DataFrame,
                 query: int,
                 rel: str = QEQ) -> pd.DataFrame:
    """
    Return a boolean array of the same shape as vectors where element
    i,j is True if and only if the byte at element i,j of vectors meets
    the requested relationship between it and the query byte.

    Parameters
    ----------
    vectors: DataFrame
        Mutation vectors
    query: int
        Byte to query; must be in range 0 - 255
    rel: str = "equals"
        Method to decide whether a byte counts, as follows:
        - "=": count only bytes equal to query
        - "<": count strict bitwise subsets of query
        - ">": count strict bitwise supersets of query
        - "{": count bitwise subsets of and bytes equal to query
        - "}": count bitwise supersets of and bytes equal to query

    Returns
    -------
    DataFrame
        Boolean type DataFrame of the same shape as vectors where each
        element is True if and only if the element at the same position
        in vectors fulfilled the relationship with query
    """
    # Validate the query byte.
    if not isinstance(query, int):
        raise TypeError(
            f"Expected query of type int, but got type {type(query).__name__}")
    try:
        query.to_bytes(1, byteorder)
    except OverflowError:
        raise ValueError(f"Expected query in range 0 - 255, but got {query}")
    if rel == QEQ:
        return np.equal(vectors, query)
    if rel == QEB:
        return np.equal(np.bitwise_or(vectors, query), query)
    if rel == QEP:
        return np.equal(np.bitwise_and(vectors, query), query)
    if rel == QSB:
        return np.logical_and(mvec_to_bvec(vectors, query, QEB),
                              np.not_equal(vectors, query))
    if rel == QSP:
        return np.logical_and(mvec_to_bvec(vectors, query, QEP),
                              np.not_equal(vectors, query))
    raise ValueError(f"Parameter 'rel' must be in {QMETHOD}, but got '{rel}'")


def sum_bits(loader: VectorLoader, *,
             coords: Iterable[tuple[int, int]] = (),
             by_pos: Sequence[tuple[int, str]] = (),
             by_vec: Sequence[tuple[int, str]] = (),
             numeric: bool = False) -> dict[tuple[int, int], tuple[dict, dict]]:
    """
    For each section, count the mutations that agree with each query by
    position and each query by vector.

    Parameters
    ----------
    loader: VectorLoader
        VectorLoader from which to load the mutation vectors
    coords: Iterable[tuple[int, int]] = ()
        Iterable of 2-tuples, each defining the 5' and 3' coordinates
        of one section over which to count the bits. If empty, then use
        the entire reference as the section.
    by_pos: Sequence[tuple[int, str]] = ()
        Queries and relationships to use for counting the matching bits
        at each position of each section (also see get_bits).
    by_vec: Sequence[tuple[int, str]] = ()
        Queries and relationships to use for counting the matching bits
        within each section in each mutation vector (also see get_bits).
    numeric: bool = False
        Whether to convert the columns from base-position strings to
        numeric (specifically, integer) values of the positions, e.g.
        ['G1', 'T2', ...] if False, [1, 2, ...] if True

    Returns
    -------
    dict[tuple[int, int], tuple[dict, dict]]
        Dictionary mapping each (5', 3') coordinate pair to a tuple of
        two dictionaries of counts. The first dictionary (item 0) maps
        each key of by_vec to a Series of counts per vector (axis 0),
        indexed by the read name. The second dictionary (item 1) maps
        the keys of by_pos to a Series of counts per position (axis 1);
        its index will be the positions as integers if numeric is True,
        otherwise strings indicating the base and position.
    """
    # If no coordinates were given, then use the entire reference.
    if not coords:
        coords = [loader.section().coord]
    # For each pair of coordinates, initialize counts of the bits by
    # vector (row, axis 0) and position (column, axis 1).
    sections = list()
    counts = dict()
    for end5, end3 in coords:
        try:
            sect = loader.section(end5, end3)
        except Exception as error:
            logger.error(f"Invalid section {end5}-{end3} of {loader}: {error}")
            continue
        # Add section to list of valid sections.
        sections.append(sect)
        counts[sect.coord] = (
            # Axis 0, rows (vectors): Initialize for each query a Series
            # that will be filled with the bit count for each vector.
            {qryrel: pd.Series(dtype=int) for qryrel in by_vec},
            # Axis 1, columns (positions): Initialize for each query a
            # Series of bit counts for each position.
            {qryrel: pd.Series(np.zeros(sect.length, dtype=int),
                               index=(sect.positions if numeric
                                      else sect.columns))
             for qryrel in by_pos}
        )
    # Iterate over all the batches.
    queries_rels = set(by_pos) | set(by_vec)
    for batch in loader.iter_batches(sects_to_pos(sections),
                                     numeric=numeric):
        # Iterate over all queries and relationships.
        for qryrel in queries_rels:
            # Compute all bits in this batch.
            bits = mvec_to_bvec(batch, *qryrel)
            if qryrel in by_vec:
                # Count the bits within each section of each vector in
                # this batch, then append to the previous batches.
                for sect in sections:
                    counts[sect.coord][0][qryrel] = pd.concat([
                        counts[sect.coord][0][qryrel],
                        bits.loc[:, (sect.positions if numeric
                                     else sect.columns)].sum(axis=1)])
            if qryrel in by_pos:
                # Count the bits in this batch at each position.
                bits_sum = bits.sum(axis=0)
                # Add the bit count to each section's count.
                for sect in sections:
                    counts[sect.coord][1][qryrel] += bits_sum.loc[
                        sect.positions if numeric else sect.columns]
    return counts


class VectorFilter(object):
    """ Filter mutation vectors. """
    __slots__ = [
        "min_ninfo_pos",
        "min_fmut_pos",
        "max_fmut_pos",
        "min_mut_gap",
        "min_finfo_read",
        "max_fmut_read",
        "total_mut",
        "total_info",
        "pos_init",
        "n_min_ninfo_pos",
        "n_min_fmut_pos",
        "n_max_fmut_pos",
        "pos_kept",
        "n_reads_init",
        "n_min_finfo_read",
        "n_max_fmut_read",
        "n_min_mut_gap",
        "reads_kept",
        "_closed"
    ]

    def __init__(self, /, *,
                 min_mut_gap: int = 0,
                 min_finfo_read: float = 0.,
                 max_fmut_read: float = 1.,
                 min_ninfo_pos: int = 0,
                 min_fmut_pos: float = 0.,
                 max_fmut_pos: float = 1.):
        """
        Parameters
        ----------
        min_mut_gap: int = 0
            Filter out reads with any two mutations separated by fewer
            than ```min_mut_gap``` positions. Adjacent mutations have a
            gap of 0. If 0, keep all. Must be in [0, length_of_section).
        min_finfo_read: float = 0.0
            Filter out reads with less than this fraction of informative
            bases (i.e. match or mutation). If 0.0, keep all. Must be in
            [0, 1].
        max_fmut_read: float = 1.0
            Filter out reads with more than this fraction of mutated
            bases. If 1.0, keep all. Must be in [0, 1].
        min_ninfo_pos: int = 0
            Filter out positions with less than this number of informative
            bases. If 0, keep all. Must be ≥ 0.
        min_fmut_pos: float = 0.0
            Filter out positions with less than this fraction of mutated
            bases. If 0, keep all. Must be in [0, 1]. Should be used
            only during clustering to exclude non-useful positions with
            very few mutations.
        max_fmut_pos: float = 1.0
            Filter out positions with more than this fraction of mutated
            bases. If 1.0, keep all. Must be in [0, 1].
        """
        self.min_mut_gap = min_mut_gap
        self.min_finfo_read = min_finfo_read
        self.max_fmut_read = max_fmut_read
        self.min_ninfo_pos = min_ninfo_pos
        self.max_fmut_pos = max_fmut_pos
        self.min_fmut_pos = min_fmut_pos
        # Counts of mutations and informative bases for each position.
        self.total_mut: pd.Series | None = None
        self.total_info: pd.Series | None = None
        # Track initial and kept positions.
        self.pos_init: pd.Index | None = None
        self.pos_kept: np.ndarray | None = None
        # Track initial and kept reads.
        self.n_reads_init = 0
        self.reads_kept: Counter[str, int] = Counter()
        # Initialize counts of filtered reads and positions.
        self.n_min_finfo_read = 0
        self.n_max_fmut_read = 0
        self.n_min_mut_gap = 0
        self.n_min_ninfo_pos = 0
        self.n_max_fmut_pos = 0
        self.n_min_fmut_pos = 0
        # Indicate that the filter is open to more reads.
        self._closed = False

    def _verify_positions(self, muts: pd.DataFrame, info: pd.DataFrame):
        """ Confirm the mutations and matches data frames are consistent
        with the previously encountered data (if any). """
        if self.pos_init is None:
            # If this is the first data encountered, initialize the
            # counts of mutations and matches at each position to zero.
            self.pos_init = muts.columns
            self.total_mut = pd.Series(np.zeros(self.pos_init.size, dtype=int),
                                       index=self.pos_init)
            self.total_info = pd.Series(np.zeros(self.pos_init.size, dtype=int),
                                        index=self.pos_init)
        else:
            # Confirm that the positions of both data frames match the
            # positions of the previously encountered data.
            if not info.columns.equals(self.total_info.index):
                raise ValueError(
                    f"Positions of input reads ({muts.columns}) and existing "
                    f"reads ({self.total_mut.index}) disagree")

    def _filter_min_finfo_read(self, muts: pd.DataFrame, info: pd.DataFrame):
        """ Filter out reads with too few informative bits. """
        if not 0. <= self.min_finfo_read <= 1.:
            raise ValueError(f"min_finfo_read must be in [0, 1], but got "
                             f"{self.min_finfo_read}")
        if self.min_finfo_read > 0.:
            # Remove in-place any reads with too little information.
            read_drop = info.index[info.mean(axis=1) < self.min_finfo_read]
            muts.drop(index=read_drop, inplace=True)
            info.drop(index=read_drop, inplace=True)
            # Count the dropped reads.
            self.n_min_finfo_read += read_drop.size

    def _filter_max_fmut_read(self, muts: pd.DataFrame, info: pd.DataFrame):
        """ Filter out reads with too many mutations. """
        if not 0. <= self.max_fmut_read <= 1.:
            raise ValueError(f"max_fmut_read must be in [0, 1], but got "
                             f"{self.max_fmut_read}")
        if self.max_fmut_read < 1.:
            # Remove in-place any reads with too many mutations.
            read_drop = muts.index[muts.sum(axis=1) / info.sum(axis=1)
                                   > self.max_fmut_read]
            muts.drop(index=read_drop, inplace=True)
            info.drop(index=read_drop, inplace=True)
            # Count the dropped reads.
            self.n_max_fmut_read += read_drop.size

    def _filter_min_mut_gap(self, muts: pd.DataFrame, info: pd.DataFrame):
        """ Filter out reads with mutations that are too close. """
        if self.min_mut_gap < 0:
            raise ValueError(
                f"min_mut_gap must be ≥ 0, but got {self.min_mut_gap}")
        if self.min_mut_gap > 0:
            # Flag reads with at least one pair of mutations that are
            # too close. Initially, flag no reads (set all to False).
            flag_min_mut_gap = pd.Series(np.zeros(muts.shape[0], dtype=bool),
                                         index=muts.index)
            # Loop through every position in the bit vectors.
            for pos5 in muts.columns:
                # Find the 3'-most position that must be checked.
                pos3 = pos5 + self.min_mut_gap
                # Get all positions that remain (e.g. after dropping Gs
                # and Us) between pos5 and pos3 (including both ends).
                # For example, if positions 1, 2, 4, 5, and 7 remain and
                # min_gap = 3, then the positions checked in each
                # iteration are [1, 2, 4], [2, 4, 5], [4, 5], [4, 5, 7].
                # Then, sum over axis 1 to count the mutations in each
                # read between positions pos5 and pos3, inclusive. If a
                # read has >1 mutation in this range, then the mutations
                # are too close. Set all such reads' flags to True.
                flag_min_mut_gap |= muts.loc[:, pos5: pos3].sum(axis=1) > 1
            # Get the names of all the reads with mutations too close.
            drop = muts.index[flag_min_mut_gap]
            # Remove the reads in-place.
            muts.drop(index=drop, inplace=True)
            info.drop(index=drop, inplace=True)
            # Count the dropped reads.
            self.n_min_mut_gap += drop.size

    def _verify_no_duplicate_reads(self):
        """ Filter out reads with duplicate names. There should be none,
        but check just in case, as duplicates mess up indexing. """
        if max(self.reads_kept.values()) > 1:
            dup_reads = [read for read, count in self.reads_kept.items()
                         if count > 1]
            raise ValueError(f"Duplicate read names: {dup_reads}")

    def _filter_min_ninfo_pos(self):
        """ Filter out positions with too few informative reads. """
        if not self.min_ninfo_pos >= 0:
            raise ValueError(
                f"min_ninfo_pos must be ≥ 0, but got {self.min_mut_gap}")
        if self.min_ninfo_pos > 0:
            drop = self.total_info.index[self.total_info < self.min_ninfo_pos]
            self.total_mut.drop(index=drop, inplace=True)
            self.total_info.drop(index=drop, inplace=True)
            # Count the dropped positions.
            self.n_min_ninfo_pos += drop.size

    def _filter_max_fmut_pos(self):
        """ Filter out positions with too many mutations. """
        if not 0. <= self.max_fmut_pos <= 1.:
            raise ValueError(f"max_fmut_pos must be in [0, 1], but got "
                             f"{self.max_fmut_pos}")
        if self.max_fmut_pos < 1.:
            fmut_pos = self.total_mut / self.total_info
            drop = self.total_mut.index[fmut_pos > self.max_fmut_pos]
            self.total_mut.drop(index=drop, inplace=True)
            self.total_info.drop(index=drop, inplace=True)
            # Count the dropped positions.
            self.n_max_fmut_pos += drop.size

    def _filter_min_fmut_pos(self):
        """ Filter out positions with too few mutations. """
        if not 0. <= self.min_fmut_pos <= 1.:
            raise ValueError(f"min_fmut_pos must be in [0, 1], but got "
                             f"{self.min_fmut_pos}")
        if self.min_fmut_pos > 0.:
            fmut_pos = self.total_mut / self.total_info
            drop = self.total_mut.index[fmut_pos < self.min_fmut_pos]
            self.total_mut.drop(index=drop, inplace=True)
            self.total_info.drop(index=drop, inplace=True)
            # Count the dropped positions.
            self.n_min_fmut_pos += drop.size

    def feed(self, muts: pd.DataFrame, refs: pd.DataFrame | None = None):
        """ Feed DataFrames of mutations and (optionally) reference
        matches to the filter. """
        if self._closed:
            raise ValueError(f"{self} is closed to new reads")
        # Confirm the mutations and matches data frames are consistent
        # with each other.
        muts, info = get_mut_info_bits(muts, refs)
        # Confirm the positions are consistent with previously fed data.
        self._verify_positions(muts, info)
        # Update the total number of reads in the initial data set.
        self.n_reads_init += muts.index.size
        # Remove reads that do not pass each filter.
        self._filter_min_finfo_read(muts, info)
        self._filter_max_fmut_read(muts, info)
        self._filter_min_mut_gap(muts, info)
        # Record the names of the reads that passed the filter.
        self.reads_kept += Counter(info.index)
        # Count the mutations and information in the kept reads.
        self.total_mut += muts.sum(axis=0)
        self.total_info += info.sum(axis=0)

    def close(self):
        """ Prevent the filter from accepting more reads, and determine
        which positions pass the filters. """
        if self._closed:
            # Already closed: nothing to do.
            return
        self._closed = True
        # If no reads have been passed, then total_mut and total_info
        # will be None, so filtering positions will fail.
        if self.pos_init is None:
            raise ValueError("No reads were passed to the filter")
        # Remove positions that do not pass each filter.
        self._filter_min_ninfo_pos()
        self._filter_max_fmut_pos()
        self._filter_min_fmut_pos()
        # Determine the positions that remain after filtering.
        if not self.total_mut.index.equals(self.total_info.index):
            raise ValueError("Positions of total_mut and total_info disagree")
        self.pos_kept = self.total_mut.index.values

    @classmethod
    def get_numeric_fields(cls):
        for key in cls.__slots__:
            try:
                field = lookup_key(key)
            except KeyError:
                continue
            else:
                if field.dtype is int or field.dtype is float:
                    yield field

    @property
    def summary(self) -> dict[str, int | float]:
        self.close()
        return {field.key: self.__getattribute__(field.key)
                for field in self.get_numeric_fields()}

    @classmethod
    def null_summary(cls):
        # Initialize a filter with all default parameters.
        filt = cls()
        try:
            # Closing the filter will raise ValueError because no reads
            # have been passed.
            filt.close()
        except ValueError:
            # Ignore the error and return a summary of the empty filter.
            return filt.summary
        # Raise an error if closing worked, which should not happen.
        raise


class BitVector(object):
    """ Compute bit vectors from mutation vectors. """

    def __init__(self, /,
                 loader: VectorLoader, *,
                 end5: int,
                 end3: int,
                 count_del: bool,
                 count_ins: bool,
                 exclude_polya: int,
                 exclude_gu: bool,
                 exclude_pos: Iterable[int],
                 filter_vec: VectorFilter | None = None):
        """
        Parameters
        ----------
        loader: VectorLoader
            Mutation vector loader
        end5: int
            Position of the 5' end of the section over which to compute
            bit vectors (1-indexed). Must be ≥ 1.
        end3: int
            Position of the 3' end of the section over which to compute
            bit vectors. Must be ≥ end5 and ≤ the sequence length.
        count_del: bool
            Whether to count deletions as mutations.
        count_ins: bool
            Whether to count insertions as mutations.
        exclude_polya: int
            Exclude stretches of consecutive A bases at least this long.
            If 0, exclude no bases. Must be ≥ 0.
        exclude_gu: bool
            Whether to exclude G and U bases.
        filter_vec: VectorFilter | None = None
            Filter out low-quality or uninformative reads and positions.
        """
        self.loader = loader
        self.count_del = count_del
        self.count_ins = count_ins
        self.exclude_polya = exclude_polya
        self.exclude_gu = exclude_gu
        self.exclude_pos = list(exclude_pos)
        self.section = loader.section(end5, end3)
        self.n_pos_init = self.section.length
        # Exclude poly(A) sequences from the section.
        filt_seq_polya, filt_pos_polya = filter_polya(self.exclude_polya,
                                                      self.section.seq,
                                                      self.section.positions)
        self.n_pos_polya = self.section.positions.size - len(filt_pos_polya)
        # Exclude Gs and Us from the section.
        if self.exclude_gu:
            filt_seq_gu, filt_pos_gu = filter_gu(filt_seq_polya,
                                                 filt_pos_polya)
        else:
            filt_seq_gu, filt_pos_gu = filt_seq_polya, filt_pos_polya
        self.n_pos_gu = len(filt_pos_polya) - len(filt_pos_gu)
        # Exclude arbitrary, user-specified positions from the section.
        filt_seq_user, filt_pos_user = filter_pos(self.exclude_pos,
                                                  filt_seq_gu,
                                                  filt_pos_gu)
        self.n_pos_user = len(filt_pos_gu) - len(filt_pos_user)
        if filter_vec is None:
            # Use all the filtered positions and all the reads.
            self.positions = np.array(filt_pos_user, dtype=int)
            self.read_names = None
            # There are no filter parameters.
            self._filter_summary = VectorFilter.null_summary()
        else:
            # Determine the reads and postions to filter.
            for batch in loader.iter_batches(filt_pos_user, numeric=True):
                # Compute the mutated and matched bits in the batch and
                # feed them through the filter.
                filter_vec.feed(self.mvec_to_muts(batch),
                                self.mvec_to_refs(batch))
            filter_vec.close()
            # Get the positions and names of the reads after filtering.
            self.positions = filter_vec.pos_kept
            self.read_names = pd.Index(filter_vec.reads_kept.keys())
            # Copy the filter parameters.
            self._filter_summary = filter_vec.summary

    @property
    def n_pos_kept(self) -> int:
        return self.positions.size

    @property
    def n_reads_kept(self) -> int:
        return self.read_names.size

    @property
    def min_mut_gap(self) -> int:
        return self._filter_summary["min_mut_gap"]

    @property
    def _query_mut(self):
        query = SUB_N
        if self.count_del:
            query |= DELET
        if self.count_ins:
            query |= INS_3
        return query, QEB

    @property
    def _query_ref(self):
        return MATCH | INS_5, QEB

    def mvec_to_muts(self, mvec: pd.DataFrame):
        """ Compute bit vectors of mutations. """
        return mvec_to_bvec(mvec, *self._query_mut)

    def mvec_to_refs(self, mvec: pd.DataFrame):
        """ Compute bit vectors of reference matches. """
        return mvec_to_bvec(mvec, *self._query_ref)

    def iter_muts(self):
        """ For each batch of mutation vectors, select the reads and
        positions that passed the filters and yield a boolean data frame
        indicating the mutations. """
        # Create a Series of uninitialized floats with its index set to
        # the names of the reads. Only the index is needed, for the
        # operation pd.concat().
        read_series = pd.Series(np.empty_like(self.read_names),
                                index=self.read_names)
        # Iterate over each batch of mutation vectors.
        for batch in self.loader.iter_batches(self.positions, numeric=True):
            if self.read_names is None:
                # Yield mutations in all reads.
                yield self.mvec_to_muts(batch)
            else:
                # Find the names of the reads in the current batch
                # that passed the filter. First, create a Series whose
                # index names all reads passing the filter (read_series)
                # and another Series whose index names all reads in the
                # current batch (pd.Series(index=batch.index)). Using
                # pd.concat(join="inner"), compute a new DataFrame whose
                # index is the intersection of the indexes of those two
                # Series. Finally, the index of that DataFrame names all
                # reads in the current batch that passed the filter.
                reads_passing = pd.concat([pd.Series(index=batch.index),
                                           read_series],
                                          axis=1, join="inner").index
                # Yield mutations in the reads that passed the filter.
                yield self.mvec_to_muts(batch.loc[reads_passing])

    def all_muts(self):
        """ Return a boolean data frame indicating the mutated positions
        of every read that passed the filters. """
        return self.mvec_to_muts(
            self.loader.all_vectors(self.positions,
                                    numeric=True).loc[self.read_names])

    @cache
    def _get_unique(self):
        """ Find the unique bit vectors and return the indexes """
        # Count each unique bit vector.
        unique, inverse, counts = np.unique(self.all_muts().values, axis=0,
                                            return_inverse=True,
                                            return_counts=True)
        # For each position, find the indexes of the unique bit vectors
        # with a mutation. Storing only the mutations requires much less
        # memory than storing the entire sparse matrix (unique) because
        # mutations are relatively rare.
        mut_idxs = tuple(map(np.flatnonzero, unique.T))
        return mut_idxs, counts, inverse

    @property
    def uniq_muts(self) -> tuple[np.ndarray, ...]:
        return self._get_unique()[0]

    @property
    def uniq_counts(self) -> np.ndarray:
        return self._get_unique()[1]

    @property
    def uniq_inverse(self) -> np.ndarray:
        return self._get_unique()[2]

    @property
    def n_uniq(self):
        return self.uniq_counts.size

    def to_dict(self):
        """ Return a dictionary containing bit vector report data. """
        data: dict[str, Any] = dict()
        # Basic information about the sample and section.
        data["sample"] = self.loader.sample
        data.update(self.section.to_dict())
        # Parameters about mutations and excluded positions.
        data.update({key: self.__getattribute__(key)
                     for key in ["count_del",
                                 "count_ins",
                                 "exclude_polya",
                                 "exclude_gu",
                                 "exclude_pos"]})
        # Parameters about filtering reads and positions.
        data.update({key: value for key, value in self._filter_summary.items()
                     if key.startswith("min") or key.startswith("max")})
        # Results about given and excluded positions.
        data.update({key: self.__getattribute__(key)
                     for key in ["n_pos_init",
                                 "n_pos_polya",
                                 "n_pos_gu",
                                 "n_pos_user"]})
        # Results about filtering positions.
        data.update({key: self._filter_summary[key]
                     for key in ["n_min_ninfo_pos",
                                 "n_min_fmut_pos",
                                 "n_max_fmut_pos"]})
        data["n_pos_kept"] = self.n_pos_kept
        # Results about filtering reads.
        data.update({key: self._filter_summary[key]
                     for key in ["n_reads_init",
                                 "n_min_finfo_read",
                                 "n_max_fmut_read",
                                 "n_min_mut_gap"]})
        data["n_reads_kept"] = self.n_reads_kept
        return data

    def __str__(self):
        return f"Bit Vectors of '{self.loader.sample}' over {self.section}"
