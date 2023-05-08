from collections import Counter
from functools import cache
import json
from logging import getLogger
from sys import byteorder
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .load import VectorLoader
from ..util import path
from ..util.sect import filter_gu, filter_polya, filter_pos, sects_to_pos
from ..util.seq import MATCH, DELET, INS_5, INS_3, SUB_N

logger = getLogger(__name__)

QEQ = "="
QSB = "<"
QSP = ">"
QEB = "{"
QEP = "}"
QMETHOD = QSB, QEB, QEQ, QEP, QSP


def get_muts_refs_info(muts: pd.DataFrame, refs: pd.DataFrame | None = None):
    """ Ensure muts and refs are both boolean and have the same axes,
    determine which bits are informative, mask uninformative bits in
    muts and refs to zero, and return the validated data frames of
    mutations, matches, and informative bits (all with the same shape).
    If refs is not given, initialize it to the logical not of muts. """
    # Ensure each data frame is boolean.
    mutsb = muts.astype(bool, copy=False)
    if refs is None:
        # If no refs data frame was given, set it to the logical not of
        # the mutations data frame.
        refsb = ~mutsb
    else:
        refsb = refs.astype(bool, copy=False)
        # Ensure the indexes and columns match.
        if not mutsb.index.equals(refsb.index):
            raise ValueError(f"Got different indexes for muts {mutsb.index} "
                             f"and refs {refsb.index}")
        if not mutsb.columns.equals(refsb.columns):
            raise ValueError(f"Got different columns for muts"
                             f"{mutsb.columns} and refs {refsb.columns}")
    # Determine which bits are informative.
    info_bits = mutsb ^ refsb
    # Mask any uninformative bits in muts and refs to zero.
    mutsb &= info_bits
    refsb &= info_bits
    # Return a boolean data frame indicating the informative bits.
    return mutsb, refsb, info_bits


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

    def __init__(self, /, *,
                 min_mut_gap: int,
                 min_finfo_read: float,
                 max_fmut_read: float,
                 min_ninfo_pos: int,
                 max_fmut_pos: float,
                 min_fmut_pos: float = 0.):
        """
        Parameters
        ----------
        min_mut_gap: int
            Filter out reads with any two mutations separated by fewer
            than ```min_mut_gap``` positions. Adjacent mutations have a
            gap of 0. If 0, keep all. Must be in [0, length_of_section).
        min_finfo_read: float
            Filter out reads with less than this fraction of informative
            bases (i.e. unambiguous match or mutation). If 0, keep all.
            Must be in [0, 1].
        max_fmut_read: float
            Filter out reads with more than this fraction of mutated
            bases (i.e. unambiguous match or mutation). If 1, keep all.
            Must be in [0, 1].
        min_ninfo_pos: int
            Filter out positions with less than this number of informative
            bases. If 0, keep all. Must be ≥ 0.
        max_fmut_pos: float
            Filter out positions with more than this fraction of mutated
            bases. If 0, keep all. Must be in [0, 1].
        min_fmut_pos: float = 0.0
            Filter out positions with less than this fraction of mutated
            bases. If 0, keep all. Must be in [0, 1]. Should be used
            only during clustering to exclude non-useful positions with
            very few mutations.
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
        self.pos_all = None
        # Track reads.
        self.n_reads_all = 0
        self.reads_use: Counter[str, int] = Counter()
        # Initialize counts of filtered reads and positions.
        self.n_min_finfo_read = 0
        self.n_max_fmut_read = 0
        self.n_min_mut_gap = 0
        self.n_min_ninfo_pos = 0
        self.n_max_fmut_pos = 0
        self.n_min_fmut_pos = 0
        # Indicate that the filter is open to more reads.
        self._closed = False

    def feed(self, muts: pd.DataFrame, refs: pd.DataFrame | None = None):
        """ Feed DataFrames of mutations and (optionally) reference
        matches to the filter. """
        if self._closed:
            raise ValueError(f"{self} is closed to new reads")
        # Confirm the mutations and matches data frames are consistent
        # with each other.
        muts, refs, info = get_muts_refs_info(muts, refs)
        self.n_reads_all += info.index.size
        # Confirm the mutations and matches data frames are consistent
        # with the previously encountered data (if any).
        if self.pos_all is None:
            # If this is the first data encountered, initialize the
            # counts of mutations and matches at each position to zero.
            self.pos_all = info.columns
            self.total_mut = pd.Series(np.zeros(self.pos_all.size, dtype=int),
                                       index=self.pos_all)
            self.total_info = pd.Series(np.zeros(self.pos_all.size, dtype=int),
                                        index=self.pos_all)
        else:
            # Confirm that the positions of both data frames match the
            # positions of the previously encountered data.
            if not info.columns.equals(self.total_info.index):
                raise ValueError(
                    f"Positions of input reads ({muts.columns}) and existing "
                    f"reads ({self.total_info.index}) disagree")

        # Filter out reads with too few informative bits.
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

        # Filter out reads with too many mutations.
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

        # Filter out reads with mutations that are too close.
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
            read_drop = muts.index[flag_min_mut_gap]
            # Remove the reads in-place.
            muts.drop(index=read_drop, inplace=True)
            info.drop(index=read_drop, inplace=True)
            # Count the dropped reads.
            self.n_min_mut_gap += read_drop.size

        # Filter out reads with duplicate names (there should be none,
        # but check just in case, as duplicates could mess up indexing).
        self.reads_use += Counter(info.index)
        if max(self.reads_use.values()) > 1:
            dup_reads = [read for read, count in self.reads_use.items()
                         if count > 1]
            raise ValueError(f"Duplicate read names: {dup_reads}")

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
        if self.pos_all is None:
            raise ValueError("No reads were passed to the filter")

        # Filter out positions with too few informative reads.
        if not self.min_ninfo_pos >= 0:
            raise ValueError(
                f"min_ninfo_pos must be ≥ 0, but got {self.min_mut_gap}")
        if self.min_ninfo_pos > 0:
            pos_drop = self.total_info.index[self.total_info
                                             < self.min_ninfo_pos]
            self.total_mut.drop(pos_drop)
            self.total_info.drop(pos_drop)
            # Count the dropped positions.
            self.n_min_ninfo_pos += pos_drop.size

        # Filter out positions with too many mutations.
        if not 0. <= self.max_fmut_pos <= 1.:
            raise ValueError(f"max_fmut_pos must be in [0, 1], but got "
                             f"{self.max_fmut_pos}")
        if self.max_fmut_pos < 1.:
            pos_drop = self.total_info.index[self.total_mut / self.total_info
                                             > self.max_fmut_pos]
            self.total_mut.drop(pos_drop)
            self.total_info.drop(pos_drop)
            # Count the dropped positions.
            self.n_max_fmut_pos += pos_drop.size

        # Filter out positions with too few mutations.
        if not 0. <= self.min_fmut_pos <= 1.:
            raise ValueError(f"min_fmut_pos must be in [0, 1], but got "
                             f"{self.min_fmut_pos}")
        if self.min_fmut_pos > 0.:
            pos_drop = self.total_info.index[self.total_mut / self.total_info
                                             < self.min_fmut_pos]
            self.total_mut.drop(pos_drop)
            self.total_info.drop(pos_drop)
            # Count the dropped positions.
            self.n_min_fmut_pos += pos_drop.size

    @property
    def n_pos_all(self) -> int:
        self.close()
        return self.pos_all.size

    @property
    def pos_kept(self) -> np.ndarray:
        self.close()
        return self.total_info.index.values

    @property
    def n_pos_kept(self) -> int:
        return self.pos_kept.size

    @property
    def n_reads_use(self) -> int:
        self.close()
        return len(self.reads_use)

    @property
    def summary_pos(self):
        self.close()
        return {
            "Positions given": self.n_pos_all,
            "Positions cut -- too few informative bits": self.n_min_ninfo_pos,
            "Positions cut -- too few mutations": self.n_min_fmut_pos,
            "Positions cut -- too many mutations": self.n_max_fmut_pos,
            "Positions kept": self.n_pos_kept,
        }

    @property
    def summary_reads(self):
        self.close()
        return {
            "Reads given": self.n_reads_all,
            "Reads cut -- too few informative bits": self.n_min_finfo_read,
            "Reads cut -- too many mutations": self.n_max_fmut_read,
            "Reads cut -- mutations too close": self.n_min_mut_gap,
            "Reads kept": self.n_reads_use,
        }

    @property
    def summary_params(self):
        return {
            "Per position -- minimum informative bits": self.min_ninfo_pos,
            "Per position -- minimum mutation fraction": self.min_fmut_pos,
            "Per position -- maximum mutation fraction": self.max_fmut_pos,
            "Per read -- minimum informative bit fraction": self.min_finfo_read,
            "Per read -- maximum mutation fraction": self.max_fmut_read,
            "Per read -- minimum gap between mutations": self.min_mut_gap,
        }


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
            self._filter_params = dict()
            self._filter_pos = dict()
            self._filter_reads = dict()
            self.min_gap = 0
        else:
            # Determine the reads and postions to filter.
            for batch in loader.iter_batches(filt_pos_user, numeric=True):
                # Compute the mutated and matched bits in the batch and
                # feed them through the filter.
                filter_vec.feed(self.mvec_to_muts(batch),
                                self.mvec_to_refs(batch))
            # Get the positions and names of the reads after filtering.
            self.positions = filter_vec.pos_kept
            self.read_names = pd.Series(filter_vec.reads_use)
            # Copy the filter parameters.
            self._filter_params = filter_vec.summary_params
            self._filter_pos = filter_vec.summary_pos
            self._filter_pos.pop("Positions given")
            self._filter_reads = filter_vec.summary_reads
            self.min_gap = filter_vec.min_mut_gap
        # Write the bit vector generation report.
        self.write_report()

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
        for batch in self.loader.iter_batches(self.positions, numeric=True):
            if self.read_names is None:
                # Select the positions passing the filters.
                mvec = batch.loc[:, self.positions]
            else:
                # Find which reads in this batch passed the filters.
                reads = pd.concat([pd.Series(index=batch.index),
                                   self.read_names],
                                  axis=1, join="inner").index
                # Select the reads and positions passing the filters.
                mvec = batch.loc[reads, self.positions]
            # Convert the batch of mutation vectors to bit vectors.
            yield self.mvec_to_muts(mvec)

    def all_muts(self):
        """ Return a boolean data frame indicating the mutated positions
        of every read that passed the filters. """
        return pd.concat(self.iter_muts(), axis=0)

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

    def report_path(self):
        """ Path to the bit vector report. """
        return path.build(path.ModSeg,
                          path.SampSeg,
                          path.RefSeg,
                          path.SectSeg,
                          path.VecRepSeg,
                          top=self.loader.out_dir,
                          module=path.MOD_VECT,
                          sample=self.loader.sample,
                          ref=self.loader.ref,
                          end5=self.section.end5,
                          end3=self.section.end3,
                          ext=path.JSON_EXT)

    @property
    def report_data(self):
        """ Return a dictionary containing bit vector report data. """
        # Basic information about the sample and section.
        data = {
            "Sample": self.loader.sample,
            "Reference": self.loader.ref,
            "5' coord": self.section.end5,
            "3' coord": self.section.end3,
            "Sequence": self.section.seq.decode(),
        }
        # Parameters about mutations and excluded positions.
        data.update({
            "Count deletions": self.count_del,
            "Count insertions": self.count_ins,
            "Exclude poly(A) sequences of length": self.exclude_polya,
            "Exclude G/U bases": self.exclude_gu,
            "Exclude user-defined positions": self.exclude_pos,
        })
        # Parameters about filtering reads and positions.
        data.update(self._filter_params)
        # Results about given and excluded positions.
        data.update({
            "Positions given": len(self.section.seq),
            "Positions cut -- poly(A) sequences": self.n_pos_polya,
            "Positions cut -- G/U bases": self.n_pos_gu,
            "Positions cut -- user-defined": self.n_pos_user
        })
        # Results about filtering positions.
        data.update(self._filter_pos)
        # Results about filtering reads.
        data.update(self._filter_reads)
        return data

    def write_report(self):
        """ Publish the bit vector report in a JSON file. """
        report_path = self.report_path()
        report_path.parent.mkdir(parents=True, exist_ok=True)
        with open(report_path, "w") as f:
            json.dump(self.report_data, f)
        logger.info(f"Wrote report of {self} to {report_path}")

    def __str__(self):
        return (f"Bit Vectors of '{self.loader.sample}' aligned to "
                f"'{self.loader.ref}': {self.section.end5}-{self.section.end3}")
