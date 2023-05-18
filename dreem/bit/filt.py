from collections import Counter
from typing import Any

import numpy as np
import pandas as pd

from .info import get_mut_info_bits


class VectorFilter(object):
    """ Filter bit vectors. """
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
        "pos_min_ninfo",
        "pos_min_fmut",
        "pos_max_fmut",
        "pos_kept",
        "n_reads_init",
        "n_reads_min_finfo",
        "n_reads_max_fmut",
        "n_reads_min_gap",
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
        # Track initial, kept, and filtered positions.
        self.pos_init: list[int] | None = None
        self.pos_min_ninfo: list[int] = list()
        self.pos_min_fmut: list[int] = list()
        self.pos_max_fmut: list[int] = list()
        self.pos_kept: list[int] = list()
        # Track initial, kept, and filtered reads.
        self.n_reads_init = 0
        self.n_reads_min_finfo = 0
        self.n_reads_max_fmut = 0
        self.n_reads_min_gap = 0
        self.reads_kept: Counter[str, int] = Counter()
        # Indicate that the filter is open to more reads.
        self._closed = False

    def _verify_positions(self, muts: pd.DataFrame, info: pd.DataFrame):
        """ Confirm the mutations and matches data frames are consistent
        with the previously encountered data (if any). """
        if self.pos_init is None:
            # If this is the first data encountered, initialize the
            # counts of mutations and matches at each position to zero.
            self.pos_init = muts.columns.to_list()
            self.total_mut = pd.Series(np.zeros(len(self.pos_init), dtype=int),
                                       index=self.pos_init)
            self.total_info = pd.Series(np.zeros(len(self.pos_init), dtype=int),
                                        index=self.pos_init)
        else:
            # Confirm that the positions of both data frames match the
            # positions of the previously encountered data.
            if not muts.columns.equals(self.total_mut.index):
                raise ValueError(
                    f"Positions of input reads ({muts.columns}) and existing "
                    f"reads ({self.total_mut.index}) disagree")
            if not info.columns.equals(self.total_info.index):
                raise ValueError(
                    f"Positions of input reads ({info.columns}) and existing "
                    f"reads ({self.total_info.index}) disagree")

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
            self.n_reads_min_finfo += read_drop.size

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
            self.n_reads_max_fmut += read_drop.size

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
            read_drop = muts.index[flag_min_mut_gap]
            # Remove the reads in-place.
            muts.drop(index=read_drop, inplace=True)
            info.drop(index=read_drop, inplace=True)
            # Count the dropped reads.
            self.n_reads_min_gap += read_drop.size

    def _verify_no_duplicate_reads(self):
        """ Check for reads with duplicate names because such duplicates
        would break the name-based indexing scheme. """
        if max(self.reads_kept.values()) > 1:
            dup_reads = [read for read, count in self.reads_kept.items()
                         if count > 1]
            raise ValueError(f"Duplicate read names: {dup_reads}")

    def _filter_min_ninfo_pos(self):
        """ Filter out positions with too few informative reads. """
        if not self.min_ninfo_pos >= 0:
            raise ValueError(
                f"min_ninfo_pos must be ≥ 0, but got {self.min_mut_gap}")
        pos_drop = self.total_info.index[self.total_info < self.min_ninfo_pos]
        self.total_mut.drop(index=pos_drop, inplace=True)
        self.total_info.drop(index=pos_drop, inplace=True)
        self.pos_min_ninfo = pos_drop.to_list()

    def _filter_min_fmut_pos(self):
        """ Filter out positions with too few mutations. """
        if not 0. <= self.min_fmut_pos <= 1.:
            raise ValueError(f"min_fmut_pos must be in [0, 1], but got "
                             f"{self.min_fmut_pos}")
        fmut_pos = self.total_mut / self.total_info
        pos_drop = self.total_mut.index[fmut_pos < self.min_fmut_pos]
        self.total_mut.drop(index=pos_drop, inplace=True)
        self.total_info.drop(index=pos_drop, inplace=True)
        self.pos_min_fmut = pos_drop.to_list()

    def _filter_max_fmut_pos(self):
        """ Filter out positions with too many mutations. """
        if not 0. <= self.max_fmut_pos <= 1.:
            raise ValueError(f"max_fmut_pos must be in [0, 1], but got "
                             f"{self.max_fmut_pos}")
        fmut_pos = self.total_mut / self.total_info
        pos_drop = self.total_mut.index[fmut_pos > self.max_fmut_pos]
        self.total_mut.drop(index=pos_drop, inplace=True)
        self.total_info.drop(index=pos_drop, inplace=True)
        self.pos_max_fmut = pos_drop.to_list()

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
        self._filter_min_fmut_pos()
        self._filter_max_fmut_pos()
        # Determine the positions that remain after filtering.
        if not self.total_mut.index.equals(self.total_info.index):
            raise ValueError("Positions of total_mut and total_info disagree")
        self.pos_kept = self.total_mut.index.to_list()

    @property
    def n_pos_min_ninfo(self):
        self.close()
        return len(self.pos_min_ninfo)

    @property
    def n_pos_min_fmut(self):
        self.close()
        return len(self.pos_min_fmut)

    @property
    def n_pos_max_fmut(self):
        self.close()
        return len(self.pos_max_fmut)

    @property
    def n_pos_kept(self):
        self.close()
        return len(self.pos_kept)

    @property
    def n_read_kept(self):
        self.close()
        return len(self.reads_kept)

    @staticmethod
    def summary_keys():
        return ["min_mut_gap", "min_finfo_read", "max_fmut_read",
                "min_ninfo_pos", "min_fmut_pos", "max_fmut_pos",
                "n_reads_init", "n_reads_min_finfo",
                "n_reads_max_fmut", "n_reads_min_gap",
                "n_pos_min_ninfo", "n_pos_min_fmut", "n_pos_max_fmut",
                "pos_min_ninfo", "pos_min_fmut", "pos_max_fmut"]

    @property
    def summary(self) -> dict[str, Any]:
        self.close()
        return {key: self.__getattribute__(key) for key in self.summary_keys()}

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
        raise RuntimeError(f"Failed to detect empty {cls.__name__}")
