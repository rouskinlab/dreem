from functools import cached_property
import json
from logging import getLogger
from pathlib import Path

import numpy as np
import pandas as pd

from ..calc.count import get_bits, QEB
from ..util import path
from ..util.sect import filter_gu, filter_polya
from ..util.seq import MATCH, DELET, INS_5, INS_3, SUB_N
from ..vector.load import VectorLoader

logger = getLogger(__name__)


class BitVector(object):
    """ Compute bit vectors from mutation vectors. """

    def __init__(self, /,
                 loader: VectorLoader, *,
                 end5: int,
                 end3: int,
                 max_polya: int,
                 include_gu: bool,
                 include_del: bool,
                 include_ins: bool,
                 max_muts_per_read: int,
                 min_mut_dist: int,
                 min_reads: int,
                 info_thresh: float,
                 signal_thresh: float):
        logger.info(f"Began computing bit vectors for {loader}, "
                    f"section {end5}-{end3}")
        self.end5 = end5
        self.end3 = end3
        self.max_polya = max_polya
        self.include_gu = include_gu
        self.include_del = include_del
        self.include_ins = include_ins
        self.max_muts_per_read = max_muts_per_read
        self.min_mut_dist = min_mut_dist
        self.min_cov = min_reads
        self.min_info = info_thresh
        self.min_sig = signal_thresh
        # Collect basic information from the report.
        self.sample = loader.sample
        self.ref = loader.ref
        self.seq = loader.section(self.end5, self.end3).seq
        self.n_pos_all = len(self.seq)
        # Find the positions of the vectors to load.
        self.positions = self._get_positions(loader)
        # Load the mutation vectors of the reads.
        reads = loader.get_all_vectors(positions=self.positions, numeric=True)
        # Numbers of positions and reads loaded.
        self.n_read_try, self.n_pos_try = reads.shape
        if self.n_pos_try != self.positions.size:
            raise ValueError(f"Number of positions loaded ({self.n_pos_try}) "
                             f"≠ number requested ({self.positions.size})")
        # Compute the bit vectors from the reads' mutation vectors.
        refs, muts = self._get_bit_vectors(reads)
        # Remove uninformative reads with too few informative bits.
        self.n_low_info = self._filter_info(refs, muts)
        # Remove low-quality reads with mutations too close.
        self.n_near_muts = self._filter_mut_dist(refs, muts)
        # Remove low-quality reads with too many mutations.
        self.n_many_muts = self._filter_max_muts(refs, muts)
        # Remove uninformative positions with low coverage.
        self.n_low_cov = self._filter_coverage(refs, muts)
        # Remove uninformative positions with low signal.
        self.n_low_signal = self._filter_signal(refs, muts)
        # Numbers of positions and mutation vectors after filtering.
        self.n_read_use, self.n_pos_use = muts.shape
        # Check that the numbers of positions and reads are consistent.
        if self.n_pos_try != (self.n_pos_use
                              + self.n_low_cov
                              + self.n_low_signal):
            raise ValueError(f"Number of positions loaded ({self.n_pos_try}) "
                             f"≠ used ({self.n_pos_use}) "
                             f"+ low-coverage ({self.n_low_cov}) "
                             f"+ low signal ({self.n_low_signal})")
        if self.n_read_try != (self.n_read_use
                               + self.n_low_info
                               + self.n_near_muts
                               + self.n_many_muts):
            raise ValueError(f"Number of reads loaded ({self.n_read_try}) "
                             f"≠ used ({self.n_read_use}) "
                             f"+ low-info ({self.n_low_info}) "
                             f"+ close mutations ({self.n_near_muts}) "
                             f"+ too many mutations ({self.n_many_muts})")
        # Save the read names.
        self.read_names = muts.index
        # Count each unique bit vector.
        bvec, self.bvec_to_read, self.counts = np.unique(muts.values.T, axis=1,
                                                         return_inverse=True,
                                                         return_counts=True)
        n_pos_use, self.n_bvec_use = bvec.shape
        if n_pos_use != self.n_pos_use:
            raise ValueError(f"Number of positions before ({self.n_pos_use}) "
                             f"and after ({n_pos_use}) finding unique bit "
                             f"vectors are not equal")
        logger.debug(f"Found {self.n_bvec_use} unique bit vectors for {self}:\n"
                     f"{bvec}")
        # For each position, find the indexes of the unique bit vectors
        # with a mutation. Storing only the mutations requires much less
        # memory than storing the complete sparse matrix (bvec) because
        # mutations are relatively rare.
        self.muts = tuple(map(np.flatnonzero, bvec))
        logger.debug(f"Found mutated bits at each position for {self}")

    def _get_positions(self, loader: VectorLoader):
        """ Filter out positions with poly(A) sequences that are too
        long or with G and U (or T) bases. """
        # Get all positions in the section of interest.
        positions = loader.section(self.end5, self.end3).positions
        # Remove poly(A) sequences longer than the maximum length.
        positions = filter_polya(loader.seq, positions, self.max_polya)
        if not self.include_gu:
            # Remove G and U (or T) bases.
            positions = filter_gu(loader.seq, positions)
        logger.debug(f"{self} has {positions.size} positions to load")
        return positions

    def _get_bit_vectors(self, mvec: pd.DataFrame):
        """ Convert mutation vectors into two sets of bit vectors:
        matches and mutations. """
        # Compute a boolean array of matches to the reference.
        ref_query = MATCH | INS_5
        refs = get_bits(mvec, ref_query, QEB)
        # Compute a boolean array of mutations to the reference.
        mut_query = SUB_N
        if self.include_del:
            mut_query |= DELET
        if self.include_ins:
            mut_query |= INS_3
        muts = get_bits(mvec, mut_query, QEB)
        if not mvec.shape == refs.shape == muts.shape:
            raise ValueError(f"Dimensions of mutation vectors ({mvec.shape}) "
                             f"and bit vectors ({refs.shape}, {muts.shape}) "
                             f"differed")
        logger.debug(f"Computed bit vectors with shape {mvec.shape} for {self}")
        return refs, muts

    def _filter_info(self, refs: pd.DataFrame, muts: pd.DataFrame):
        """ Remove bit vectors with too few informative positions. """
        if self.min_info > 0.0:
            # Count the information (matches + mutations) in each read.
            n_info_per_read = refs.sum(axis=1) + muts.sum(axis=1)
            f_info_per_read = n_info_per_read / self.positions.size
            # Remove in-place any reads with too few informative bits.
            too_few_info = muts.index[f_info_per_read < self.min_info]
            refs.drop(index=too_few_info, inplace=True)
            muts.drop(index=too_few_info, inplace=True)
        else:
            # No read can have too few informative bits.
            too_few_info = pd.Index([])
        # Return the number of reads with too few informative bits.
        removed = too_few_info.size
        logger.debug(f"Removed {removed} reads with an information content "
                     f"< {self.min_info} from {self}")
        return removed

    def _filter_mut_dist(self, refs: pd.DataFrame, muts: pd.DataFrame):
        """ Remove bit vectors with a pair of mutations that are too
        close to each other. """
        if self.min_mut_dist > 1:
            # Flag reads with at least one pair of mutations that are
            # too close. Initially, flag no reads (set all to False).
            flag_muts_too_close = pd.Series(np.zeros(muts.shape[0], dtype=bool),
                                            index=muts.index)
            # Loop through every position in the bit vectors.
            for pos5 in muts.columns:
                # Find the 3'-most position that must be checked. Since
                # Pandas indexing includes the last index, subtract 1 so
                # that min_mut_dist positions are in the range.
                pos3 = pos5 + self.min_mut_dist - 1
                # Get all positions that remain (e.g. after dropping Gs
                # and Us) between pos5 and pos3 (including both ends).
                # For example, if positions 1, 2, 4, 5, and 7 remain and
                # min_mut_dist = 4, then the positions checked in each
                # iteration are [1, 2, 4], [2, 4, 5], [4, 5], [4, 5, 7].
                muts_pos53 = muts.loc[:, pos5: pos3]
                # Then, sum over axis 1 to count the mutations in each
                # read between positions pos5 and pos3, inclusive. If a
                # read has >1 mutation in this range, then the mutations
                # are too close. Set all such reads' flags to True.
                flag_muts_too_close |= muts_pos53.sum(axis=1) > 1
            # Get the names of all the reads with mutations too close.
            muts_too_close = muts.index[flag_muts_too_close]
            # Remove the reads in-place.
            refs.drop(index=muts_too_close, inplace=True)
            muts.drop(index=muts_too_close, inplace=True)
        else:
            # If min_mut_dist <= 1, then no mutations can be too close
            # because even two adjacent mutations have a distance of 1.
            muts_too_close = pd.Index([])
        # Return the number of reads with mutations that were too close.
        removed = muts_too_close.size
        logger.debug(f"Removed {removed} reads with mutations closer than "
                     f"{self.min_mut_dist} nt apart from {self}")
        return removed

    def _filter_max_muts(self, refs: pd.DataFrame, muts: pd.DataFrame):
        """ Remove bit vectors with too many mutations. """
        # Count the matches and mutations in each read.
        # n_refs_per_read = refs.sum(axis=1)
        n_muts_per_read = muts.sum(axis=1)
        # Calculate the maximum permitted number of mutations in a read.
        if self.max_muts_per_read <= 0:
            mean = np.nanmean(n_muts_per_read)
            if np.isnan(mean):
                logger.warning(f"Mean number of mutations is undefined")
                mean = 0.0
            std = np.nanstd(n_muts_per_read)
            if np.isnan(std):
                logger.warning(f"Standard deviation of mutations is undefined")
                std = 0.0
            # Permit up to 3 standard deviations above the mean.
            self.max_muts_per_read = int(mean + 3.0 * std)
        # Remove in-place any reads with too many mutations.
        too_many_muts = refs.index[n_muts_per_read > self.max_muts_per_read]
        refs.drop(index=too_many_muts, inplace=True)
        muts.drop(index=too_many_muts, inplace=True)
        # Return how many reads have too many mutations.
        removed = too_many_muts.size
        logger.debug(f"Removed {removed} reads with "
                     f"> {self.max_muts_per_read} mutations from {self}")
        return removed

    def _filter_coverage(self, refs: pd.DataFrame, muts: pd.DataFrame):
        """ Remove positions with low coverage. """
        if self.min_cov > 0:
            # Calculate the coverage at each position.
            coverage = refs.sum(axis=0) + muts.sum(axis=0)
            # Remove in-place any positions with low coverage.
            too_low_coverage = refs.columns[coverage < self.min_cov]
            refs.drop(columns=too_low_coverage, inplace=True)
            muts.drop(columns=too_low_coverage, inplace=True)
            self.positions = muts.columns.values
        else:
            # No positions can have low coverage.
            too_low_coverage = pd.Index([])
        # Return how many positions have low coverage.
        removed = too_low_coverage.size
        logger.debug(f"Removed {removed} positions with coverage "
                     f"< {self.min_cov} from {self}")
        return removed

    def _filter_signal(self, refs: pd.DataFrame, muts: pd.DataFrame):
        """ Remove positions with low signal. """
        if self.min_sig > 0:
            # Calculate the signal (fraction mutated) for each position.
            n_refs_per_pos = refs.sum(axis=0)
            n_muts_per_pos = muts.sum(axis=0)
            signal = n_muts_per_pos / (n_refs_per_pos + n_muts_per_pos)
            # Remove in-place any positions with low signal.
            too_low_signal = refs.columns[signal < self.min_sig]
            refs.drop(columns=too_low_signal, inplace=True)
            muts.drop(columns=too_low_signal, inplace=True)
            self.positions = muts.columns.values
        else:
            # No positions can have low signal.
            too_low_signal = pd.Index([])
            # Return how many positions have low signal.
        removed = too_low_signal.size
        logger.debug(f"Removed {removed} positions with signal "
                     f"< {self.min_sig} from {self}")
        return removed

    def bv_report_path(self, out_dir: Path):
        """ Path to the bit vector report. """
        return path.build(path.ModSeg,
                          path.SampSeg,
                          path.RefSeg,
                          path.SectSeg,
                          path.VecRepSeg,
                          top=out_dir,
                          module=path.MOD_CLUST,
                          sample=self.sample,
                          ref=self.ref,
                          end5=self.end5,
                          end3=self.end3,
                          ext=path.JSON_EXT)

    @property
    def bv_report_data(self):
        """ Return a dictionary containing bit vector report data. """
        polya = "A" * (self.max_polya + 1)
        return {
            "Sample": self.sample,
            "Reference": self.ref,
            "5' coord": self.end5,
            "3' coord": self.end3,
            "Sequence": self.seq.decode(),
            "Positions total": len(self.seq),
            f"Positions cut, G/U/{polya}": len(self.seq) - self.n_pos_try,
            f"Positions cut, coverage < {self.min_cov}": self.n_low_cov,
            f"Positions cut, signal < {self.min_sig}": self.n_low_signal,
            "Positions used": self.n_pos_use,
            "Reads total": self.n_read_try,
            f"Reads cut, info < {self.min_info}": self.n_low_info,
            f"Reads cut, muts < {self.min_mut_dist} nt apart": self.n_near_muts,
            f"Reads cut, > {self.max_muts_per_read} muts": self.n_many_muts,
            "Reads used": self.n_read_use,
            "Unique reads used": self.n_bvec_use,
        }

    def publish_preprocessing_report(self, out_dir: Path):
        """ Publish the bit vector report in a JSON file. """
        report_path = self.bv_report_path(out_dir)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        with open(report_path, "w") as f:
            json.dump(self.bv_report_data, f)
        logger.info(f"Wrote preprocessing report of {self} to {report_path}")

    @cached_property
    def description(self):
        return (f"BitVectors for sample '{self.sample}', ref '{self.ref}', "
                f"section {self.end5}-{self.end3}")

    def __str__(self):
        return self.description
