from logging import getLogger
from typing import Iterable

import numpy as np
import pandas as pd

from .load import load_read_names_batch
from .report import CallReport
from .write import write_batch
from ..core.bit import BitCaller, BitCounter, BitVectorSet
from ..core.sect import mask_gu, mask_polya, mask_pos, Section
from ..relate.load import RelVecLoader

logger = getLogger(__name__)


class BitFilter(object):
    """ Filter bit vectors. """

    def __init__(self, /,
                 loader: RelVecLoader,
                 bit_caller: BitCaller,
                 section: Section | None = None, *,
                 exclude_polya: int = 0,
                 exclude_gu: bool = False,
                 exclude_pos: Iterable[int] = (),
                 min_mut_gap: int = 0,
                 min_finfo_read: float = 0.,
                 max_fmut_read: float = 1.,
                 min_ninfo_pos: int = 0,
                 max_fmut_pos: float = 1.):
        """
        Parameters
        ----------
        loader: RelVecLoader
            Relation vector loader
        bit_caller: BitCaller
            Bit caller
        section: Section | None = None
            Section of the reference sequence to use. If omitted, use
            the entire sequence.
        exclude_polya: int = 0
            Exclude stretches of consecutive A bases at least this long.
            If 0, exclude no bases. Must be ≥ 0.
        exclude_gu: bool = False
            Whether to exclude G and U bases.
        exclude_pos: Iterable[int] = ()
            Additional, arbitrary positions to exclude.
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
            bases. Must be ≥ 0.
        max_fmut_pos: float = 1.0
            Filter out positions with more than this fraction of mutated
            reads. Must be in [0, 1].
        """
        self.out_dir = loader.out_dir
        self.sample = loader.sample
        self.ref = loader.ref
        self.section = loader.get_section() if section is None else section
        self.bit_caller = bit_caller
        self.loader_str = str(loader)
        logger.info(f"Began {self}")
        # Define parameters for excluding positions from the section.
        self.exclude_polya = exclude_polya
        self.exclude_gu = exclude_gu
        self.exclude_pos = np.array(exclude_pos, dtype=int)
        # Exclude positions that meet those parameters.
        self.pos_mask = np.zeros_like(self.section.positions, dtype=bool)
        self.pos_polya = self._mask_pos(self.section.positions,
                                        mask_polya(self.section.seq,
                                                   self.exclude_polya))
        logger.debug(f"Excluded {len(self.pos_polya)} positions with poly(A) "
                     f"sequences of ≥ {self.exclude_polya} nt from {self}")
        self.pos_gu = self._mask_pos(self.section.positions,
                                     mask_gu(self.section.seq,
                                             self.exclude_gu))
        logger.debug(f"Excluded {len(self.pos_gu)} positions with G or U "
                     f"bases from {self}")
        self.pos_user = self._mask_pos(self.section.positions,
                                       mask_pos(self.section.positions,
                                                self.exclude_pos))
        logger.debug(f"Excluded {len(self.pos_user)} positions that were "
                     f"pre-specified for exclusion from {self}")
        # Determine which positions remain after excluding the above.
        self.pos_load = self.section.positions[np.logical_not(self.pos_mask)]
        # Set the parameters for filtering reads.
        self.min_mut_gap = min_mut_gap
        self.min_finfo_read = min_finfo_read
        self.max_fmut_read = max_fmut_read
        # Load every batch of mutation vectors, count the informative
        # and mutated bits in every vector and position, and drop reads
        # that do not pass the filters.
        self.bit_counter = BitCounter(bit_caller,
                                      loader.iter_batches(self.pos_load),
                                      filters=[self._mask_min_finfo_read,
                                               self._mask_max_fmut_read,
                                               self._mask_min_mut_gap])
        # Write batches of read names and record their checksums.
        self.checksums = [write_batch(read_names, self._batch_path(batch))
                          for batch, read_names
                          in enumerate(self.bit_counter.read_batches)]
        # Set the parameters for filtering positions.
        self.min_ninfo_pos = min_ninfo_pos
        self.max_fmut_pos = max_fmut_pos
        # Remove positions that do not pass the filters.
        self.pos_mask = np.zeros_like(self.pos_load, dtype=bool)
        if not self.min_ninfo_pos >= 0:
            raise ValueError(
                f"min_ninfo_pos must be ≥ 0, but got {self.min_ninfo_pos}")
        self.pos_min_ninfo = self._mask_pos(self.pos_load,
                                            (self.bit_counter.ninfo_per_pos
                                             < self.min_ninfo_pos))
        logger.debug(f"Dropped {len(self.pos_min_ninfo)} positions with "
                     f"< {self.min_ninfo_pos} informative reads from {self}")
        if not 0. <= self.max_fmut_pos <= 1.:
            raise ValueError(f"max_fmut_pos must be in [0, 1], but got "
                             f"{self.max_fmut_pos}")
        self.pos_max_fmut = self._mask_pos(self.pos_load,
                                           (self.bit_counter.fmuts_per_pos
                                            > self.max_fmut_pos))
        logger.debug(f"Dropped {len(self.pos_max_fmut)} positions with mutated "
                     f"fractions > {self.max_fmut_pos} from {self}")
        self.pos_kept = self.pos_load[np.logical_not(self.pos_mask)]
        logger.info(f"Ended filtering positions and reads from {loader} "
                    f"over {section} with {bit_caller}: "
                    f"kept {self.pos_kept.size} positions and "
                    f"{self.n_reads_kept} reads")

    @property
    def n_reads_init(self):
        return self.bit_counter.nvec_given

    @property
    def n_reads_min_finfo(self):
        return self.bit_counter.nvec_filtered[0]

    @property
    def n_reads_max_fmut(self):
        return self.bit_counter.nvec_filtered[1]

    @property
    def n_reads_min_gap(self):
        return self.bit_counter.nvec_filtered[2]

    @property
    def n_reads_kept(self):
        return self.bit_counter.nvec

    @property
    def n_batches(self):
        return len(self.checksums)

    def _mask_pos(self, all_pos: np.ndarray, exclude: np.ndarray) -> np.ndarray:
        """ Exclude positions before counting bits. """
        # Determine which positions to exclude are new (have not been
        # excluded already).
        new_excluded = np.logical_and(exclude, np.logical_not(self.pos_mask))
        # Update the record of all positions that have been excluded.
        self.pos_mask = np.logical_or(new_excluded, self.pos_mask)
        # Return the newly excluded positions.
        return all_pos[new_excluded]

    def _mask_min_finfo_read(self, bvec: BitVectorSet) -> np.ndarray:
        """ Mask reads with too few informative bits. """
        if not 0. <= self.min_finfo_read <= 1.:
            raise ValueError(f"min_finfo_read must be in [0, 1], but got "
                             f"{self.min_finfo_read}")
        if self.min_finfo_read == 0.:
            # Nothing to mask.
            return np.zeros_like(bvec.reads, dtype=bool)
        return bvec.finfo_per_vec.values < self.min_finfo_read

    def _mask_max_fmut_read(self, bvec: BitVectorSet) -> np.ndarray:
        """ Mask reads with too many mutations. """
        if not 0. <= self.max_fmut_read <= 1.:
            raise ValueError(f"max_fmut_read must be in [0, 1], but got "
                             f"{self.max_fmut_read}")
        if self.max_fmut_read == 1.:
            # Nothing to mask.
            return np.zeros_like(bvec.reads, dtype=bool)
        return bvec.fmuts_per_vec.values > self.max_fmut_read

    def _mask_min_mut_gap(self, bvec: BitVectorSet) -> np.ndarray:
        """ Mask reads with mutations that are too close. """
        if not self.min_mut_gap >= 0:
            raise ValueError(
                f"min_mut_gap must be ≥ 0, but got {self.min_mut_gap}")
        # Initially, mask no reads (set all to False).
        mask = np.zeros_like(bvec.reads, dtype=bool)
        if self.min_mut_gap == 0:
            # Nothing to mask.
            return mask
        # Make a mask for the current window of positions.
        window = pd.Series(np.zeros_like(bvec.positions, dtype=bool),
                           index=bvec.positions)
        # Loop through every position in the bit vectors.
        for pos5 in bvec.positions:
            # Find the 3'-most position that must be checked.
            pos3 = pos5 + self.min_mut_gap
            # Mask all positions between pos5 and pos3, inclusive.
            # For example, if positions 1, 2, 4, 5, and 7 remain and
            # min_gap = 3, then the positions checked in each
            # iteration are [1, 2, 4], [2, 4, 5], [4, 5], [4, 5, 7].
            window.loc[pos5: pos3] = True
            # Sum over axis 1 to count the mutations in each read
            # between positions pos5 and pos3, inclusive. If a read
            # has > 1 mutation in this range, then it has mutations
            # that are too close. Set all such reads' flags to True.
            mask |= bvec.muts.loc[:, window.values].sum(axis=1) > 1
            # Erase the mask over the current window.
            window.loc[pos5: pos3] = False
            if pos3 >= bvec.positions[-1]:
                # The end of the sequence has been reached. Stop.
                break
        return mask

    def _batch_path(self, batch: int):
        """ Path to a file of read names in a batch. """
        return CallReport.build_batch_path(self.out_dir, batch,
                                           sample=self.sample,
                                           ref=self.ref,
                                           sect=self.section.name)

    def _load_batch(self, batch: int):
        """ Load the names of the reads in one batch from a file. """
        return load_read_names_batch(self._batch_path(batch))

    def create_report(self):
        """ Create a FilterReport from a BitFilter object. """
        return CallReport(
            out_dir=self.out_dir,
            sample=self.sample,
            ref=self.ref,
            sect=self.section.name,
            end5=self.section.end5,
            end3=self.section.end3,
            checksums=self.checksums,
            n_batches=self.n_batches,
            count_refs=self.bit_caller.ref_caller,
            count_muts=self.bit_caller.mut_caller,
            exclude_gu=self.exclude_gu,
            exclude_polya=self.exclude_polya,
            exclude_pos=self.exclude_pos,
            min_ninfo_pos=self.min_ninfo_pos,
            max_fmut_pos=self.max_fmut_pos,
            n_pos_init=self.section.length,
            n_pos_gu=len(self.pos_gu),
            n_pos_polya=len(self.pos_polya),
            n_pos_user=len(self.pos_user),
            n_pos_min_ninfo=len(self.pos_min_ninfo),
            n_pos_max_fmut=len(self.pos_max_fmut),
            n_pos_kept=len(self.pos_kept),
            pos_gu=self.pos_gu,
            pos_polya=self.pos_polya,
            pos_user=self.pos_user,
            pos_min_ninfo=self.pos_min_ninfo,
            pos_max_fmut=self.pos_max_fmut,
            pos_kept=self.pos_kept,
            min_finfo_read=self.min_finfo_read,
            max_fmut_read=self.max_fmut_read,
            min_mut_gap=self.min_mut_gap,
            n_reads_init=self.n_reads_init,
            n_reads_min_finfo=self.n_reads_min_finfo,
            n_reads_max_fmut=self.n_reads_max_fmut,
            n_reads_min_gap=self.n_reads_min_gap,
            n_reads_kept=self.n_reads_kept,
        )

    def __str__(self):
        return f"Filter {self.loader_str} {self.section} with {self.bit_caller}"


def filter_sect(loader: RelVecLoader,
                section: Section,
                count_del: bool,
                count_ins: bool,
                discount: Iterable[str], *,
                rerun: bool,
                **kwargs):
    """ Filter a section of a set of bit vectors. """
    # Check if the report file already exists.
    report_file = CallReport.build_path(loader.out_dir,
                                        sample=loader.sample,
                                        ref=loader.ref,
                                        sect=section.name)
    if rerun or not report_file.is_file():
        # Create the bit caller.
        bit_caller = BitCaller.from_counts(count_del, count_ins, discount)
        # Create and apply the filter.
        bit_filter = BitFilter(loader, bit_caller, section, **kwargs)
        # Output the results of filtering.
        report = bit_filter.create_report()
        report.save()
    else:
        logger.warning(f"File exists: {report_file}")
    return report_file
