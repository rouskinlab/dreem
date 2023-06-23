"""
DREEM -- Mask Module
"""

from logging import getLogger
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .report import MaskReport
from ..core import path
from ..core.bitc import BitCaller
from ..core.bitv import BitBatch, BitCounter
from ..core.files import digest_file
from ..core.sect import mask_gu, mask_polya, mask_pos, Section
from ..relate.load import RelateLoader

logger = getLogger(__name__)

MASK_FINFO = "finfo"
MASK_FMUT = "fmut"
MASK_GAP = "gap"


def write_batch(read_names: Iterable[str],
                out_dir: Path, sample: str, ref: str, sect: str, batch: int):
    """ Write the names of the reads in one batch to a file. """
    # Determine the path to the batch file.
    batch_file = MaskReport.build_batch_path(out_dir, batch,
                                             sample=sample, ref=ref, sect=sect,
                                             ext=path.CSVZIP_EXT)
    # Write the read names to the batch file.
    read_data = pd.Series(read_names)
    read_data.to_csv(batch_file, index=False)
    logger.debug(f"Wrote {read_data.size} read names to {batch_file}:\n"
                 f"{read_data}")
    # Compute the checksum of the batch file.
    return digest_file(batch_file)


class BitMasker(object):
    """ Mask bit vectors. """

    def __init__(self, /,
                 loader: RelateLoader,
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
        loader: RelateLoader
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
        # Set the general parameters.
        self.out_dir = loader.out_dir
        self.sample = loader.sample
        self.ref = loader.ref
        self.section = (Section(loader.ref, loader.seq,
                                end5=1, end3=len(loader.seq)) if section is None
                        else section)
        self.bit_caller = bit_caller
        self.loader_str = str(loader)
        logger.info(f"Began {self}")
        # Set the parameters for excluding positions from the section.
        self.exclude_polya = exclude_polya
        self.exclude_gu = exclude_gu
        self.exclude_pos = np.array(exclude_pos, dtype=int)
        # Set the parameters for filtering reads.
        self.min_mut_gap = min_mut_gap
        self.min_finfo_read = min_finfo_read
        self.max_fmut_read = max_fmut_read
        # Set the parameters for filtering positions.
        if not min_ninfo_pos >= 0:
            raise ValueError(
                f"min_ninfo_pos must be ≥ 0, but got {min_ninfo_pos}")
        self.min_ninfo_pos = min_ninfo_pos
        if not 0. <= max_fmut_pos < 1.:
            raise ValueError(
                f"max_fmut_pos must be in [0, 1), but got {max_fmut_pos}")
        self.max_fmut_pos = max_fmut_pos
        # Exclude positions based on the  parameters.
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
        # Load only those positions from each batch of mutation vectors,
        # count the informative and mutated bits in every vector and
        # position, and drop reads that do not pass the filters.
        rel_batches = loader.iter_batches_processed(positions=self.pos_load)
        masks = {MASK_FINFO: self._mask_min_finfo_read,
                 MASK_FMUT: self._mask_max_fmut_read,
                 MASK_GAP: self._mask_min_mut_gap}
        self.counter = BitCounter(bit_caller.iter(rel_batches, masks))
        # Write batches of read names and record their checksums.
        self.checksums = [write_batch(read_names, self.out_dir, self.sample,
                                      self.ref, self.section.name, batch)
                          for batch, read_names
                          in enumerate(self.counter.read_batches)]
        # Check if any reads were counted.
        if self.counter.nreads > 0:
            # Remove positions that do not pass the filters.
            self.pos_mask = np.zeros_like(self.pos_load, dtype=bool)
            self.pos_min_ninfo = self._mask_pos(self.pos_load,
                                                (self.counter.nall_per_pos
                                                 < self.min_ninfo_pos))
            logger.debug(f"Dropped {len(self.pos_min_ninfo)} positions with < "
                         f"{self.min_ninfo_pos} informative reads from {self}")
            self.pos_max_fmut = self._mask_pos(self.pos_load,
                                               (self.counter.fyes_per_pos
                                                > self.max_fmut_pos))
            logger.debug(f"Dropped {len(self.pos_max_fmut)} positions with "
                         f"mutated fractions > {self.max_fmut_pos} from {self}")
        else:
            # If no reads were counted, then set dummy values for the
            # remaining parameters so that a report can be written.
            logger.critical(f"No reads passed through {self}")
            self.pos_mask = np.ones_like(self.pos_load, dtype=bool)
            self.pos_min_ninfo = self.pos_load
            self.pos_max_fmut = np.array([], dtype=int)
        self.pos_kept = self.pos_load[np.logical_not(self.pos_mask)]
        # Check if any positions passed the filters.
        if self.pos_kept.size == 0:
            logger.critical(f"No positions passed through {self}")
        logger.info(f"Ended filtering positions and reads from {loader} "
                    f"over {section} with {bit_caller}: "
                    f"kept {self.pos_kept.size} positions and "
                    f"{self.n_reads_kept} reads")

    @property
    def n_reads_init(self):
        return self.counter.nreads_given

    @property
    def n_reads_min_finfo(self):
        return self.counter.nmasked[MASK_FINFO]

    @property
    def n_reads_max_fmut(self):
        return self.counter.nmasked[MASK_FMUT]

    @property
    def n_reads_min_gap(self):
        return self.counter.nmasked[MASK_GAP]

    @property
    def n_reads_kept(self):
        return self.counter.nreads

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

    def _mask_min_finfo_read(self, batch: BitBatch) -> np.ndarray:
        """ Mask reads with too few informative bits. """
        if not 0. <= self.min_finfo_read <= 1.:
            raise ValueError(f"min_finfo_read must be in [0, 1], but got "
                             f"{self.min_finfo_read}")
        if self.min_finfo_read == 0.:
            # Nothing to mask.
            return np.zeros_like(batch.reads, dtype=bool)
        return batch.fall_per_read.values < self.min_finfo_read

    def _mask_max_fmut_read(self, batch: BitBatch) -> np.ndarray:
        """ Mask reads with too many mutations. """
        if not 0. <= self.max_fmut_read <= 1.:
            raise ValueError(f"max_fmut_read must be in [0, 1], but got "
                             f"{self.max_fmut_read}")
        if self.max_fmut_read == 1.:
            # Nothing to mask.
            return np.zeros_like(batch.reads, dtype=bool)
        return batch.fyes_per_read.values > self.max_fmut_read

    def _mask_min_mut_gap(self, batch: BitBatch) -> np.ndarray:
        """ Mask reads with mutations that are too close. """
        if not self.min_mut_gap >= 0:
            raise ValueError(
                f"min_mut_gap must be ≥ 0, but got {self.min_mut_gap}")
        # Initially, mask no reads (set all to False).
        mask = np.zeros_like(batch.reads, dtype=bool)
        if self.min_mut_gap == 0:
            # Nothing to mask.
            return mask
        # Make a mask for the current window of positions.
        window = pd.Series(np.zeros_like(batch.pos, dtype=bool),
                           index=batch.pos)
        # Loop through every position in the bit vectors.
        for pos5 in batch.pos:
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
            mask |= batch.yes.loc[:, window.values].sum(axis=1) > 1
            # Erase the mask over the current window.
            window.loc[pos5: pos3] = False
            if pos3 >= batch.pos[-1]:
                # The end of the sequence has been reached. Stop.
                break
        return mask

    def create_report(self):
        """ Create a FilterReport from a BitFilter object. """
        return MaskReport(
            out_dir=self.out_dir,
            sample=self.sample,
            ref=self.ref,
            sect=self.section.name,
            end5=self.section.end5,
            end3=self.section.end3,
            checksums=self.checksums,
            n_batches=self.n_batches,
            count_refs=self.bit_caller.nos_call,
            count_muts=self.bit_caller.yes_call,
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
        return f"Mask {self.loader_str} {self.section} with {self.bit_caller}"


def mask_section(loader: RelateLoader,
                 section: Section,
                 count_del: bool,
                 count_ins: bool,
                 discount: Iterable[str], *,
                 rerun: bool,
                 **kwargs):
    """ Filter a section of a set of bit vectors. """
    # Check if the report file already exists.
    report_file = MaskReport.build_path(loader.out_dir,
                                        sample=loader.sample,
                                        ref=loader.ref,
                                        sect=section.name)
    if rerun or not report_file.is_file():
        # Create the bit caller.
        bit_caller = BitCaller.from_counts(count_del, count_ins, discount)
        # Create and apply the mask.
        bit_masker = BitMasker(loader, bit_caller, section, **kwargs)
        # Output the results of filtering.
        report = bit_masker.create_report()
        report.save()
    else:
        logger.warning(f"File exists: {report_file}")
    return report_file
