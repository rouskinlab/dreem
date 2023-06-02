from functools import cached_property
from logging import getLogger
from pathlib import Path

import numpy as np
import pandas as pd

from .report import MaskReport
from ..core.bitc import BitCaller
from ..core.bitv import BitCounter, BitMonolith
from ..core.load import DataLoader
from ..core.sect import seq_pos_to_index, Section
from ..relate.load import RelateLoader
from ..relate.report import RelateReport

logger = getLogger(__name__)


def load_reads_batch(batch_file: Path) -> pd.Series:
    """ Load the names of the reads in one batch from a file. """
    read_data = pd.read_csv(batch_file).squeeze()
    logger.debug(f"Loaded {read_data.size} read names from {batch_file}:\n"
                 f"{read_data}")
    return read_data


class MaskLoader(DataLoader):
    """ Load batches of masked relation vectors. """

    def __init__(self, report: MaskReport):
        super().__init__(report)
        self._sect = report.sect
        self.n_batches = report.n_batches
        self.pos_kept = report.pos_kept
        self.min_mut_gap = report.min_mut_gap
        self.count_refs = report.count_refs
        self.count_muts = report.count_muts
        self.get_batch_path = report.get_batch_path

    @classmethod
    def report_type(cls):
        return MaskReport

    @property
    def sect(self) -> str:
        return self._sect

    @cached_property
    def section(self):
        return Section(self.ref, self.rel_load.seq,
                       end5=self.end5, end3=self.end3, name=self.sect)

    @property
    def seq(self):
        return self.section.seq

    @cached_property
    def index_kept(self):
        """ Indexes of the positions that were kept. """
        return seq_pos_to_index(self.section.seq,
                                self.pos_kept,
                                start=self.section.end5)

    @cached_property
    def rel_load(self):
        return RelateLoader.open(RelateReport.build_path(self.out_dir,
                                                         sample=self.sample,
                                                         ref=self.ref))

    def load_rel_batch(self, batch: int) -> pd.DataFrame:
        """ Load and filter relation vectors in one batch. """
        # Load the names of the selected reads.
        reads = load_reads_batch(self.get_batch_path(batch))
        # Load the selected positions and reads of the relation vectors.
        return self.rel_load.load_rel_batch(batch, self.pos_kept).loc[reads]

    def iter_rel_batches(self):
        return map(self.load_rel_batch, range(self.n_batches))

    def get_read_names(self):
        try:
            # Concatenate all indexes, which have the read names.
            return np.hstack(load_reads_batch(self.get_batch_path(batch)).values
                             for batch in range(self.n_batches))
        except ValueError:
            # If there are no batches, return an empty array.
            return np.array([], dtype=str)

    @cached_property
    def bit_caller(self):
        """ Get the BitCaller associated with the mask. """
        return BitCaller(self.count_refs, self.count_muts)

    def iter_bit_batches(self):
        """ Yield every batch of filtered bit vectors. """
        return map(self.bit_caller.call, self.iter_rel_batches())

    def get_bit_counts(self):
        """ Count the bits and return them as a BitCounter. """
        return BitCounter(self.iter_bit_batches())

    def get_bit_monolith(self):
        """ Return a BitMonolith object of all filtered bit vectors. """
        return BitMonolith(self.iter_bit_batches())

    @classmethod
    def open(cls, report_file: Path):
        return cls(MaskReport.open(report_file))

    def __str__(self):
        return f"Bit Vectors of {self.sample}@{self.section.ref_sect}"
