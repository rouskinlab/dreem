from functools import cached_property
from logging import getLogger
from pathlib import Path

import numpy as np
import pandas as pd

from .report import MaskReport
from ..core.bitc import BitCaller
from ..core.bitv import BitCounter, BitMonolith
from ..core.load import SectBatchChainLoader
from ..core.sect import seq_pos_to_index
from ..relate.load import RelateLoader

logger = getLogger(__name__)


class MaskLoader(SectBatchChainLoader):
    """ Load batches of masked relation vectors. """

    def __init__(self, report: MaskReport):
        super().__init__(report)
        self.pos_kept = report.pos_kept
        self.min_mut_gap = report.min_mut_gap
        self.count_refs = report.count_refs
        self.count_muts = report.count_muts

    @classmethod
    def get_report_type(cls):
        return MaskReport

    @classmethod
    def get_upstream_type(cls):
        return RelateLoader

    def _get_refseq(self):
        return self._upload.seq

    @cached_property
    def index_kept(self):
        """ Indexes of the positions that were kept. """
        return seq_pos_to_index(self.section.seq,
                                self.pos_kept,
                                start=self.section.end5)

    def load_data(self, batch_file: Path):
        # Load the names of the reads in the batch from the file.
        read_data = pd.read_csv(batch_file).squeeze()
        logger.debug(f"Loaded {read_data.size} read names from {batch_file}:\n"
                     f"{read_data}")
        return read_data

    def _iter_upstream_batches(self):
        # Load only the positions that were kept after masking.
        yield from super()._iter_upstream_batches(positions=self.pos_kept)

    def _process_batch(self, rel_batch: pd.DataFrame, read_batch: pd.Series):
        # Load only the reads that were kept after masking, and call the
        # informative and mutated bits.
        return self.bit_caller.call(rel_batch.loc[read_batch])

    def get_read_names(self):
        try:
            # Concatenate all indexes, which have the read names.
            return np.hstack(self.load_data(batch_file).values
                             for batch_file in self._rep.iter_batch_paths())
        except ValueError:
            # If there are no batches, return an empty array.
            return np.array([], dtype=str)

    @cached_property
    def bit_caller(self):
        """ Get the BitCaller associated with the mask. """
        return BitCaller(self.count_refs, self.count_muts)

    def get_bit_counts(self):
        """ Count the bits and return them as a BitCounter. """
        return BitCounter(self.iter_batches())

    def get_bit_monolith(self):
        """ Return a BitMonolith object of all filtered bit vectors. """
        return BitMonolith(self.iter_batches())
