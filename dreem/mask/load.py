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
    def get_import_type(cls):
        return RelateLoader

    def get_refseq(self):
        return self.import_loader.seq

    @cached_property
    def index_kept(self):
        """ Indexes of the positions that were kept. """
        return seq_pos_to_index(self.section.seq,
                                self.pos_kept,
                                start=self.section.end5)

    def iter_batches_import(self):
        # Modify the ancestral function to load only the positions kept.
        yield from self.import_loader.iter_batches_public(self.pos_kept)

    def _load_data_private(self, batch_file: Path):
        # Load the names of the reads in the batch from the file.
        read_names = pd.read_csv(batch_file).squeeze()
        logger.debug(f"Loaded {read_names.size} read names from {batch_file}:\n"
                     f"{read_names}")
        return read_names

    def _publish_batch(self, private_batch: pd.Series,
                       imported_batch: pd.DataFrame,
                       *bit_callers: BitCaller):
        # Load only the relation vectors that were kept after masking.
        relvecs_kept = imported_batch.loc[private_batch]
        if not bit_callers:
            # If no bit callers were given, then use the mask's own bit
            # caller to call the bits.
            bit_callers = self.bit_caller,
        # Yield a BitBatch from the relation vectors for each BitCaller.
        for bit_caller in bit_callers:
            yield bit_caller.call(relvecs_kept)

    def get_read_names(self):
        """ Return an array naming all reads that were kept. """
        try:
            # Concatenate all indexes, which have the read names.
            return np.hstack(self._load_batch_private(batch_file).values
                             for batch_file in self._report.iter_batch_paths())
        except ValueError:
            # If there are no batches, return an empty array.
            return np.array([], dtype=str)

    @cached_property
    def bit_caller(self):
        """ Get the BitCaller associated with the mask. """
        return BitCaller(self.count_refs, self.count_muts)

    def get_bit_counts(self):
        """ Return a BitCounter of all masked bit vectors. """
        return BitCounter(self.iter_batches_public())

    def get_bit_monolith(self):
        """ Return a BitMonolith of all masked bit vectors. """
        return BitMonolith(self.iter_batches_public())
