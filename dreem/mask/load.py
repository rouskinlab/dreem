from functools import cached_property
from logging import getLogger
from pathlib import Path

import numpy as np
import pandas as pd

from .report import MaskReport
from ..core.bitc import BitCaller
from ..core.load import SectBatchChainLoader, no_kwargs
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

    @cached_property
    def bit_caller(self):
        """ Get the BitCaller associated with the mask. """
        return BitCaller(self.count_refs, self.count_muts)

    @no_kwargs
    def load_data_personal(self, batch_file: Path):
        # This method accepts no keyword arguments.
        # Load the names of the reads in the batch from the file.
        read_names = pd.read_csv(batch_file).squeeze()
        logger.debug(f"Loaded {read_names.size} read names from {batch_file}:\n"
                     f"{read_names}")
        return read_names

    def iter_batches_personal(self):
        yield from super().iter_batches_personal()

    def iter_read_batches(self):
        """ Yield the read names that were kept in every batch. """
        for batch in self.iter_batches_personal():
            yield batch.values

    def get_read_names(self):
        """ Return an array naming all reads that were kept. """
        try:
            # Concatenate all indexes, which have the read names.
            return np.hstack(self.iter_read_batches())
        except ValueError:
            # If there are no batches, return an empty array.
            return np.array([], dtype=str)

    def process_batch(self, imported_batch: pd.DataFrame,
                      private_batch: pd.Series, *,
                      bit_caller: BitCaller | None = None,
                      merge: bool = False, invert: bool = False):
        if bit_caller is None:
            # If no BitCaller was given, then use the mask's bit caller.
            bit_caller = self.bit_caller
        else:
            # Otherwise, intersect the given bit caller with the mask's
            # bit caller so that a type of relationships is counted only
            # if both bit callers accept it.
            bit_caller = BitCaller.inter(bit_caller, self.bit_caller,
                                         merge=merge, invert=invert)
        # Select only the relation vectors that were kept after masking
        # using imported_batch.loc[private_batch], then call the bits
        # using the BitCaller.
        return bit_caller.call(imported_batch.loc[private_batch])

    def iter_batches_processed(self, *, bit_caller: BitCaller | None = None,
                               merge: bool = False, invert: bool = False):
        yield from super().iter_batches_processed(positions=self.pos_kept,
                                                  bit_caller=bit_caller,
                                                  merge=merge, invert=invert)
