from functools import cached_property
from logging import getLogger
from pathlib import Path

import numpy as np
import pandas as pd

from .report import CallReport
from ..relate.load import RelaVecLoader
from ..relate.report import RelateReport
from ..core.bit import BitCaller, BitCounter, BitVectorSet
from ..core.sect import seq_pos_to_cols, Section

logger = getLogger(__name__)


def load_read_names_batch(batch_file: Path) -> pd.Series:
    """ Load the names of the reads in one batch from a file. """
    read_data = pd.read_csv(batch_file).squeeze()
    logger.debug(f"Loaded {read_data.size} read names from {batch_file}:\n"
                 f"{read_data}")
    return read_data


class CallVecLoader(object):
    """ Load batches of filtered mutation vectors. Wrapper around
    FilterReport that exposes only the attributes of the report that are
    required for loading batches of mutation vectors. """

    def __init__(self, report: CallReport):
        self.out_dir = report.out_dir
        self.sample = report.sample
        self.ref = report.ref
        self.sect = report.sect
        self.end5 = report.end5
        self.end3 = report.end3
        self.n_batches = report.n_batches
        self.pos_kept = report.pos_kept
        self.min_mut_gap = report.min_mut_gap
        self.count_refs = report.count_refs
        self.count_muts = report.count_muts
        self.get_batch_path = report.get_batch_path

    @cached_property
    def section(self):
        return Section(self.ref, self.relavec_loader.seq,
                       end5=self.end5, end3=self.end3, name=self.sect)

    @cached_property
    def pos_init(self):
        return np.arange(self.end5, self.end3 + 1)

    @property
    def n_pos_init(self):
        return self.end3 - self.end5 + 1

    @cached_property
    def cols_kept(self):
        return seq_pos_to_cols(self.section.seq, self.pos_kept)

    @cached_property
    def relavec_loader(self):
        return RelaVecLoader.open(RelateReport.build_path(self.out_dir,
                                                          sample=self.sample,
                                                          ref=self.ref))

    def load_relavec_batch(self, batch: int) -> pd.DataFrame:
        """ Load and filter relation vectors in one batch. """
        # Load the names of the selected reads.
        reads = load_read_names_batch(self.get_batch_path(batch))
        # Load the selected positions and reads of the relation vectors.
        return self.relavec_loader.load_batch(batch, self.pos_kept).loc[reads]

    def iter_relavec_batches(self):
        return map(self.load_relavec_batch, range(self.n_batches))

    def get_read_names(self):
        try:
            # Concatenate all indexes, which have the read names.
            return np.hstack(bat.index.values
                             for bat in self.iter_relavec_batches())
        except ValueError:
            # If there are no batches, return an empty array.
            return np.array([], dtype=str)

    def get_bit_caller(self):
        return BitCaller(self.count_refs, self.count_muts)

    def get_bit_counter(self):
        """ Return a BitCounter object of all filtered bit vectors. """
        return BitCounter(self.get_bit_caller(), self.iter_relavec_batches())

    def get_bit_vectors(self):
        """ Return a BitVectorSet object of all filtered bit vectors. """
        return BitVectorSet(self.get_bit_caller(), self.iter_relavec_batches())

    @classmethod
    def open(cls, report_file: Path):
        return cls(CallReport.open(report_file))
