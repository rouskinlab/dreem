from functools import cached_property
from logging import getLogger
from pathlib import Path

import numpy as np
import pandas as pd

from .report import MaskReport
from ..relate.load import RelVecLoader
from ..relate.report import RelateReport
from dreem.core.bit import BitCaller, BitCounter, BitVectorSet
from ..core.sect import seq_pos_to_index, Section

logger = getLogger(__name__)


def load_read_names_batch(batch_file: Path) -> pd.Series:
    """ Load the names of the reads in one batch from a file. """
    read_data = pd.read_csv(batch_file).squeeze()
    logger.debug(f"Loaded {read_data.size} read names from {batch_file}:\n"
                 f"{read_data}")
    return read_data


class BitVecLoader(object):
    """ Load batches of filtered bit vectors. Wrapper around CallReport
    that exposes only the attributes of the report that are required for
    loading batches of filtered bit vectors. """

    def __init__(self, report: MaskReport):
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
        return Section(self.ref, self.relvec.seq,
                       end5=self.end5, end3=self.end3, name=self.sect)

    @property
    def seq(self):
        return self.section.seq

    @property
    def positions(self):
        return self.section.positions

    @property
    def pos_init(self):
        return self.positions

    @property
    def n_pos_init(self):
        return self.section.length

    @cached_property
    def index_kept(self):
        return seq_pos_to_index(self.section.seq,
                                self.pos_kept,
                                start=self.section.end5)

    @cached_property
    def relvec(self):
        return RelVecLoader.open(RelateReport.build_path(self.out_dir,
                                                         sample=self.sample,
                                                         ref=self.ref))

    def load_relvec_batch(self, batch: int) -> pd.DataFrame:
        """ Load and filter relation vectors in one batch. """
        # Load the names of the selected reads.
        reads = load_read_names_batch(self.get_batch_path(batch))
        # Load the selected positions and reads of the relation vectors.
        return self.relvec.load_batch(batch, self.pos_kept).loc[reads]

    def iter_relvec_batches(self):
        return map(self.load_relvec_batch, range(self.n_batches))

    def get_read_names(self):
        try:
            # Concatenate all indexes, which have the read names.
            return np.hstack(bat.index.values
                             for bat in self.iter_relvec_batches())
        except ValueError:
            # If there are no batches, return an empty array.
            return np.array([], dtype=str)

    def get_bit_caller(self):
        return BitCaller(self.count_refs, self.count_muts)

    def get_bit_counter(self):
        """ Return a BitCounter object of all filtered bit vectors. """
        return BitCounter(self.get_bit_caller(), self.iter_relvec_batches())

    def get_bit_vectors(self):
        """ Return a BitVectorSet object of all filtered bit vectors. """
        return BitVectorSet(self.get_bit_caller(), self.iter_relvec_batches())

    @classmethod
    def open(cls, report_file: Path):
        return cls(MaskReport.open(report_file))

    def __str__(self):
        return f"Bit Vectors of {self.sample}@{self.section.ref_sect}"
