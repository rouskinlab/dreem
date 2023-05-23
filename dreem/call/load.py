from functools import cached_property
from logging import getLogger
from pathlib import Path

import numpy as np
import pandas as pd

from .report import CallReport
from ..relate.load import RelaVecLoader
from ..relate.report import RelateReport
from ..core.bit import BitCaller, BitCounter, BitVectorSet
from ..core.sect import seq_pos_to_cols

logger = getLogger(__name__)


def load_batch(batch_file: Path) -> pd.Series:
    """ Load the names of the reads in one batch from a file. """
    read_data = pd.read_csv(batch_file).squeeze()
    logger.debug(f"Loaded {read_data.size} read names from {batch_file}:\n"
                 f"{read_data}")
    return read_data


class CallLoader(object):
    """ Load batches of filtered mutation vectors. Wrapper around
    FilterReport that exposes only the attributes of the report that are
    required for loading batches of mutation vectors. """

    def __init__(self, report: CallReport):
        self._rep = report

    @property
    def out_dir(self):
        return self._rep.out_dir

    @property
    def sample(self):
        return self._rep.sample

    @property
    def ref(self):
        return self._rep.ref

    @property
    def sect(self):
        return self._rep.sect

    @property
    def seq(self):
        return self._rep.seq

    @property
    def n_batches(self):
        return self._rep.n_batches

    @property
    def end5(self):
        return self._rep.end5

    @property
    def end3(self):
        return self._rep.end3

    @cached_property
    def pos_init(self):
        return np.arange(self.end5, self.end3 + 1)

    @property
    def n_pos_init(self):
        return self.end3 - self.end5 + 1

    @property
    def pos_kept(self):
        return np.array(self._rep.pos_kept)

    @cached_property
    def cols_kept(self):
        return seq_pos_to_cols(self.seq, self.pos_kept)

    @property
    def n_pos_kept(self):
        return self._rep.n_pos_kept

    @property
    def min_mut_gap(self):
        return self._rep.min_mut_gap

    @property
    def mvec_report_file(self):
        return RelateReport.build_path(self.out_dir,
                                       sample=self.sample,
                                       ref=self.ref)

    @cached_property
    def mvec_loader(self):
        return RelaVecLoader.open(self.mvec_report_file)

    def load_batch(self, batch: int) -> pd.DataFrame:
        """ Load and filter mutation vectors in one batch. """
        # Load the names of the selected reads.
        reads = load_batch(self._rep.get_batch_path(batch))
        # Load the selected positions and reads of the mutation vectors.
        return self.mvec_loader.load_batch(batch, self._rep.pos_kept).loc[reads]

    def iter_batches(self):
        return map(self.load_batch, range(self.n_batches))

    def get_read_names(self):
        try:
            # Concatenate all indexes, which have the read names.
            return np.hstack(bat.index.values for bat in self.iter_batches())
        except ValueError:
            # If there are no batches, return an empty array.
            return np.array([], dtype=str)

    def get_bit_caller(self):
        return BitCaller(self._rep.count_refs, self._rep.count_muts)

    def get_bit_counter(self):
        """ Return a BitCounter object of all filtered bit vectors. """
        return BitCounter(self.get_bit_caller(), self.iter_batches())

    def get_bit_vectors(self):
        """ Return a BitVectorSet object of all filtered bit vectors. """
        return BitVectorSet(self.get_bit_caller(), self.iter_batches())

    @classmethod
    def open(cls, report_file: Path):
        return cls(CallReport.open(report_file))
