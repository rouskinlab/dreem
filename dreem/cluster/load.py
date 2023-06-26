from functools import cached_property
from logging import getLogger
from pathlib import Path

import pandas as pd

from .indexes import ORD_CLS_NAME
from .report import ClustReport
from ..core import path
from ..core.bitvect import BitBatch, ClusterBitBatch
from ..core.load import BatchChainLoader, no_kwargs
from ..mask.load import MaskLoader

logger = getLogger(__name__)


def get_cluster_index(max_order: int):
    """ Get the index for clusters up to a given order. """
    return pd.MultiIndex.from_tuples((order, cluster)
                                     for order in range(1, max_order + 1)
                                     for cluster in range(1, order + 1))


class ClustLoader(BatchChainLoader):
    """ Load clustering results. """

    def __init__(self, report: ClustReport):
        super().__init__(report)
        self.best_order = report.best_order

    @classmethod
    def get_report_type(cls):
        return ClustReport

    @classmethod
    def get_import_type(cls):
        return MaskLoader

    @property
    def import_path_fields(self):
        return {**super().import_path_fields, path.SECT: self.sect}

    def get_refseq(self):
        return self.import_loader.get_refseq()

    @property
    def section(self):
        return self.import_loader.section

    @property
    def min_mut_gap(self) -> int:
        return self.import_loader.min_mut_gap

    @cached_property
    def clusters(self):
        """ Order and number of each cluster. """
        return get_cluster_index(self.best_order)

    @no_kwargs
    def load_data_personal(self, batch_file: Path):
        # Load the cluster memberships of the reads in the batch.
        data = pd.read_csv(batch_file,
                           index_col=[0],
                           header=list(range(len(ORD_CLS_NAME))))
        # Cast the cluster numbers from strings to integers.
        data.columns = pd.MultiIndex.from_tuples([(int(order), int(k))
                                                  for order, k in data.columns],
                                                 names=ORD_CLS_NAME)
        return data

    def iter_batches_personal(self):
        yield from super().iter_batches_personal()

    @no_kwargs
    def process_batch(self, imported_batch: BitBatch,
                      personal_batch: pd.DataFrame):
        return ClusterBitBatch(imported_batch.section,
                               imported_batch.info,
                               imported_batch.affi,
                               personal_batch)

    def iter_batches_processed(self, **kwargs):
        yield from super().iter_batches_processed(**kwargs)
