from functools import cached_property
from logging import getLogger
from pathlib import Path

import pandas as pd

from .indexes import ORD_CLS_NAME
from .report import ClustReport
from ..core.bitv import BitBatch, ClusterBitBatch
from ..core.load import SectBatchChainLoader, no_kwargs
from ..mask.load import MaskLoader

logger = getLogger(__name__)


def get_cluster_index(max_order: int):
    """ Get the index for clusters up to a given order. """
    return pd.MultiIndex.from_tuples((order, cluster)
                                     for order in range(1, max_order + 1)
                                     for cluster in range(1, order + 1))


class ClustLoader(SectBatchChainLoader):
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
    def import_kwargs(self):
        return dict()

    def get_refseq(self):
        return self.import_loader.get_refseq()

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
    def process_batch(self,
                      imported_batch: BitBatch,
                      personal_batch: pd.DataFrame):
        return ClusterBitBatch(imported_batch.all, imported_batch.yes,
                               personal_batch)

    def iter_batches_processed(self, **kwargs):
        yield from super().iter_batches_processed(**kwargs)

    '''

    @cached_property
    def f_obs(self):
        """ Return the fraction observed for every cluster as a Series
        with dimension (clusters). """
        return calc_f_obs_df(self.mus, self.section, self.min_mut_gap)

    @cached_property
    def props(self) -> pd.Series:
        """ Return the bias-corrected cluster proportions as a Series
        with dimension (clusters). """
        # Compute the observed cluster proportions, then divide by the
        # cluster denominators to adjust for drop-out.
        props = self.resps.mean(axis=0) / self.f_obs
        # Convert the index into a MultiIndex so that the clusters can
        # be grouped by their number of clusters.
        props.index = parse_names(props.index)
        # Normalize the proportions so that those in each cluster number
        # sum to unity.
        props /= props.groupby(level=[IDX_NCLUSTERS]).sum()
        # Convert the index back into string labels.
        props.index = format_names(props.index)
        return props

'''
