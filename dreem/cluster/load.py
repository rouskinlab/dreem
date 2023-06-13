from functools import cached_property
import re
from typing import Iterable

import numpy as np
import pandas as pd

from .report import ClustReport
from ..mask.load import MaskLoader
from ..mask.report import MaskReport
from ..core import path
from ..core.load import DataLoader
from ..core.mu import calc_mu_adj_df, calc_f_obs_df

IDX_NCLUSTERS = "NumClusters"
IDX_CLUSTER = "Cluster"
IDXS_CLUSTERS = IDX_NCLUSTERS, IDX_CLUSTER

# Cluster fields
CLUST_FORMAT = "Cluster {k}-{c}"
CLUST_PATTERN = re.compile("Cluster ([0-9]+)-([0-9]+)")
CLUST_NAME_IDX = "Cluster"
CLUST_PROP_COL = "Proportion"


def parse_names(names: Iterable[str]):
    """ Parse an iterable of cluster names formatted according to
    CLUST_PATTERN and return a MultiIndex of the number of clusters and
    the number of each cluster. """
    kc_pairs = [(int((kc := CLUST_PATTERN.match(n).groups())[0]), int(kc[1]))
                for n in names]
    return pd.MultiIndex.from_tuples(kc_pairs, names=IDXS_CLUSTERS)


def format_names(kc_pairs: Iterable[tuple[int, int]]):
    """ Given an iterable of tuples of the number of clusters and the
    number of the cluster, return for each tuple the name of the cluster
    formatted according to CLUST_PATTERN. """
    return pd.Index(CLUST_FORMAT.format(k=k, c=c) for k, c in kc_pairs)


def format_names_k(k: int):
    """ Return the name of every degree-k cluster. """
    return format_names((k, c) for c in range(1, k + 1))


def format_names_ks(ks: Iterable[int]):
    return pd.Index([name for k in ks for name in format_names_k(k)])


class ClustLoader(DataLoader):
    """ Load clustering results. """

    def __init__(self, report: ClustReport):
        super().__init__(report)
        self.n_clust = report.n_clust
        fields = dict(sample=report.sample, ref=report.ref, sect=report.sect)
        self._mask_load = MaskLoader.open(MaskReport.build_path(report.out_dir,
                                                                **fields))

    @classmethod
    def report_type(cls):
        return ClustReport

    @property
    def sect(self):
        return self._mask_load.sect

    @property
    def section(self):
        return self._mask_load.section

    @property
    def seq(self):
        return self._mask_load.seq

    @property
    def positions(self):
        return self._mask_load.positions

    @property
    def min_mut_gap(self) -> int:
        return self._mask_load.min_mut_gap

    @cached_property
    def bitvec(self):
        return self._mask_load.get_bit_monolith()

    def get_resps_file(self, k):
        return path.build(path.ModSeg, path.SampSeg, path.RefSeg,
                          path.SectSeg, path.ClustTabSeg,
                          top=self.out_dir, module=path.MOD_CLUST,
                          sample=self.sample, ref=self.ref, sect=self.sect,
                          table=path.CLUST_RESP_RUN_TABLE, k=k, run=0,
                          ext=path.CSVZIP_EXT)

    @cached_property
    def resps(self):
        """ Return the responsibilities (i.e. cluster memberships) as a
        DataFrame with dimensions (clusters x reads). """
        # Read the responsibilities from the CSV file for every k.
        resps_list = list()
        for k in range(1, self.n_clust + 1):
            # The responsibilities are stored as logs, so use np.exp().
            rk = np.exp(pd.read_csv(self.get_resps_file(k), index_col=0))
            # Add the number of clusters to the columns.
            rk.columns = format_names(zip([k] * k,
                                          rk.columns.astype(int),
                                          strict=True))
            # Add the responsibilities for this k to the list.
            resps_list.append(rk)
        # Concatenate the responsibilities of all k values.
        return pd.concat(resps_list, axis=1)

    @property
    def clusters(self):
        return self.resps.columns

    @cached_property
    def ninfo_per_pos(self):
        return self.bitvec.info.T @ self.resps

    @cached_property
    def nmuts_per_pos(self):
        return self.bitvec.muts.T @ self.resps

    @cached_property
    def mus(self):
        return calc_mu_adj_df(self.nmuts_per_pos / self.ninfo_per_pos,
                              self.section,
                              self.min_mut_gap)

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
