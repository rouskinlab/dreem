from functools import cache, cached_property
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .report import ClusterReport
from ..call.load import BitVecLoader
from ..call.report import CallReport
from ..core import path
from ..core.mu import calc_mu_df, calc_f_obs

IDX_NCLUSTERS = "NumClusters"
IDX_CLUSTER = "Cluster"
IDXS_CLUSTERS = IDX_NCLUSTERS, IDX_CLUSTER


def kc_pairs(ks: Iterable[int]):
    """
    Return the multi-index for a data frame with each number of clusters
    in ```ks```. The first level of the multi-index is the total number
    of clusters, ```k```; the second level is the number of the cluster,
    ```c```. For example, ```(3, 2)``` signifies cluster number 2 from
    the best run of EM clustering with ```k``` = 3 clusters in total.
    Note that the length of the multi-index equals ```sum(ks)```.

    Examples
    --------
    >>> kc_pairs([1, 2, 3]).tolist()
    [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3)]
    """
    return pd.MultiIndex.from_tuples([(k, c) for k in ks
                                      for c in range(1, k + 1)],
                                     names=IDXS_CLUSTERS)


class ClusterLoader(object):
    """ Load clustering results. """

    def __init__(self, *,
                 out_dir: Path,
                 sample: str,
                 ref: str,
                 sect: str,
                 n_clust: int):
        self._bv_load = BitVecLoader.open(CallReport.build_path(out_dir,
                                                                sample=sample,
                                                                ref=ref,
                                                                sect=sect))
        self.n_clust = n_clust

    @property
    def out_dir(self):
        return self._bv_load.out_dir

    @property
    def sample(self):
        return self._bv_load.sample

    @property
    def ref(self):
        return self._bv_load.ref

    @property
    def sect(self):
        return self._bv_load.sect

    @property
    def section(self):
        return self._bv_load.section

    @property
    def seq(self):
        return self._bv_load.seq

    @property
    def positions(self):
        return self._bv_load.positions

    @property
    def min_mut_gap(self) -> int:
        return self._bv_load.min_mut_gap

    def get_resps_file(self, k):
        return path.build(path.ModSeg, path.SampSeg, path.RefSeg,
                          path.SectSeg, path.ClustTabSeg,
                          top=self.out_dir, module=path.MOD_CLUST,
                          sample=self.sample, ref=self.ref, sect=self.sect,
                          table=path.CLUST_RESP_RUN_TABLE, k=k, run=0,
                          ext=path.CSVZIP_EXT)

    @cache
    def _get_resps(self):
        """ Return the responsibilities (i.e. cluster memberships). """
        # Read the responsibilities from the CSV file for every k.
        resps_list = list()
        for k in range(1, self.n_clust + 1):
            # The responsibilities are stored as logs, so use np.exp().
            rk = np.exp(pd.read_csv(self.get_resps_file(k), index_col=0))
            if rk.shape[1] != k:
                raise ValueError(f"Expected resps to have {k} columns, "
                                 f"but got {rk.shape[1]}")
            # Add the number of clusters to the columns.
            rk.columns = pd.MultiIndex.from_arrays([np.full(k, k, dtype=int),
                                                    rk.columns.astype(int)],
                                                   names=IDXS_CLUSTERS)
            # Add the responsibilities for this k to the list.
            resps_list.append(rk)
        # Concatenate the responsibilities of all k values.
        resps = pd.concat(resps_list, axis=1)
        # Return the values, index, and columns separately.
        return resps.values, resps.index, resps.columns

    @property
    def resps(self):
        """ Return the responsibilities (a.k.a. cluster memberships) as
        a DataFrame with dimensions (clusters x reads). """
        # Get the values, index, and columns of the resps DataFrame.
        values, index, columns = self._get_resps()
        # Reassemble and return the DataFrame.
        # The purpose of this disassembly-reassembly process is to let
        # multiple functions each use a unique instance of the DataFrame
        # while sharing the same instance of the data (i.e. the values
        # variable, which is cached by the function self._get_resps()
        # and so does not need to be re-loaded from the file) to minimize
        # the memory and I/O requirements.
        # For example, the table-writing functions need two different
        # types of columns for the DataFrame and thus cannot share the
        # same instance of it without conflicts.
        # This disassemble-cache-reassemble method solves that problem.
        return pd.DataFrame(values, index=index, columns=columns, copy=False)

    @property
    def clusters(self):
        return self.resps.columns

    @cached_property
    def bitvec(self):
        return self._bv_load.get_bit_vectors()

    @cached_property
    def ninfo_per_pos(self):
        return self.bitvec.info.T @ self.resps

    @cached_property
    def nmuts_per_pos(self):
        return self.bitvec.muts.T @ self.resps

    @cached_property
    def mus(self):
        return calc_mu_df(self.nmuts_per_pos / self.ninfo_per_pos,
                          self.section,
                          self.min_mut_gap)

    @cached_property
    def f_obs(self):
        """ Return the fraction observed for every cluster as a Series
        with dimension (clusters). """
        return pd.Series(calc_f_obs(self.mus.fillna(0.).values, self.min_mut_gap),
                         index=self.clusters)

    @cached_property
    def props(self) -> pd.Series:
        """ Return the bias-corrected cluster proportions as a Series
        with dimension (clusters). """
        # Compute the observed cluster proportions, then divide by the
        # cluster denominators to adjust for drop-out.
        props = self.resps.mean(axis=0) / self.f_obs
        # Normalize the proportions so that those in each cluster number
        # sum to unity.
        print(self.resps)
        print(props)
        return props / props.groupby(level=[IDX_NCLUSTERS]).sum()

    @classmethod
    def open(cls, report_file: Path):
        """ Create a ClusterLoader from a clustering report file. """
        report = ClusterReport.open(report_file)
        return cls(out_dir=report.out_dir, sample=report.sample,
                   ref=report.ref, sect=report.sect, n_clust=report.n_clust)
