from functools import cached_property
from pathlib import Path

import pandas as pd

from .report import ClusterReport
from ..call.load import BitVecLoader
from ..call.report import CallReport
from ..cluster.write import IDX_NCLUSTERS
from ..core import path
from ..core.mu import calc_mus, denom


class ClusterLoader(object):
    """ Load clustering results. """

    def __init__(self, *,
                 out_dir: Path,
                 sample: str,
                 ref: str,
                 sect: str):
        self._bv_load = BitVecLoader.open(CallReport.build_path(out_dir,
                                                                sample=sample,
                                                                ref=ref,
                                                                sect=sect))

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

    @cached_property
    def resps_file(self):
        return path.build(path.ModSeg, path.SampSeg, path.RefSeg,
                          path.SectSeg, path.ClustTabSeg,
                          top=self.out_dir, module=path.MOD_CLUST,
                          sample=self.sample, ref=self.ref, sect=self.sect,
                          table=path.CLUST_RESP_RUN_TABLE, run=0,
                          ext=path.CSVZIP_EXT)

    @cached_property
    def resps(self):
        """ Return the responsibilities (a.k.a. cluster memberships) as
        a DataFrame with dimensions (clusters x reads). """
        # The index has two levels: number of clusters and cluster
        levels = list(range(2))
        # Read the responsibilities from the CSV file.
        resps = pd.read_csv(self.resps_file, index_col=0, header=levels).T
        # Convert the index from str (the default) to int.
        idx_arrays = [list(map(int, resps.index.get_level_values(i)))
                      for i in levels]
        resps.index = pd.MultiIndex.from_arrays(idx_arrays,
                                                names=resps.index.names)
        return resps

    @property
    def clusters(self):
        return self.resps.index

    @cached_property
    def bitvec(self):
        return self._bv_load.get_bit_vectors()

    @cached_property
    def ninfo_per_pos(self):
        return self.bitvec.info.T @ self.resps.T

    @cached_property
    def nmuts_per_pos(self):
        return self.bitvec.muts.T @ self.resps.T

    @cached_property
    def mus(self):
        return calc_mus(self.nmuts_per_pos / self.ninfo_per_pos,
                        self.section,
                        self.min_mut_gap)

    @cached_property
    def denoms(self):
        """ Return the denominators of all clusters as a Series with
        dimension (clusters). """
        return pd.Series(denom(self.mus.fillna(0.).values, self.min_mut_gap),
                         index=self.clusters)

    @cached_property
    def props(self) -> pd.Series:
        """ Return the bias-corrected cluster proportions as a Series
        with dimension (clusters). """
        # Compute the observed cluster proportions, then divide by the
        # cluster denominators to adjust for drop-out.
        props = self.resps.mean(axis=1) / self.denoms
        # Normalize the proportions so that those in each cluster number
        # sum to unity.
        return props / props.groupby(level=[IDX_NCLUSTERS]).sum()

    @classmethod
    def open(cls, report_file: Path):
        """ Create a ClusterLoader from a clustering report file. """
        report = ClusterReport.open(report_file)
        return cls(out_dir=report.out_dir, sample=report.sample,
                   ref=report.ref, sect=report.sect)
