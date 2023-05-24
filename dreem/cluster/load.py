from functools import cached_property
from pathlib import Path

import pandas as pd

from .report import ClusterReport
from ..call.load import CallVecLoader
from ..call.report import CallReport
from ..cluster.write import IDX_NCLUSTERS
from ..core import path
from ..core.mu import calc_mu, denom


class ClusterLoader(object):
    """ Load clustering results. """

    def __init__(self, *,
                 out_dir: Path,
                 sample: str,
                 ref: str,
                 sect: str):
        self.out_dir = out_dir
        self.sample = sample
        self.ref = ref
        self.sect = sect

    @cached_property
    def callvec(self):
        """ Return the CallLoader for this ClusterLoader. """
        return CallVecLoader.open(CallReport.build_path(self.out_dir,
                                                        sample=self.sample,
                                                        ref=self.ref,
                                                        sect=self.sect))

    @property
    def section(self):
        return self.callvec.section

    @property
    def min_mut_gap(self) -> int:
        return self.callvec.min_mut_gap

    @cached_property
    def resps_file(self):
        return path.build(path.ModSeg, path.SampSeg, path.RefSeg,
                          path.SectSeg, path.ClustTabSeg,
                          top=self.out_dir, module=path.MOD_CLUST,
                          sample=self.sample, ref=self.ref, sect=self.sect,
                          table=path.CLUST_RESP_TABLE, run=0,
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
    def mus(self):
        """ Return the bias-corrected mutation rates as a DataFrame with
        dimensions (positions x clusters). """
        bitvec = self.callvec.get_bit_vectors()
        # Compute the responsibility-weighted sums of the mutations and
        # matches at each position for each cluster.
        fmut = ((self.resps @ bitvec.muts) / (self.resps @ bitvec.info)).T
        return calc_mu(fmut, self.section, self.min_mut_gap)

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
