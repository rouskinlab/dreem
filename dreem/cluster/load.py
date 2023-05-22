from functools import cached_property
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .report import ClusterReport
from dreem.quant.vect import BitVector
from ..mvec.load import MutVecLoader
from ..mvec.write import MutVecReport
from ..quant.bias import denom
from ..quant.mu import mus_obs_to_real
from ..util import path
from ..util.sect import Section
from ..util.seq import DNA


class ClusterLoader(object):
    """ Load clustering results. """

    def __init__(self, *,
                 out_dir: Path,
                 sample: str,
                 ref: str,
                 seq: DNA,
                 sect: str,
                 end5: int,
                 end3: int,
                 count_del: bool,
                 count_ins: bool,
                 exclude_polya: int,
                 exclude_gu: bool,
                 exclude_pos: Iterable[int],
                 min_mut_gap: int):
        self.out_dir = out_dir
        self.sample = sample
        self.ref = ref
        self.seq = seq
        self.sect = sect
        self.end5 = end5
        self.end3 = end3
        self.count_del = count_del
        self.count_ins = count_ins
        self.exclude_polya = exclude_polya
        self.exclude_gu = exclude_gu
        self.exclude_pos = exclude_pos
        self.min_mut_gap = min_mut_gap

    @classmethod
    def open(cls, report_file: Path):
        """ Create a ClusterLoader from a clustering report file. """
        rep = ClusterReport.open(report_file)
        # Exclude positions that were excluded during clustering, except
        # for those with a mutation rate below the minimum because those
        # positions are only dropped within the EM algorithm to speed up
        # clustering.
        exclude_pos = rep.exclude_pos + rep.pos_min_ninfo + rep.pos_max_fmut
        return cls(out_dir=path.parse(report_file,
                                      *ClusterReport.path_segs())[path.TOP],
                   sample=rep.sample, ref=rep.ref, seq=rep.seq,
                   sect=rep.sect, end5=rep.end5, end3=rep.end3,
                   count_del=rep.count_del, count_ins=rep.count_ins,
                   exclude_polya=rep.exclude_polya,
                   exclude_gu=rep.exclude_gu,
                   exclude_pos=exclude_pos,
                   min_mut_gap=rep.min_mut_gap)

    @cached_property
    def vector_loader(self):
        loader = MutVecLoader.open(MutVecReport.build_path(self.out_dir,
                                                           sample=self.sample,
                                                           ref=self.ref))
        loader_seq = loader.section(self.end5, self.end3).seq
        if loader_seq != self.seq:
            raise ValueError(f"Sequences differ in cluster report ({self.seq}) "
                             f"and vector report ({loader_seq})")
        return loader

    @cached_property
    def section(self):
        return Section(self.ref, self.vector_loader.seq,
                       end5=self.end5, end3=self.end3, name=self.sect)

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

    @cached_property
    def mus(self):
        """ Return the bias-corrected mutation rates as a DataFrame with
        dimensions (positions x clusters). """
        bvec = BitVector(self.vector_loader, self.section,
                         count_del=self.count_del,
                         count_ins=self.count_ins,
                         exclude_polya=self.exclude_polya,
                         exclude_gu=self.exclude_gu,
                         exclude_pos=self.exclude_pos)
        # Compute the responsibility-weighted sums of the mutations and
        # matches at each position for each cluster.
        mut_sums = (self.resps @ bvec.all_muts().loc[self.resps.columns]).T
        ref_sums = (self.resps @ bvec.all_refs().loc[self.resps.columns]).T
        # Determine the positions and clusters.
        pos_all = pd.Index(self.section.columns)
        pos_use = mut_sums.index
        pos_cut = pos_all.drop(pos_use)
        clusters = self.resps.index
        # Initialize to zero the observed mutation rate of each position
        # in each cluster.
        mus_obs = pd.DataFrame(0., index=pos_all, columns=clusters)
        # Set the observed mutation rates of the used positions.
        mus_obs.loc[pos_use] = mut_sums / (mut_sums + ref_sums)
        # Correct the bias in the observed mutation rates.
        mus_real = mus_obs_to_real(mus_obs, self.min_mut_gap)
        # Mask the unused positions with NaN.
        mus_real.loc[pos_cut] = np.nan
        return mus_real

    @cached_property
    def denoms(self):
        """ Return the denominators of all clusters as a Series with
        dimension (clusters). """
        return pd.Series(denom(self.mus.fillna(0.).values, self.min_mut_gap),
                         index=self.resps.index)

    @cached_property
    def props(self) -> pd.Series:
        """ Return the bias-corrected cluster proportions as a Series
        with dimension (clusters). """
        # Compute the observed cluster proportions, then divide by the
        # cluster denominators to adjust for drop-out.
        props = self.resps.mean(axis=1) / self.denoms
        # Normalize the proportions so that those in each cluster number
        # sum to unity.
        return props / props.groupby(level=["NumClusters"]).sum()
