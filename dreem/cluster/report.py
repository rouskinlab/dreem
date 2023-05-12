from abc import ABC
from itertools import chain
from pathlib import Path

from .emalgo import EmClustering
from .metrics import (calc_bics, get_common_clusts, get_converged,
                      get_log_likes, get_log_like_mean, get_log_like_std,
                      get_var_info, find_best_k)
from ..util import path
from ..util.epu import get_common_attrib
from ..util.report import Report


class ClusterReport(Report, ABC):
    __slots__ = [
        "sample",
        "ref",
        "seq",
        "sect",
        "end5",
        "end3",
        "max_clust",
        "num_runs",
        "min_iter",
        "max_iter",
        "conv_thresh",
        "count_del",
        "count_ins",
        "exclude_gu",
        "exclude_polya",
        "exclude_pos",
        "min_ninfo_pos",
        "min_fmut_pos",
        "max_fmut_pos",
        "min_mut_gap",
        "min_finfo_read",
        "max_fmut_read",
        "n_pos_init",
        "n_pos_polya",
        "n_pos_gu",
        "n_pos_user",
        "n_min_ninfo_pos",
        "n_min_fmut_pos",
        "n_max_fmut_pos",
        "n_pos_kept",
        "n_reads_init",
        "n_min_finfo_read",
        "n_max_fmut_read",
        "n_min_mut_gap",
        "n_reads_kept",
        "n_uniq_reads_kept",
        "converged",
        "log_likes",
        "log_like_mean",
        "log_like_std",
        "var_info",
        "bic",
        "n_clust",
    ]

    @classmethod
    def from_clusters(cls, /,
                      clusters: dict[int, list[EmClustering]],
                      max_clusters: int,
                      num_runs: int):
        # Get common attributes of all the clusters.
        common = get_common_clusts(clusters)
        bvec = get_common_attrib(chain(*clusters.values()), "bvec", key=id)
        # Compute the log likelihoods
        log_likes = get_log_likes(clusters)
        # Initialize a new ClusterReport.
        return cls(max_clust=max_clusters,
                   num_runs=num_runs,
                   min_iter=get_common_attrib(common, "min_iter"),
                   max_iter=get_common_attrib(common, "max_iter"),
                   conv_thresh=get_common_attrib(common, "conv_thresh"),
                   n_clust=find_best_k(clusters),
                   bic=calc_bics(clusters),
                   converged=get_converged(clusters),
                   log_likes=log_likes,
                   log_like_mean=get_log_like_mean(log_likes),
                   log_like_std=get_log_like_std(log_likes),
                   var_info=get_var_info(clusters),
                   n_uniq_reads_kept=bvec.n_uniq,
                   **bvec.to_dict())

    def get_path(self, out_dir: Path):
        return path.build(path.ModSeg, path.SampSeg, path.RefSeg,
                          path.SectSeg, path.ClustRepSeg,
                          top=out_dir,
                          module=path.MOD_CLUST,
                          sample=self.sample,
                          ref=self.ref,
                          end5=self.end5,
                          end3=self.end3,
                          ext=path.JSON_EXT)
