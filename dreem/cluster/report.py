from abc import ABC
from itertools import chain

from .emalgo import EmClustering
from .metrics import (calc_bics, get_common_clusts, get_converged,
                      get_log_likes, get_log_like_mean, get_log_like_std,
                      get_var_info, find_best_k)
from ..util import path
from ..util.epu import get_common_attrib
from ..util.report import Report


class ClusterReport(Report, ABC):
    __slots__ = [
        # Sample, reference, and section information.
        "sample", "ref", "seq", "sect", "end5", "end3",
        # Types of mutations to count.
        "count_del", "count_ins",
        # Position filtering parameters.
        "exclude_gu", "exclude_polya", "exclude_pos",
        "min_ninfo_pos", "min_fmut_pos", "max_fmut_pos",
        # Position filtering results.
        "n_pos_init",
        "n_pos_gu", "n_pos_polya", "n_pos_user",
        "n_pos_min_ninfo", "n_pos_min_fmut", "n_pos_max_fmut",
        "n_pos_kept",
        "pos_min_ninfo", "pos_min_fmut", "pos_max_fmut",
        # Read filtering parameters.
        "min_finfo_read", "max_fmut_read", "min_mut_gap",
        # Read filtering results.
        "n_reads_init", "n_reads_min_finfo", "n_reads_max_fmut",
        "n_reads_min_gap", "n_reads_kept", "n_uniq_reads_kept",
        # Clustering parameters.
        "max_clust", "num_runs", "min_iter", "max_iter", "conv_thresh",
        # Clustering results.
        "converged", "log_likes", "log_like_mean", "log_like_std", "var_info",
        "bic", "n_clust",
    ]

    @classmethod
    def path_segs(cls):
        return super().path_segs() + [path.SectSeg, path.ClustRepSeg]

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.MOD: path.MOD_CLUST}

    @classmethod
    def from_clusters(cls, /,
                      clusters: dict[int, list[EmClustering]],
                      max_clusters: int,
                      num_runs: int):
        """ Create a ClusterReport from EmClustering objects. """
        # Get common attributes of all the clusters.
        common = get_common_clusts(clusters)
        bvec = get_common_attrib(chain(*clusters.values()), "bvec", key=id)
        # Compute the log likelihoods.
        log_likes = get_log_likes(clusters)
        # Initialize a new ClusterReport.
        return cls(max_clust=max_clusters,
                   num_runs=num_runs,
                   n_clust=find_best_k(clusters),
                   bic=calc_bics(clusters),
                   converged=get_converged(clusters),
                   log_likes=log_likes,
                   log_like_mean=get_log_like_mean(log_likes),
                   log_like_std=get_log_like_std(log_likes),
                   var_info=get_var_info(clusters),
                   n_uniq_reads_kept=bvec.n_uniq,
                   **{attrib: get_common_attrib(common, attrib)
                      for attrib in ["min_iter", "max_iter", "conv_thresh"]},
                   **bvec.to_dict())
