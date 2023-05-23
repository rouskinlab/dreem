from .emalgo import EmClustering
from .metric import (calc_bics, get_converged, get_log_likes,
                     get_log_like_mean, get_log_like_std,
                     get_var_info, find_best_k)
from ..call.load import CallLoader
from ..core import path
from ..core.bit import UniqMutBits
from ..core.report import Report


class ClusterReport(Report):
    __slots__ = (
        # Sample, reference, and section information.
        "sample", "ref", "seq", "sect", "end5", "end3", "n_uniq_reads_kept",
        # Clustering parameters.
        "max_clust", "num_runs", "min_iter", "max_iter", "conv_thresh",
        # Clustering results.
        "converged", "log_likes", "log_like_mean", "log_like_std", "var_info",
        "bic", "n_clust",
    )

    @classmethod
    def path_segs(cls):
        return super().path_segs() + (path.SectSeg, path.ClustRepSeg)

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.MOD: path.MOD_CLUST}

    @classmethod
    def from_clusters(cls, /,
                      loader: CallLoader,
                      uniq_muts: UniqMutBits,
                      clusters: dict[int, list[EmClustering]],
                      max_clusters: int,
                      num_runs: int, *,
                      min_iter: int,
                      max_iter: int,
                      conv_thresh: float):
        """ Create a ClusterReport from EmClustering objects. """
        # Get the log likelihoods.
        log_likes = get_log_likes(clusters)
        # Initialize a new ClusterReport.
        return cls(out_dir=loader.out_dir,
                   sample=loader.sample,
                   ref=loader.ref,
                   seq=loader.seq,
                   sect=loader.sect,
                   end5=loader.end5,
                   end3=loader.end3,
                   n_uniq_reads_kept=uniq_muts.n_uniq,
                   max_clust=max_clusters,
                   num_runs=num_runs,
                   min_iter=min_iter,
                   max_iter=max_iter,
                   conv_thresh=conv_thresh,
                   converged=get_converged(clusters),
                   log_likes=log_likes,
                   log_like_mean=get_log_like_mean(log_likes),
                   log_like_std=get_log_like_std(log_likes),
                   var_info=get_var_info(clusters),
                   bic=calc_bics(clusters),
                   n_clust=find_best_k(clusters))
