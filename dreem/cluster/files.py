from abc import ABC
from itertools import chain
from logging import getLogger
from pathlib import Path
from typing import Callable, Iterable

import pandas as pd

from .emalgo import EmClustering
from .metrics import (calc_bics, get_common_clusts, get_converged,
                      get_log_likes, get_log_like_mean, get_log_like_std,
                      get_var_info, find_best_k)
from ..util import path
from ..util.epu import get_common_attrib
from ..util.report import Report, DECIMAL_PRECISION
from ..vector.bits import BitVector

logger = getLogger(__name__)

IDX_NCLUSTERS = "NumClusters"
IDX_CLUSTER = "Cluster"
IDXS_CLUSTERS = IDX_NCLUSTERS, IDX_CLUSTER


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
        "n_clust",
        "bic",
        "converged",
        "log_likes",
        "log_like_mean",
        "log_like_std",
        "var_info",
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


def write_results(clusters: dict[int, list[EmClustering]],
                  out_dir: Path):
    """ Write CSV files of the proportions, mutation rates, and read
    responsibilities for each cluster. Return the file paths. """
    bvec = get_common_attrib(chain(*clusters.values()), "bvec", key=id)
    # Proportions: proportion of each cluster in the ensemble
    props = write_tables(build_tables(clusters, EmClustering.output_props),
                         bvec, out_dir, path.CLUST_PROP_TABLE)
    # Mutation rates: fraction of mutated bits at each position
    mus = write_tables(build_tables(clusters, EmClustering.output_mus),
                       bvec, out_dir, path.CLUST_MUS_TABLE)
    # Responsibilities: likelihood that each read came from each cluster
    resps = write_tables(build_tables(clusters, EmClustering.output_resps),
                         bvec, out_dir, path.CLUST_RESP_TABLE, gzip=True)
    return props, mus, resps


def build_tables(clusters: dict[int, list[EmClustering]],
                 output_func: Callable[[EmClustering], pd.DataFrame]):
    """ Build a DataFrame of one attribute of one or more clusters. """
    tables: dict[int, pd.DataFrame] = dict()
    for k, runs in clusters.items():
        for r, run in enumerate(runs):
            ktable = output_func(run)
            if r not in tables:
                tables[r] = pd.DataFrame(index=ktable.index,
                                         columns=kc_pairs(clusters),
                                         dtype=float)
            for c, column in ktable.items():
                tables[r].loc[:, (k, c)] = column
    return tables


def write_tables(tables: dict[int, pd.DataFrame],
                 bvec: BitVector,
                 out_dir: Path,
                 name: str,
                 gzip: bool = False) -> dict[int, Path]:
    return {run: write_table(table, bvec, out_dir, name, run, gzip)
            for run, table in tables.items()}


def write_table(table: pd.DataFrame,
                bvec: BitVector,
                out_dir: Path,
                name: str,
                run: int,
                gzip: bool = False):
    """ Write a DataFrame of one clustering attribute to a CSV file. """
    file = path.build(path.ModSeg, path.SampSeg, path.RefSeg,
                      path.SectSeg, path.ClustTabSeg,
                      top=out_dir, module=path.MOD_CLUST,
                      sample=bvec.loader.sample, ref=bvec.loader.ref,
                      end5=bvec.section.end5, end3=bvec.section.end3,
                      table=name, run=run,
                      ext=(path.CSVZIP_EXT if gzip else path.CSV_EXT))
    table.round(DECIMAL_PRECISION).to_csv(file, header=True, index=True)
    logger.info(f"Wrote {name} run {run} of {bvec} to {file}")
    return file


def load_cluster_resps(cluster_resps_file: Path):
    """ Load the responsibilities of the reads for each cluster. """
    return pd.read_csv(cluster_resps_file,
                       index_col=[0],
                       header=list(range(len(IDXS_CLUSTERS))))
