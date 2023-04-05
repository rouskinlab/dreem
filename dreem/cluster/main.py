from click import command
import pandas as pd

from ..cluster.bitvector import BitVector
from ..cluster.clusteringAnalysis import ClusteringAnalysis
from ..cluster.EMclustering import EMclustering
from ..util.cli import (opt_report, opt_out_dir,
                        opt_max_procs,
                        opt_max_clusters, opt_num_runs, opt_signal_thresh,
                        opt_include_gu, opt_include_del, opt_polya_max,
                        opt_min_iter, opt_max_iter, opt_convergence_cutoff, opt_min_reads, opt_verbose)

from ..util import docdef, path


params = [
    # Input/output directories
    opt_report,
    opt_out_dir,
    # Parallelization
    opt_max_procs,
    # Clustering options
    opt_max_clusters,
    opt_num_runs,
    opt_signal_thresh,
    opt_include_gu,
    opt_include_del,
    opt_polya_max,
    opt_min_iter,
    opt_max_iter,
    opt_convergence_cutoff,
    opt_min_reads
]


@command("cluster", params=params)
def cli(*args, **kwargs):
    return run(*args, **kwargs)


@docdef.auto()
def run(mp_report: tuple[str], *,
        out_dir: str,
        max_procs: int,
        # Clustering options
        max_clusters: int,
        num_runs: int,
        signal_thresh: float,
        include_gu: bool,
        include_del: bool,
        polya_max: int,
        min_iter: int,
        max_iter: int,
        convergence_cutoff: float,
        min_reads: int):
    """
    Run the clustering module. 
    """

    cluster_out_files = list()

    # Create a clustering args 
    clustering_args = dict(
        max_clusters=max_clusters,
        min_iter=min_iter,
        max_iter=max_iter,
        signal_thresh=signal_thresh,
        include_gu=include_gu,
        include_del=include_del,
        min_reads=min_reads,
        convergence_cutoff=convergence_cutoff,
        num_runs=num_runs,
        max_procs=max_procs
    )

    # Get the bitvector files in the input directory and all of its subdirectories
    for report_file in mp_report:
        report_path = path.MutVectorReportFilePath.parse(report_file)

        # Run the clustering algorithm
        bitvector = BitVector(path=report_file, signal_thresh=signal_thresh, include_gu=include_gu, min_reads=min_reads,
                              include_del=include_del)
        ca = ClusteringAnalysis(bitvector, max_clusters, num_runs, clustering_args)
        clusters = ca.run()

        # Compute the likelihood of each read for each cluster
        columns = [(k, c + 1) for k in clusters for c in range(k)]
        reads_best_cluster = pd.DataFrame(dtype=float,
                                          index=pd.Index(bitvector.read_names, name="Read Name"),
                                          columns=pd.MultiIndex.from_tuples(columns,
                                                                            names=["K", "i"]))
        for k, clusters_k in clusters.items():
            em = EMclustering(bitvector.bv, k, bitvector.read_hist, bitvector.base_to_keep, bitvector.sequence,
                              **clustering_args)
            likelihood_reads_best_cluster, _, _ = em.expectation(clusters_k[0]['mu'], clusters_k[0]['pi'])
            for c in range(k):
                reads_best_cluster.loc[:, (k, c + 1)] = likelihood_reads_best_cluster[bitvector.read_inverse, c]

        # Save the results
        cluster_out_file = str(report_path.path.with_name("clusters.csv.gz"))
        reads_best_cluster.to_csv(cluster_out_file, header=True, index=True,
                                  compression="gzip")
        cluster_out_files.append(cluster_out_file)

    return cluster_out_files
