import json
import os
from click import command, pass_obj

from ..cluster.bitvector import BitVector
from ..cluster.clusteringAnalysis import ClusteringAnalysis
from ..cluster.EMclustering import EMclustering
from ..util.cli import (DreemCommandName, dreem_command,
                        opt_report, opt_out_dir,
                        opt_max_procs,
                        opt_max_clusters, opt_num_runs, opt_signal_thresh,
                        opt_include_gu, opt_include_del,
                        opt_min_iter, opt_max_iter, opt_convergence_cutoff, opt_min_reads, opt_verbose)

from ..util import docdef

@command(DreemCommandName.CLUSTER.value, params=[
    opt_verbose,
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
    opt_min_iter,
    opt_max_iter,
    opt_convergence_cutoff,
    opt_min_reads
])
# Pass context object
@pass_obj
# Turn into DREEM command
@dreem_command(imports=("mp_report"))
def cli(*args, **kwargs):
    return run(*args, **kwargs)

@docdef.auto()
def run(*,
        mp_report: tuple[str],
        out_dir: str,
        max_procs: int,
        # Clustering options
        max_clusters: int,
        num_runs: int,
        signal_thresh: float,
        include_gu: bool,
        include_del: bool,
        min_iter: int,
        max_iter: int,
        convergence_cutoff: float,
        min_reads: int,
        verbose: bool):
    """
    Run the clustering module. 
    """

    # Create the output file
    best_clusters_samples = {}

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
        max_procs=max_procs,
        verbose=verbose
    )


    # Get the bitvector files in the input directory and all of its subdirectories
    for i, report_file in enumerate(mp_report):  
        section = report_file.split('/')[-2]
        print("\n\nSTARTING SAMPLE", i, '|', section)

        # Run the clustering algorithm
        bitvector = BitVector(path=report_file, signal_thresh=signal_thresh, include_gu=include_gu, min_reads=min_reads, include_del=include_del)
        ca = ClusteringAnalysis(bitvector, max_clusters, num_runs, clustering_args)
        clusters = ca.run()

        # Compute the likelihood of each read for each cluster
        reads_best_cluster = {}
        for k in clusters:
            em = EMclustering(bitvector.bv, int(k[1]), bitvector.read_hist, bitvector.base_to_keep, bitvector.sequence,
                              **clustering_args)
            likelihood_reads_best_cluster, _, _ = em.expectation(clusters[k][0]['mu'], clusters[k][0]['pi'])
            reads_best_cluster[k] = bitvector.associate_reads_with_likelihoods(likelihood_reads_best_cluster)

        best_clusters_samples[section] = reads_best_cluster

    # Save the results
    with open(os.path.join(out_dir, 'best_cluster_reads.json'), 'w') as f:
        json.dump(best_clusters_samples, f, indent=4)

    return 1