import json
import os
from click import command, pass_obj

from ..cluster.bitvector import BitVector
from ..cluster.clusteringAnalysis import ClusteringAnalysis
from ..cluster.EMclustering import EMclustering
from ..util.cli import (DreemCommandName, dreem_command,
                        opt_report, opt_out_dir, opt_temp_dir,
                        opt_rerun, opt_resume, opt_save_temp,
                        opt_parallel, opt_max_procs,
                        opt_max_clusters, opt_num_runs, opt_signal_thresh,
                        opt_info_thresh, opt_include_gu, opt_include_del,
                        opt_min_iter, opt_convergence_cutoff, opt_min_reads)


@command(DreemCommandName.CLUSTER.value, params=[
    # Input files
    opt_report,
    # Output directories
    opt_out_dir,
    opt_temp_dir,
    # File generation
    opt_rerun,
    opt_resume,
    opt_save_temp,
    # Parallelization
    opt_parallel,
    opt_max_procs,
    # Clustering options
    opt_max_clusters,
    opt_num_runs,
    opt_signal_thresh,
    opt_info_thresh,
    opt_include_gu,
    opt_include_del,
    opt_min_iter,
    opt_convergence_cutoff,
    opt_min_reads,
])
# Pass context object
@pass_obj
# Turn into DREEM command
@dreem_command(imports=("mp_report",))
def cli(*args, **kwargs):
    return run(*args, **kwargs)


def run(report_files: tuple[str],
        n_cpus: int,
        out_dir: str,
        max_clusters: int,
        min_iter: int,
        signal_thresh: float,
        info_thresh: float,
        include_gu: bool,
        include_del: bool,
        min_reads: int,
        convergence_cutoff: float,
        num_runs: int,
        verbose: bool):
    """
    Run the clustering module.

    Clusters the reads of all given bitvectors and outputs the likelihoods of the clusters as `name`.json in the directory `output_path`, using `temp_dir` as a temp directory.
    Each bitvector is a file containing the reads of a construct. Bitvectors from the same sample should be grouped in a folder and the path to the folder should be given as `bv_dir`.
    `name` is the name of the output file, and should be the sample name.

    Parameters from args:
    ---------------------

    input_dir: str
        Path to the bit vector folder or list of paths to the bit vector folders.
    out_dir: str
        Path to the output folder.
    max_clusters: int
        Maximum number of clusters.
    min_iter: int
        Minimum number of iteration per EM exectution.
    signal_thresh: float
        Signal threshold
    info_thresh: float
        Float from 0 to 1, where 1 means that all bases are unvalid and 0 means that all bases are valid (valid:= just 0s and 1s in the bitvector). If info_thresh of a read is above this threshold, it is considered unvalid and isn't used.
    include_g_u: bool
        Include G and U
    include_del: bool
        Include deletions
    min_reads: int
        Minimum number of reads per cluster.
    convergence_cutoff: float
        Convergence cutoff
    num_runs: int
        Number of runs
    n_cpus: int
        Number of cpus
    verbose: bool
        Verbose
        
    Returns
    -------
    
    clustering: dict
        The structure is the following:
        {reference_1}: # name of the bitvector file (from fasta file)
            {section_1}: # name of the clustered section (from library)
                {read_1}:  # read number (from bitvector)
                    K1_1: likelihood that this read belongs to cluster 1 when using 1 cluster
                    K2_1: likelihood that this read belongs to cluster 1 when using 2 clusters
                    K2_2: likelihood that this read belongs to cluster 2 when using 2 clusters
                    K3_1: likelihood that this read belongs to cluster 1 when using 3 clusters
                    K3_2: likelihood that this read belongs to cluster 2 when using 3 clusters
                    K3_3: likelihood that this read belongs to cluster 3 when using 3 clusters   
    
    """

    # Create the output file
    best_clusters_samples = {}

    # Create a clustering args 
    clustering_args = dict(
        max_clusters=max_clusters,
        min_iter=min_iter,
        signal_thresh=signal_thresh,
        info_thresh=info_thresh,
        include_g_u=include_gu,
        include_del=include_del,
        min_reads=min_reads,
        convergence_cutoff=convergence_cutoff,
        num_runs=num_runs,
        n_cpus=n_cpus,
        verbose=verbose
    )






    # Get the bitvector files in the input directory and all of its subdirectories
    for i, report_file in enumerate(report_files):
        
        ## PLACEHOLDER FOR THE FLAKE8
        f_in = None #TODO
        ############################
        
        section = f_in.split('/')[-2]
        print("\n\nSTARTING SAMPLE", i, '|', section)
        bitvector = BitVector(path=f_in)
        bitvector.publish_preprocessing_report(path=os.path.join(out_dir, section + '_preprocessing_report.txt'))
        ca = ClusteringAnalysis(bitvector, max_clusters, num_runs, clustering_args)
        clusters = ca.run()
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
