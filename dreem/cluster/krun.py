from concurrent.futures import Future, ProcessPoolExecutor
from logging import getLogger
from math import inf
from pathlib import Path

from .emalgo import EmClustering
from .files import ClusterReport, write_results
from .metrics import find_best_k
from ..util.parallel import get_num_parallel
from ..vector.bits import BitVector, VectorFilter
from ..vector.load import VectorLoader

logger = getLogger(__name__)


def cluster_sect(loader: VectorLoader, end5: int, end3: int, *,
                 exclude_gu: bool, include_del: bool, include_ins: bool,
                 exclude_polya: int, max_muts_per_read: int, min_gap: int,
                 min_reads: int, info_thresh: float, signal_thresh: float,
                 max_clusters: int, n_runs: int, conv_thresh: float,
                 min_iter: int, max_iter: int, max_procs: int, out_dir: Path):
    """ Cluster a section of a set of bit vectors. """
    # Generate the bit vectors from the mutation vector report.
    try:
        filter_vec = VectorFilter(min_mut_gap=min_gap,
                                  min_finfo_read=info_thresh,
                                  max_fmut_read=max_muts_per_read,
                                  min_ninfo_pos=min_reads,
                                  min_fmut_pos=signal_thresh,
                                  max_fmut_pos=1.0)  # FIXME: add this parameter to CLI/API
        bvec = BitVector(loader, end5=end5, end3=end3,
                         count_del=include_del,
                         count_ins=include_ins,
                         exclude_polya=exclude_polya,
                         exclude_gu=exclude_gu,
                         exclude_pos=list(),  # FIXME: add this parameter to CLI/API
                         filter_vec=filter_vec)
    except Exception as error:
        logger.critical(f"Failed to generate bit vectors: {error}")
        return
    # Cluster the bit vectors.
    logger.info(f"Began clustering {bvec} with up to k = {max_clusters} "
                f"clusters and {n_runs} replicate runs per number of clusters")
    try:
        clusts = run_max_clust(bvec, max_clusters, n_runs,
                               conv_thresh=conv_thresh,
                               min_iter=min_iter,
                               max_iter=max_iter,
                               max_procs=max_procs)
    except Exception as error:
        logger.critical(f"Failed to cluster {bvec}: {error}")
        return
    logger.info(f"Ended clustering {bvec}: got {find_best_k(clusts)} clusters")
    # Output the results of clustering.
    try:
        report = ClusterReport.from_clusters(clusts, max_clusters, n_runs)
        report_file = report.save(out_dir)
        write_results(clusts, out_dir)
    except Exception as error:
        logger.critical(f"Failed to write clusters for {bvec}: {error}")
        return
    # Return the path to the file of read responsibilities.
    return report_file


def run_max_clust(bvec: BitVector, max_clusters: int, n_runs: int, *,
                  conv_thresh: float, min_iter: int, max_iter: int,
                  max_procs: int):
    """
    Find the optimal number of clusters for EM, up to max_clusters.
    """
    if n_runs < 1:
        logger.warning(f"Number of EM runs must be ≥ 1: setting to 1")
        n_runs = 1
    logger.info(f"Began clustering {bvec} with up to k = {max_clusters} "
                f"clusters and {n_runs} replicate runs per number of clusters")
    # Store the clustering runs for each number of clusters.
    runs: dict[int, list[EmClustering]] = dict()
    # Run clustering for each number of clusters, up to max_clusters.
    while len(runs) < max_clusters:
        # Get the number of clusters, k.
        k = len(runs) + 1
        # Run EM clustering n_runs times with different starting points.
        runs[k] = run_n_clust(bvec, k,
                              n_runs=(n_runs if k > 1 else 1),
                              conv_thresh=(conv_thresh if k > 1 else inf),
                              min_iter=(min_iter if k > 1 else 1),
                              max_iter=(max_iter if k > 1 else 2),
                              max_procs=max_procs)
        # Find the best (smallest) BIC obtained from clustering.
        best_bic = runs[k][0].bic
        logger.debug(f"The best BIC with {k} cluster(s) is {best_bic}")
        if k > 1:
            prev_bic = runs[k - 1][0].bic
            # Check if the best BIC is better (smaller) than the best
            # BIC from clustering with one cluster fewer.
            if best_bic < prev_bic:
                # If the best BIC is < the previous best BIC, then the
                # current model is the best.
                logger.debug(f"The BIC decreased from {best_bic} "
                             f"(k = {k - 1}) to {best_bic} "
                             f"(k = {k})")
            else:
                # If the best BIC is ≥ than the previous best BIC, then
                # this model is worse than the previous model. Stop.
                logger.info(f"The BIC increased from {prev_bic} "
                            f"(k = {k - 1}) to {best_bic} "
                            f"(k = {k}): stopping")
                break
    return runs


def run_n_clust(bvec: BitVector, n_clusters: int, n_runs: int, *,
                conv_thresh: float, min_iter: int, max_iter: int,
                max_procs: int) -> list[EmClustering]:
    """ Run EM with a specific number of clusters. """
    logger.info(f"Began EM with {n_clusters} clusters ({n_runs} runs)")
    # Initialize one EmClustering object for each replicate run.
    reps = [EmClustering(bvec, n_clusters, rep,
                         min_iter=min_iter, max_iter=max_iter,
                         conv_thresh=conv_thresh)
            for rep in range(1, n_runs + 1)]
    # Determine how many tasks to run in parallel.
    n_tasks, _ = get_num_parallel(len(reps), max_procs,
                                  parallel=True, hybrid=False)
    if n_tasks > 1:
        futures: list[Future] = list()
        with ProcessPoolExecutor(max_workers=n_tasks) as pool:
            for rep in reps:
                futures.append(pool.submit(rep.run))
        runs = [future.result() for future in futures]
    else:
        runs = [rep.run() for rep in reps]
    logger.info(f"Ended EM with {n_clusters} clusters ({n_runs} runs)")
    # Sort the replicate runs of EM clustering in ascending order
    # by BIC so that the run with the best (smallest) BIC is first.
    return sorted(runs, key=lambda x: x.bic)
