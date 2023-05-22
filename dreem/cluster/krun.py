from concurrent.futures import Future, ProcessPoolExecutor
from logging import getLogger
from math import inf

from .emalgo import EmClustering
from .metric import find_best_k
from .report import ClusterReport
from .write import write_results
from ..filt.load import FilterLoader
from ..util.bit import UniqMutBits
from ..util.parallel import get_num_parallel

logger = getLogger(__name__)


def cluster_filt(loader: FilterLoader, max_clusters: int, n_runs: int, *,
                 min_iter: int, max_iter: int, conv_thresh: float,
                 max_procs: int):
    """ Run all processes of clustering reads from one filter. """
    logger.info(f"Began clustering {loader} with up to k = {max_clusters} "
                f"clusters and {n_runs} replicate runs per number of clusters")
    # Get the unique bit vectors.
    uniq_muts = loader.get_bit_vectors().get_unique_muts()
    try:
        clusts = run_max_clust(loader, uniq_muts,
                               max_clusters, n_runs,
                               min_iter=min_iter,
                               max_iter=max_iter,
                               conv_thresh=conv_thresh,
                               max_procs=max_procs)
    except Exception as error:
        raise
        logger.critical(f"Failed to cluster {bvec}: {error}")
        return
    logger.info(f"Ended clustering {loader}: {find_best_k(clusts)} clusters")
    # Output the results of clustering.
    try:
        report = ClusterReport.from_clusters(loader, uniq_muts, clusts,
                                             max_clusters, n_runs,
                                             min_iter=min_iter,
                                             max_iter=max_iter,
                                             conv_thresh=conv_thresh)
        report_file = report.save()
        write_results(loader, clusts)
    except Exception as error:
        raise
        logger.critical(f"Failed to write clusters for {bvec}: {error}")
        return
    # Return the path of the clustering report file.
    return report_file


def run_max_clust(loader: FilterLoader, uniq_muts: UniqMutBits,
                  max_clusters: int, n_runs: int, *,
                  min_iter: int, max_iter: int, conv_thresh: float,
                  max_procs: int):
    """
    Find the optimal number of clusters for EM, up to max_clusters.
    """
    if n_runs < 1:
        logger.warning(f"Number of EM runs must be ≥ 1: setting to 1")
        n_runs = 1
    logger.info(f"Began clustering {loader} with up to {max_clusters} clusters "
                f"and {n_runs} runs per number of clusters")
    # Store the clustering runs for each number of clusters.
    runs: dict[int, list[EmClustering]] = dict()
    # Run clustering for each number of clusters, up to max_clusters.
    while len(runs) < max_clusters:
        # Get the number of clusters, k.
        k = len(runs) + 1
        # Run EM clustering n_runs times with different starting points.
        runs[k] = run_n_clust(loader, uniq_muts, k,
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


def run_n_clust(loader: FilterLoader,
                uniq_muts: UniqMutBits,
                n_clusters: int,
                n_runs: int, *,
                min_iter: int, max_iter: int, conv_thresh: float,
                max_procs: int) -> list[EmClustering]:
    """ Run EM with a specific number of clusters. """
    logger.info(f"Began EM with {n_clusters} clusters ({n_runs} runs)")
    # Initialize one EmClustering object for each replicate run.
    reps = [EmClustering(loader, uniq_muts, n_clusters, rep,
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
