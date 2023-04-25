from concurrent.futures import Future, ProcessPoolExecutor
from logging import getLogger
from pathlib import Path

from .emalgo import EmClustering
from .bvec import BitVector
from .write import write_results
from ..util.parallel import get_num_parallel
from ..vector.load import VectorLoader

logger = getLogger(__name__)


def cluster_sect(loader: VectorLoader, end5: int, end3: int, *,
                 include_gu: bool, include_del: bool, include_ins: bool,
                 max_polya: int, max_muts_per_read: int, min_mut_dist: int,
                 min_reads: int, info_thresh: float, signal_thresh: float,
                 max_clusters: int, n_runs: int, conv_thresh: float,
                 min_iter: int, max_iter: int, max_procs: int, out_dir: Path):
    """ Cluster a section of a set of bit vectors. """
    # Generate the bit vectors from the mutation vector report.
    try:
        bitvectors = BitVector(loader, end5=end5, end3=end3,
                               include_gu=include_gu,
                               include_del=include_del,
                               include_ins=include_ins,
                               max_polya=max_polya,
                               max_muts_per_read=max_muts_per_read,
                               min_mut_dist=min_mut_dist,
                               min_reads=min_reads,
                               info_thresh=info_thresh,
                               signal_thresh=signal_thresh)
        bitvectors.publish_preprocessing_report(out_dir)
    except Exception as error:
        logger.critical(f"Failed to generate bit vectors: {error}")
        return
    # Cluster the bit vectors.
    logger.info(f"Began clustering {bitvectors} with up to k = {max_clusters} "
                f"clusters and {n_runs} replicate runs per number of clusters")
    try:
        clusters = run_max_clust(bitvectors, max_clusters, n_runs,
                                 min_iter=min_iter, max_iter=max_iter,
                                 conv_thresh=conv_thresh, max_procs=max_procs)
    except Exception as error:
        logger.critical(f"Failed to cluster {bitvectors}: {error}")
        return
    logger.info(
        f"Ended clustering {bitvectors} and obtained {len(clusters)} clusters")
    # Output the results of clustering.
    try:
        resps_file = write_results(bitvectors, clusters, out_dir)
    except Exception as error:
        logger.critical(f"Failed to write clusters for {bitvectors}: {error}")
        return
    # Return the path to the file of read responsibilities.
    return resps_file


def run_max_clust(bitvectors: BitVector, max_clusters: int, n_runs: int, *,
                  conv_thresh: float, min_iter: int, max_iter: int,
                  max_procs: int):
    """
    """
    if n_runs < 1:
        logger.warning(f"Number of EM runs must be ≥ 1: setting to 1")
        n_runs = 1
    logger.info(f"Began clustering {bitvectors} with up to k = {max_clusters} "
                f"clusters and {n_runs} replicate runs per number of clusters")
    # Store the best clustering replicate for each number of clusters.
    # The 0th index is filled with None as a placeholder.
    best_reps: dict[int, EmClustering | None] = {0: None}
    # Run clustering for each number of clusters, up to max_clusters.
    while len(best_reps) <= max_clusters:
        # Run EM clustering n_runs times with different starting points.
        reps_n = run_n_clust(bitvectors, len(best_reps),
                             n_runs=(n_runs if len(best_reps) > 1 else 1),
                             conv_thresh=(conv_thresh if len(best_reps) > 1
                                          else float("inf")),
                             min_iter=(min_iter if len(best_reps) > 1 else 1),
                             max_iter=(max_iter if len(best_reps) > 1 else 2),
                             max_procs=max_procs)
        # Find the best (smallest) BIC obtained from clustering.
        best_rep = reps_n[0]
        logger.debug(
            f"The best BIC with {len(best_reps)} cluster(s) is {best_rep.bic}")
        if (best_rep_prev := best_reps[len(best_reps) - 1]) is not None:
            # Check if the best BIC is worse (larger) than the best BIC
            # from clustering with one less cluster.
            if best_rep.bic > best_rep_prev.bic:
                # If the best BIC is larger, then the model with one
                # less cluster was better than the current model. Stop.
                logger.info(f"The BIC increased from {best_rep_prev.bic} "
                            f"(k = {len(best_reps) - 1}) to {best_rep.bic} "
                            f"(k = {len(best_reps)}): discarding results with "
                            f"{len(best_reps)} clusters and stopping")
                break
            else:
                logger.debug(f"The BIC decreased from {best_rep_prev.bic} "
                             f"(k = {len(best_reps) - 1}) to {best_rep.bic} "
                             f"(k = {len(best_reps)}): keeping results with "
                             f"{len(best_reps)} clusters")
        # Add the best replicate with this number of clusters to the
        # list of best replicates with all numbers of clusters.
        best_reps[len(best_reps)] = best_rep
    # Remove the 0th index placeholder.
    best_reps.pop(0)
    return best_reps


def run_n_clust(bitvector: BitVector, n_clusters: int, n_runs: int, *,
                conv_thresh: float, min_iter: int, max_iter: int,
                max_procs: int) -> list[EmClustering]:
    logger.info(f"Began EM with {n_clusters} clusters ({n_runs} runs)")
    # Initialize one EmClustering object for each replicate run.
    reps = [EmClustering(bitvector, n_clusters, rep,
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
    # Sort the replicate runs of EM clustering in ascending order
    # by BIC so that the run with the best (smallest) BIC is first.
    sorted_runs = sorted(runs, key=lambda x: x.bic)
    logger.info(f"Ended EM with {n_clusters} clusters ({n_runs} runs)")
    return sorted_runs
