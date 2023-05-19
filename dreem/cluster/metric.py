from itertools import chain, combinations

import numpy as np
import pandas as pd

from .emalgo import EmClustering


def get_common_clusts(clusters: dict[int, list[EmClustering]]
                      ) -> list[EmClustering]:
    """ Return the EM runs with all attributes in common. If all runs
    have just one cluster, then return those runs. Otherwise, return
    all runs except those with one cluster. """
    n = len(clusters)
    return list(chain(*[clusters[n] for n in range(min(n, 2), n + 1)]))


def calc_bics(clusters: dict[int, list[EmClustering]]) -> dict[int, float]:
    """ For each number of clusters, return the best (smallest) BIC. """
    return {k: min(run.bic for run in runs) for k, runs in clusters.items()}


def find_best_k(clusters: dict[int, list[EmClustering]]) -> int:
    """ Find the number of clusters with the best (smallest) BIC. """
    return sorted(calc_bics(clusters).items(), key=lambda k_bic: k_bic[1])[0][0]


def get_converged(clusters: dict[int, list[EmClustering]]
                  ) -> dict[int, [list[int]]]:
    """ For each number of clusters, return a list of the number of
    iterations required for each run to converge, or 0 for each run
    that did not converge within the maximum number of iterations. """
    return {k: [run.iter if run.converged else 0 for run in runs]
            for k, runs in clusters.items()}


def get_log_likes(clusters: dict[int, list[EmClustering]]
                  ) -> dict[int, [list[float]]]:
    """ For each number of clusters, return a list of the log likelihood
    of each EM run. """
    return {k: [run.log_like for run in runs] for k, runs in clusters.items()}


def get_log_like_mean(log_likes: dict[int, [list[float]]]) -> dict[int, float]:
    """ For each number of clusters, find the mean log likelihood. """
    return {k: np.mean(log_like_k) for k, log_like_k in log_likes.items()}


def get_log_like_std(log_likes: dict[int, [list[float]]]) -> dict[int, float]:
    """ For each number of clusters, find the stdev log likelihood. """
    return {k: np.std(log_like_k) for k, log_like_k in log_likes.items()}


def get_var_info(clusters: dict[int, list[EmClustering]]) -> dict[int, float]:
    """ For each number of clusters, find variation of information. """
    return {k: calc_exp_var_info(runs) for k, runs in clusters.items()}


def calc_var_info(p: np.ndarray, q: np.ndarray, r: np.ndarray):
    """
    Calculate the variation of information for two partitions, X and Y,
    of the same set A. For more details and the source of the formula,
    see https://en.wikipedia.org/wiki/Variation_of_information.

    Parameters
    ----------
    p: ndarray
        A 1-dimensional array of the fraction of the total elements in
        each partition of X. All must be in (0, 1] and sum to 1.
    q: ndarray
        A 1-dimensional array of the fraction of the total elements in
        each partition of Y. All must be in (0, 1] and sum to 1.
    r: ndarray
        A 2-dimensional array (len(p) x len(q)) of the fraction of the
        total elements that are in each pair of partitions from X and Y.
        All must be in (0, 1] and sum to 1.

    Returns
    -------
    float
        Variation of information
    """
    # Verify dimensions.
    if p.ndim != 1:
        raise ValueError(f"p must be 1-dimensional, but got {p.ndim} dims")
    if q.ndim != 1:
        raise ValueError(f"q must be 1-dimensional, but got {q.ndim} dims")
    if r.ndim != 2:
        raise ValueError(f"r must be 2-dimensional, but got {r.ndim} dim(s)")
    # Verify bounds.
    if np.any(p <= 0.) or np.any(p > 1.):
        raise ValueError(f"All values in p must be in (0, 1], but got {p}")
    if np.any(q <= 0.) or np.any(q > 1.):
        raise ValueError(f"All values in q must be in (0, 1], but got {q}")
    if np.any(r <= 0.) or np.any(r > 1.):
        raise ValueError(f"All values in r must be in (0, 1], but got {r}")
    # Verify sums.
    if not np.isclose(p.sum(), 1.):
        raise ValueError(f"p must sum to 1, but got {p.sum()}")
    if not np.isclose(q.sum(), 1.):
        raise ValueError(f"q must sum to 1, but got {q.sum()}")
    if not np.isclose(r.sum(), 1.):
        raise ValueError(f"r must sum to 1, but got {r.sum()}")
    # Compute the variation of information.
    return np.sum(r * (np.log(p)[:, np.newaxis] + np.log(q)[np.newaxis, :]
                       - 2 * np.log(r)))


def calc_var_info_runs(run1: pd.DataFrame, run2: pd.DataFrame):
    """ Calculate the variation of information between two EM runs. """
    # Find the number of reads and clusters.
    n_reads, n_clusts = run1.shape
    if run2.shape != (n_reads, n_clusts):
        raise ValueError(f"Dimensions of run1 {run1.shape} "
                         f"and run2 {run2.shape} differ")
    if n_clusts <= 0:
        raise ValueError(f"Number of clusters must be ≥ 1, but got {n_clusts}")
    if n_reads == 0:
        # There is zero variation if there are zero reads.
        return 0.
    # Verify that run1 and run2 have the same reads (indexes) and
    # clusters (columns).
    if not run1.index.equals(run2.index):
        raise ValueError("Indexes of resps1 and resps2 differ")
    if not run1.columns.equals(run2.columns):
        raise ValueError("Clusters of run1 and run2 differ")
    # Verify that the probability that each read belongs to any cluster
    # equals 1.
    if not np.allclose(run1.sum(axis=1), 1.):
        raise ValueError("Probabilities of reads in run 1 did not all sum to 1")
    if not np.allclose(run2.sum(axis=1), 1.):
        raise ValueError("Probabilities of reads in run 2 did not all sum to 1")
    # For each run, compute the proportion of reads in each cluster.
    props1 = run1.mean(axis=0).values
    props2 = run2.mean(axis=0).values
    # For each pair of clusters, compute the proportion of reads in both
    # clusters.
    props12 = np.array([[np.vdot(run1[c1], run2[c2]) / n_reads
                         for c2 in run2.columns]
                        for c1 in run1.columns])
    # Compute the variation of information.
    return calc_var_info(props1, props2, props12)


def calc_exp_var_info(runs: list[EmClustering]):
    """ Calculate the expected variation of information among ≥ 2 runs
    of EM clustering. """
    # List every pair of EM runs.
    pairs = list(combinations(range(len(runs)), 2))
    if not pairs:
        # Variation of information defaults to 0 if no pairs exist.
        return 0.
    # Find the mean variation of information among pairs of EM runs.
    # FIXME: This part can be made more computationally efficient by caching each resps and each props1 and props2.
    return sum(calc_var_info_runs(runs[i1].output_resps(),
                                  runs[i2].output_resps())
               for i1, i2 in pairs) / len(pairs)
