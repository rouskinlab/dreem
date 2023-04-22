from concurrent.futures import Future, ProcessPoolExecutor

from .em import EMclustering
from .bitvector import BitVector


def run(bitvector: BitVector, max_clusters: int, num_runs: int, min_iter: int,
        max_iter: int, convergence_cutoff: float, max_procs: int):
    """
    Run the clustering algorithm many times to iterate over K_max the number of clusters and num_runs the number of runs.

    Parameters
    ----------

    bitvector: object
        bitvector object. Conatins the processed bitvector and the count of reads per bitvector.

    max_clusters: int
        Number of clusters. The algorithm will iterate from 1 to max_clusters.

    num_runs: int
        Number of runs. The algorithm will be launched num_runs times.

    Returns
    -------

    The output is a dictionary containing the results of EMclustering for each K and each run.
    The first level of the dictionary is the number of clusters.
    The second level of the dictionary is a list containing the results of EMclustering for each run.
    The second levels of the dictionary are sorted by performance (log-likelihood).

    Example:
    --------

    results = {
        K1: [EMclustering1, EMclustering2, ...],
        K2: [EMclustering1, EMclustering2, ...],
        ...
    }

    The EMclustering output dictionary contains the following attributes:
        - mu: array (K x D)
            Mean of each cluster.
        - pi: array (K)
            Probability of each cluster.
        - log_likelihood: float
            Log-likelihood of the model.
    """

    # global dT # !! For test_input !!
    em = EMclustering(bitvector, 1, 10, 10, float("inf"))
    results = {1: [em.run()]}
    for k in range(2, max_clusters + 1):
        em = EMclustering(bitvector, k, min_iter, max_iter, convergence_cutoff)
        futures: list[Future] = list()
        with ProcessPoolExecutor(max_workers=min(max_procs, num_runs)) as pool:
            futures.append(pool.submit(em.run))
        results[k] = sorted([future.result() for future in futures],
                            key=lambda res: res['log_likelihood'], reverse=True)
    return results
