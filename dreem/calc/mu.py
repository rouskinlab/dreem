from functools import partial

import numpy as np
import pandas as pd
from scipy.optimize import newton_krylov

MIN_MUT_DIST = 4


def cluster_mus(refbits: pd.DataFrame,
                mutbits: pd.DataFrame,
                members: pd.DataFrame,
                min_dist: int = MIN_MUT_DIST) -> pd.DataFrame:
    """
    Calculate the mutation rate at each position in a mutational profile
    for one or more clusters.

    Parameters
    ----------
    refbits: DataDrame
        Reference bits: each index (i) is the name of a read, each
        column (j) a position in the profile, and each value (i, j) a
        boolean where True means that position (j) of read (i) matches
        the reference sequence and False means that it does not.
    mutbits: DataDrame
        Mutated bits: each index (i) is the name of a read, each column
        (j) a position in the profile, and each value (i, j) a boolean
        where True means that position (j) of read (i) is mutated and
        False means that it is not.
    members: DataFrame
        Cluster membership: each index (i) is the name of a read, each
        column (k) the name of a cluster, and each value (i, k) the
        likelihood that read (i) came from cluster (k). Any read absent
        from members is ignored when calculating mutations.
    min_dist: int
        Minimum distance between mutations. This parameter is used to
        correct the bias in observed mutation rates that is caused by
        dropout of reads with nearby mutations, but NOT to remove any
        rows of mutbits or refbits with nearby mutations. Instead, the
        'members' parameter is responsible for filtering such vectors,
        because it comes from the clustering step in which vectors with
        nearby mutations were already removed.

    Returns
    -------
    DataFrame
        Mutation rates: each index (j) is a position in the profile,
        each column (k) the name of a cluster, and each value (j, k)
        the mutation fraction at position (j) in cluster (k).
    """
    # Ensure dimensions and axes match.
    if refbits.shape != mutbits.shape:
        raise ValueError("refbits and mutbits must have the same shape, "
                         f"but got {refbits.shape} ≠ {mutbits.shape}")
    if np.any(np.not_equal(refbits.index, mutbits.index)):
        raise ValueError("refbits and mutbits do not have the same read names")
    if np.any(np.not_equal(refbits.columns, mutbits.columns)):
        raise ValueError("refbits and mutbits do not have the same positions")
    # Weighted count of matched bits at each position in each cluster.
    refsums = refbits.loc[members.index].T.dot(members)
    # Weighted count of mutated bits at each position in each cluster.
    mutsums = mutbits.loc[members.index].T.dot(members)
    # The observed cluster mutation rate of each position is the number
    # of mutated bits divided by the number of matched or mutated bits.
    mus = mutsums / (refsums + mutsums)
    if min_dist >= 2:
        # Apply the bias correction to each cluster (column).
        mus = mus.apply(partial(bv_means_to_mu, min_dist))
    return mus


def _calc_denom_s2(i: int,
                   mu: np.ndarray,
                   denom_probs: np.ndarray,
                   s2_probs: np.ndarray,
                   min_dist: int):
    """
    calc_denom from DREEM
    """
    if i >= mu.size:
        # Base case.
        return 1.0, s2_probs
    if denom_probs[i] >= 0.0:
        # Already encountered.
        return denom_probs[i], s2_probs
    # Probability that the next base (if any) is not mutated.
    s1 = _calc_denom_s2(i + 1, mu, denom_probs, s2_probs, min_dist)[0]
    # Probability that
    s2 = ((1.0 - mu[i + 1: i + min_dist]).prod() *
          _calc_denom_s2(i + min_dist, mu, denom_probs, s2_probs, min_dist)[0])
    denom_probs[i] = ((1.0 - mu[i]) * s1) + (mu[i] * s2)
    s2_probs[i] = s2
    return denom_probs[i], s2_probs


def calc_denom_s2(mu: np.ndarray | pd.Series, min_dist: int):
    if isinstance(mu, np.ndarray):
        if len(mu.shape) != 1:
            raise ValueError("Expected mu to be a 1-dimensional array, "
                             f"but got a {len(mu.shape)}-dimensional array")
        if mu.dtype is not float:
            raise ValueError("Expected mu to be an array of type 'float', "
                             f"but got an array of type '{mu.dtype.__name__}'")
        return _calc_denom_s2(0, mu,
                              np.full_like(mu, -1.0, dtype=float),
                              np.full_like(mu, -1.0, dtype=float),
                              min_dist)
    if isinstance(mu, pd.Series):
        return pd.Series(calc_denom_s2(mu.values, min_dist), index=mu.index)
    raise TypeError("Expected mu to be of type 'np.ndarray' or 'pd.Series', "
                    f"but got '{type(mu).__name__}'")


def calc_denom(mu: np.ndarray | pd.Series, min_dist: int):
    return calc_denom_s2(mu, min_dist)[0]


def get_mu_der(mu, bv_means, min_dist: int):
    """
    :param mu: (corrected) mutation rates
    :param bv_means: uncorrected bv means
    :param min_dist: minimum distance between mutations
    :return:
    """
    return mu_to_bv_means(mu, min_dist) - bv_means


def bv_means_to_mu(bv_means, min_dist: int):
    return newton_krylov(lambda mu_k: get_mu_der(mu_k, bv_means, min_dist),
                         bv_means)


def mu_to_bv_means(mu, min_dist: int):
    denom, s2 = calc_denom_s2(mu, min_dist)
    denom_rev, s2_rev = calc_denom_s2(mu[::-1], min_dist)
    return np.array([
        (mu[i] * s2[i] * s2_rev[len(mu) - i - 1] / denom)
        for i in range(len(mu))
    ])
