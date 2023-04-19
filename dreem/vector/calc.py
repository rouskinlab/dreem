from functools import partial
from logging import getLogger
from sys import byteorder
from typing import Iterable, Sequence

import numpy as np
import pandas as pd
from scipy.optimize import newton_krylov

from .load import VectorLoader
from ..util.sect import sections_to_positions


logger = getLogger(__name__)

QEQ = "="
QSB = "<"
QSP = ">"
QEB = "{"
QEP = "}"
QMETHOD = QSB, QEB, QEQ, QEP, QSP

MIN_MUT_DIST = 4


def get_bits(vectors: pd.DataFrame,
             query: int,
             rel: str = QEQ) -> pd.DataFrame:
    """
    Return a boolean array of the same shape as vectors where element
    i,j is True if and only if the byte at element i,j of vectors meets
    the requested relationship between it and the query byte

    Parameters
    ----------
    vectors: DataFrame
        Mutation vectors
    query: int
        Byte to query; must be in range 0 - 255
    rel: str = "equals"
        Method to decide whether a byte counts, as follows:
        - "=": count only bytes equal to query
        - "<": count strict bitwise subsets of query
        - ">": count strict bitwise supersets of query
        - "{": count bitwise subsets of and bytes equal to query
        - "}": count bitwise supersets of and bytes equal to query

    Returns
    -------
    DataFrame
        Boolean type DataFrame of the same shape as vectors where each
        element is True if and only if the element at the same position
        in vectors fulfilled the relationship with query
    """
    # Validate the query byte.
    if not isinstance(query, int):
        raise TypeError(
            f"Expected query of type int, but got type {type(query).__name__}")
    try:
        query.to_bytes(1, byteorder)
    except OverflowError:
        raise ValueError(f"Expected query in range 0 - 255, but got {query}")
    if rel == QEQ:
        return np.equal(vectors, query)
    if rel == QEB:
        return np.equal(np.bitwise_or(vectors, query), query)
    if rel == QEP:
        return np.equal(np.bitwise_and(vectors, query), query)
    if rel == QSB:
        return np.logical_and(get_bits(vectors, query, QEB),
                              np.not_equal(vectors, query))
    if rel == QSP:
        return np.logical_and(get_bits(vectors, query, QEP),
                              np.not_equal(vectors, query))
    raise ValueError(f"Parameter 'rel' must be in {QMETHOD}, but got '{rel}'")


def sum_bits(loader: VectorLoader, *,
             coords: Iterable[tuple[int, int]] = (),
             by_pos: Sequence[tuple[int, str]] = (),
             by_vec: Sequence[tuple[int, str]] = (),
             numeric: bool = False) -> dict[tuple[int, int], tuple[dict, dict]]:
    """
    Return the number of mutations that agree with the query for each
    position and/or vector.

    Parameters
    ----------
    loader: VectorLoader
        VectorLoader from which to load the mutation vectors
    coords: Iterable[tuple[int, int]] = ()
        Iterable of 2-tuples, each defining the 5' and 3' coordinates
        of one section over which to count the bits. If None, then use
        the entire reference as the section.
    by_pos: Sequence[tuple[int, str]] = ()
        Queries and relationships to use for counting the matching bits
        at each position of each section (also see get_bits).
    by_vec: Sequence[tuple[int, str]] = ()
        Queries and relationships to use for counting the matching bits
        within each section in each mutation vector (also see get_bits).
    numeric: bool = False
        Whether to convert the columns from base-position strings to
        numeric (specifically, integer) values of the positions, e.g.
        ['G1', 'T2', ...] if False, [1, 2, ...] if True

    Returns
    -------
    dict[tuple[int, int], tuple[dict, dict]]
        Dictionary mapping each (5', 3') coordinate pair to a tuple of
        two dictionaries of counts. The first dictionary (item 0) maps
        each key of by_vec to a Series of counts per vector (axis 0),
        indexed by the read name. The second dictionary (item 1) maps
        the keys of by_pos to a Series of counts per position (axis 1);
        its index will be the positions as integers if numeric is True,
        otherwise strings indicating the base and position.
    """
    # If no coordinates were given, then use the entire reference.
    if not coords:
        coords = [(1, len(loader.seq))]
    # For each pair of coordinates, initialize counts of the bits by
    # vector (row, axis 0) and position (column, axis 1).
    sections = list()
    counts = dict()
    for end5, end3 in coords:
        try:
            sect = loader.section(end5, end3)
        except Exception as error:
            logger.error(f"Invalid section {end5}-{end3} of {loader}: {error}")
            continue
        # Add section to list of valid sections.
        sections.append(sect)
        counts[sect.coord] = (
            # Axis 0, rows (vectors): Initialize for each query a Series
            # that will be filled with the bit count for each vector.
            {qryrel: pd.Series(dtype=int) for qryrel in by_vec},
            # Axis 1, columns (positions): Initialize for each query a
            # Series of bit counts for each position.
            {qryrel: pd.Series(np.zeros(sect.length, dtype=int),
                               index=(sect.positions if numeric
                                      else sect.columns))
             for qryrel in by_pos}
        )
    # Iterate over all the batches.
    queries_rels = set(by_pos) | set(by_vec)
    for batch in loader.iter_batches(sections_to_positions(sections),
                                     numeric=numeric):
        # Iterate over all queries and relationships.
        for qryrel in queries_rels:
            # Compute all bits in this batch.
            bits = get_bits(batch, *qryrel)
            if qryrel in by_vec:
                # Count the bits within each section of each vector in
                # this batch, then append to the previous batches.
                for sect in sections:
                    counts[sect.coord][0][qryrel] = pd.concat([
                        counts[sect.coord][0][qryrel],
                        bits.loc[:, (sect.positions if numeric
                                     else sect.columns)].sum(axis=1)])
            if qryrel in by_pos:
                # Count the bits in this batch at each position.
                bits_sum = bits.sum(axis=0)
                # Add the bit count to each section's count.
                for sect in sections:
                    counts[sect.coord][1][qryrel] += bits_sum.loc[
                        sect.positions if numeric else sect.columns]
    return counts


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
