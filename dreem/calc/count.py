from logging import getLogger
from sys import byteorder
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from dreem.vector.load import VectorLoader
from dreem.util.sect import sections_to_positions


logger = getLogger(__name__)

QEQ = "="
QSB = "<"
QSP = ">"
QEB = "{"
QEP = "}"
QMETHOD = QSB, QEB, QEQ, QEP, QSP


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
