from logging import getLogger
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from dreem.bit.call import muts_to_bits
from dreem.util.sect import sects_to_pos, Section
from dreem.mut.load import VectorLoader

logger = getLogger(__name__)


def sum_bits(loader: VectorLoader,
             sections: Iterable[Section] = (), *,
             by_pos: Sequence[tuple[int, str]] = (),
             by_vec: Sequence[tuple[int, str]] = (),
             numeric: bool = False) -> dict[str, tuple[dict, dict]]:
    """
    For each section, count the mutations that agree with each query by
    position and each query by vector.

    Parameters
    ----------
    loader: VectorLoader
        VectorLoader from which to load the mutation vectors
    sections: Iterable[tuple[int, int]] = ()
        Iterable of 2-tuples, each defining the 5' and 3' coordinates
        of one section over which to count the bits. If empty, then use
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
    # If no sections were given, then use the entire reference.
    if not sections:
        sections = [loader.section()]
    # For each pair of coordinates, initialize counts of the bits by
    # vector (row, axis 0) and position (column, axis 1).
    counts = dict()
    for section in sections:
        counts[section.name] = (
            # Axis 0, rows (vectors): Initialize for each query a Series
            # that will be filled with the bit count for each vector.
            {qryrel: pd.Series(dtype=float) for qryrel in by_vec},
            # Axis 1, columns (positions): Initialize for each query a
            # Series of bit counts for each position.
            {qryrel: pd.Series(np.zeros(section.length, dtype=float),
                               index=(section.positions if numeric
                                      else section.columns))
             for qryrel in by_pos}
        )
    # Iterate over all the batches.
    queries_rels = set(by_pos) | set(by_vec)
    for batch in loader.iter_batches(sects_to_pos(sections), numeric=numeric):
        # Iterate over all queries and relationships.
        for qryrel in queries_rels:
            # Compute all bits in this batch.
            bits = muts_to_bits(batch, *qryrel)
            if qryrel in by_vec:
                # Count the bits within each section of each vector in
                # this batch, then append to the previous batches.
                for section in sections:
                    counts[section.name][0][qryrel] = pd.concat([
                        counts[section.name][0][qryrel],
                        bits.loc[:, (section.positions if numeric
                                     else section.columns)].sum(axis=1)])
            if qryrel in by_pos:
                # Count the bits in this batch at each position.
                bits_sum = bits.sum(axis=0)
                # Add the bit count to each section's count.
                for section in sections:
                    counts[section.name][1][qryrel] += bits_sum.loc[
                        section.positions if numeric else section.columns]
    return counts
