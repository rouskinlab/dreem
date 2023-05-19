from logging import getLogger
from typing import Iterable

import numpy as np
import pandas as pd

from dreem.bit.call import BitCaller
from dreem.util.sect import sects_to_pos, Section
from dreem.mut.load import VectorLoader

logger = getLogger(__name__)


def sum_bits(loader: VectorLoader,
             sections: Iterable[Section] = (), *,
             by_pos: dict[str, BitCaller] | None = None,
             by_vec: dict[str, BitCaller] | None = None,
             ) -> dict[str, tuple[dict, dict]]:
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
    by_pos: dict[str, BitCaller] | None = None
        Bit callers to use for counting the matching bits at each
        position of each section.
    by_vec: dict[str, BitCaller] | None = None
        Bit callers to use for counting the matching bits within each
        section in each mutation vector.

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
            {key: pd.Series(dtype=float) for key in by_vec},
            # Axis 1, columns (positions): Initialize for each query a
            # Series of bit counts for each position.
            {key: pd.Series(np.zeros(section.length, dtype=float),
                            index=section.columns)
             for key in by_pos},
        )
    # Iterate over all the batches of mutation vectors.
    for muts in loader.iter_batches(sects_to_pos(sections)):
        print("MUTS")
        print(muts)
        # Iterate over all bit callers.
        for key, bit_caller in (by_vec | by_pos).items():
            print("KEY")
            print(key)
            # Call all bits in this batch of mutation vectors.
            bits = bit_caller.call(muts)
            if key in by_vec:
                # Count the bits within each section of each vector in
                # this batch, then append to the previous batches.
                for section in sections:
                    counts[section.name][0][key] = pd.concat([
                        counts[section.name][0][key],
                        bits.loc[:, section.columns].sum(axis=1)])
            if key in by_pos:
                # Count the bits in this batch at each position.
                bsum = bits.sum(axis=0)
                # Add the bit count to each section's count.
                for section in sections:
                    counts[section.name][1][key] += bsum.loc[section.columns]
    return counts
