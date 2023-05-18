import numpy as np
import pandas as pd


def get_mut_info_bits(muts: pd.DataFrame, refs: pd.DataFrame | None = None):
    """ Ensure muts and refs are both boolean and have the same axes,
    determine which bits are informative, mask uninformative bits in
    muts and refs to zero, and return the validated data frames of
    mutations, matches, and informative bits (all with the same shape).
    If refs is not given, initialize it to the logical not of muts. """
    # Ensure each data frame is boolean.
    mut_bits = muts.astype(bool, copy=True)
    if refs is None:
        # If no refs data frame was given, assume that it is the logical
        # not of muts, so all bits are informative.
        info_bits = pd.DataFrame(np.ones_like(mut_bits, dtype=bool),
                                 index=mut_bits.index,
                                 columns=mut_bits.columns)
    else:
        # Ensure the indexes and columns of muts and refs match.
        if not muts.index.equals(refs.index):
            raise ValueError(f"Got different indexes for muts {muts.index} "
                             f"and refs {refs.index}")
        if not muts.columns.equals(refs.columns):
            raise ValueError(f"Got different columns for muts"
                             f"{mut_bits.columns} and refs {refs.columns}")
        # Determine which bits are informative.
        info_bits = mut_bits ^ refs.astype(bool, copy=False)
        # Mask any uninformative bits in mut_bits to zero.
        mut_bits &= info_bits
    # Return boolean data frames of the mutated and informative bits.
    return mut_bits, info_bits
