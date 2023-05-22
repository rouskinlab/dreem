import pandas as pd

from .bias import obs_to_real


def mus_obs_to_real(mus_obs: pd.DataFrame, min_gap: int):
    """
    Parameters
    ----------
    mus_obs: DataFrame
        A (positions x clusters) array of the observed mutation rate at
        each position in each cluster.
    min_gap: int
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
        A (positions x clusters) array of the real mutation rate at each
        position in each cluster.
    """
    # Estimate the real mutation rates.
    return pd.DataFrame(obs_to_real(mus_obs.values, min_gap),
                        index=mus_obs.index, columns=mus_obs.columns)


def reads_to_mus_obs(members: pd.DataFrame,
                     muts: pd.DataFrame,
                     refs: pd.DataFrame | None = None) -> pd.DataFrame:
    """
    Calculate the mutation rate at each position in a mutational profile
    for one or more clusters.

    Parameters
    ----------
    members: DataFrame
        Cluster membership: each index (i) is the name of a read, each
        column (k) the name of a cluster, and each value (i, k) the
        likelihood that read (i) came from cluster (k). Any read absent
        from members is ignored when calculating mutations.
    muts: DataDrame
        Mutated bits: each index (i) is the name of a read, each column
        (j) a position in the profile, and each value (i, j) a boolean
        where True means that position (j) of read (i) is mutated and
        False means that it is not.
    refs: DataDrame | None = None
        Reference bits: each index (i) is the name of a read, each
        column (j) a position in the profile, and each value (i, j) a
        boolean where True means that position (j) of read (i) matches
        the reference sequence and False means that it does not.
        If omitted, defaults to the logical not of ```mutbits```.

    Returns
    -------
    DataFrame
        Mutation rates: each index (j) is a position in the profile,
        each column (k) the name of a cluster, and each value (j, k)
        the mutation fraction at position (j) in cluster (k).
    """
    # Ensure dimensions and axes match.
    muts, info = get_mut_info_bits(muts, refs)
    # The observed cluster mutation rate of each position is the number
    # of mutated bits divided by the number of informative bits.
    return (muts.loc[members.index].T.dot(members) /
            info.loc[members.index].T.dot(members))


def reads_to_mus_real(members: pd.DataFrame, min_gap: int,
                      muts: pd.DataFrame, refs: pd.DataFrame | None = None):
    """ A convenience function to chain ```reads_to_mus_obs``` and
    ```mus_obs_to_real```. """
    return mus_obs_to_real(reads_to_mus_obs(members, muts, refs), min_gap)
