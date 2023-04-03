import os.path

import numpy as np
import pandas as pd

from ..util.seq import (EVERY, MATCH, DELET, INS_5, INS_3,
                        SUB_A, SUB_C, SUB_G, SUB_T, SUB_N)

from ..vector.profile import VectorReader


def generate_mut_profile_from_bit_vector(vectors: VectorReader,
                                         clustering_file: str = "",
                                         verbose=False):
    """
    Generate a mutation profile from mutation vectors.

    Parameters
    ----------
    vectors: VectorReader
        Vector reader object

    Returns
    -------
    df_row: dict
        - info: the count per residue of bases for which we have information.
        - cov: the count per residue of bases for which we have coverage.
        - sub_N: the count per residue of bases for which we have mutations.
        - sub_A: the count per residue of bases who are modified to A.
        - sub_C: the count per residue of bases who are modified to C.
        - sub_G: the count per residue of bases who are modified to G.
        - sub_T: the count per residue of bases who are modified to T.
        - sub_rate: the Mutation fraction per residue.
        - num_aligned: the number of aligned reads.
    
    """
    # Convert to a mutation profile
    out = dict()
    out['sequence'] = vectors.seq.decode()
    out['num_aligned'] = vectors.n_vectors
    out["match"] = vectors.count_muts_by_pos(MATCH | INS_5, subsets=True)
    out["sub_A"] = vectors.count_muts_by_pos(SUB_A)
    out["sub_C"] = vectors.count_muts_by_pos(SUB_C)
    out["sub_G"] = vectors.count_muts_by_pos(SUB_G)
    out["sub_T"] = vectors.count_muts_by_pos(SUB_T)
    out["sub_N"] = vectors.count_muts_by_pos(SUB_N, subsets=True)
    out["del"] = vectors.count_muts_by_pos(DELET)
    out["ins"] = vectors.count_muts_by_pos(INS_3, supersets=True)
    # Count non-blank bytes
    out["cov"] = vectors.count_muts_by_pos(EVERY, subsets=True)
    # Unambiguously matching or substituted (informative)
    out["info"] = out["match"] + out["sub_N"]
    # Mutation fraction (fraction mutated among all unambiguously matching/mutated)
    out["sub_rate"] = out["sub_N"] / out["info"]

    muts_per_vector = vectors.count_muts_by_vec(SUB_N, subsets=True)
    hist_margin = 0.5
    hist_min = 0
    hist_max = muts_per_vector.max()
    if np.isnan(hist_max):
        hist_max = 0
    hist_bins = hist_max - hist_min + 1
    out['sub_hist'] = np.histogram(muts_per_vector,
                                   bins=np.linspace(hist_min - hist_margin,
                                                    hist_max + hist_margin,
                                                    hist_bins + 1))[0]

    out['min_cov'] = min(out['cov'])


    # Add clusters
    if clustering_file:
        ext = os.path.splitext(clustering_file)[1]
        if ext == ".json":
            membership = pd.read_json(clustering_file, orient="columns")
        elif ext == ".csv":
            membership = pd.read_csv(clustering_file)
        else:
            raise ValueError(clustering_file)
        cluster_mus = vectors.get_cluster_mus(membership, SUB_N, subsets=True)
        for cluster in cluster_mus.columns:
            out[f"cluster_{cluster}"] = cluster_mus[cluster]

    out.pop("match")
    out.pop('sequence')
    for k in out:
        if isinstance(out[k], pd.Series):
            out[k] = list(out[k])
    return out
