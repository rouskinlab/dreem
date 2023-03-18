
from ..util.util import *
from ..util.seq import *
import pyarrow.orc as po
import pyarrow as pa

from ..vector.profile import VectorReader


def generate_mut_profile_from_bit_vector(vectors: VectorReader, clustering_file, verbose=False):
    """
    Generate a mutation profile from mutation vectors.

    Parameters
    ----------
    vectors: VectorReader
        Vector reader object
    verbose : bool
        If True, print progress.

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
        - sub_rate: the mutation rate per residue.
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
    out["del"]   = vectors.count_muts_by_pos(DELET)
    out["ins"]   = vectors.count_muts_by_pos(INS_3, supersets=True)
    # Count non-blank bytes
    out["cov"] = vectors.count_muts_by_pos(EVERY, subsets=True)
    # Unambiguously matching or substituted (informative)
    out["info"] = out["match"] + out["sub_N"]
    # Mutation rate (fraction mutated among all unambiguously matching/mutated)
    try:
        out["sub_rate"] = out["sub_N"] / out["info"]
    except ZeroDivisionError:
        out["sub_rate"] = [m/i if i != 0 else None for m, i in zip(out["sub_N"], out["info"])]

    muts_per_vector = vectors.count_muts_by_vec(SUB_N, subsets=True)
    hist_margin = 0.5
    hist_min = 0
    hist_max = vectors.length
    # hist_max = muts_per_vector.max()  # FIXME: I think muts_per_vector.max() would be a better hist_max because vectors.length will produce a lot of zeros on the right
    hist_bins = hist_max - hist_min + 1
    out['sub_hist'] = np.histogram(muts_per_vector,
                                   bins=np.linspace(hist_min - hist_margin,
                                                    hist_max + hist_margin,
                                                    hist_bins + 1))[0]
    
    out['min_cov'] = min(out['cov'])
    
    out.pop("match")
    out.pop('sequence')
    for k in out:
        if isinstance(out[k], pd.Series):
            out[k] = list(out[k])
    return out
