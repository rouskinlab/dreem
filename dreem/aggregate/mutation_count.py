
from ..util.util import *
from ..util.seq import *
import pyarrow.orc as po
import pyarrow as pa

def generate_mut_profile_from_bit_vector(bit_vector, clustering_file, verbose=False):
    """
    Generate a mutation profile from a bit vector.

    Parameters
    ----------
    bit_vector : str
        Path to the bit vector.
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
    # Read in the bit vector 
    bv = pa.concat_tables([po.read_table(os.path.join(bit_vector,b)) for b in os.listdir(bit_vector) if b.endswith(".orc")])
    bv = bv.drop(['__index_level_0__'])
    muts = np.array(bv, dtype=np.uint8).T
    # Convert to a mutation profile
    out = dict()
    out['sequence'] = ''.join([c[0] for c in bv.column_names])
    out['num_aligned'] = muts.shape[0]
    out["match"] = query_muts(muts, MATCH[0] | INS_5[0], set_type='subset')
    out["sub_A"] = query_muts(muts, SUB_A[0], set_type='subset')
    out["sub_C"] = query_muts(muts, SUB_C[0], set_type='subset')
    out["sub_G"] = query_muts(muts, SUB_G[0], set_type='subset')
    out["sub_T"] = query_muts(muts, SUB_T[0], set_type='subset')
    out["sub_N"] = query_muts(muts, SUB_N[0], set_type='subset')
    out["del"]   = query_muts(muts, DELET[0], set_type='subset')
    out["ins"]   = query_muts(muts, INS_3[0], set_type='superset')
    # Can have any mutation, but not a match
    out["cov"] = muts.astype(bool).sum(axis=0)  # i.e. not BLANK
    # Unambiguously matching or mutated (informative)
    out["info"] = out["match"] + out["sub_N"]
    # Mutation rate (fraction mutated among all unambiguously matching/mutated)
    try:
        out["sub_rate"] = out["sub_N"] / out["info"]
    except ZeroDivisionError:
        out["sub_rate"] = [m/i if i != 0 else None for m, i in zip(out["sub_N"], out["info"])]
    
    out['sub_hist'] = np.histogram(query_muts(muts, SUB_N[0], axis=1), bins=range(0, muts.shape[1]))[0] # query_muts(muts, SUB_N[0] | DELET[0] | INS_3[0], axis=1)
    
    out['min_cov'] = min(out['cov'])
    
    out.pop("match")
    for k in out:
        if isinstance(out[k], np.ndarray):
            out[k] = list(out[k])
    return out
