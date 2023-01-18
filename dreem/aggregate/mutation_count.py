
from dreem.util.util import *
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
        - info_bases: the count per residue of bases for which we have information.
        - cov_bases: the count per residue of bases for which we have coverage.
        - mut_bases: the count per residue of bases for which we have mutations.
        - mod_bases_A: the count per residue of bases who are modified to A.
        - mod_bases_C: the count per residue of bases who are modified to C.
        - mod_bases_G: the count per residue of bases who are modified to G.
        - mod_bases_T: the count per residue of bases who are modified to T.
        - mut_rates: the mutation rate per residue.
        - num_aligned: the number of aligned reads.
        
    out['sequence'] = ''.join([c[-1] for c in bv.column_names])
    out['num_aligned'] = muts.shape[0]
    out["match_bases"] = query_muts(muts, MATCH[0])
    out["mod_bases_A"] = query_muts(muts, SUB_A[0])
    out["mod_bases_C"] = query_muts(muts, SUB_C[0])
    out["mod_bases_G"] = query_muts(muts, SUB_G[0])
    out["mod_bases_T"] = query_muts(muts, SUB_T[0])
    out["mod_bases_N"] = query_muts(muts, SUB_N[0])
    out["del_bases"]   = query_muts(muts, DELET[0])
    out["ins_bases"]   = query_muts(muts, INS_5[0] | INS_3[0])
    # Can have any mutation, but not a match
    out["mut_bases"] = query_muts(muts, SUB_N[0] | DELET[0] | INS_5[0] | INS_3[0])
    out["cov_bases"] = muts.astype(bool).sum(axis=0)  # i.e. not BLANK
    # Unambiguously matching or mutated (informative)
    out["info_bases"] = out["match_bases"] + out["mut_bases"]
    # Mutation rate (fraction mutated among all unambiguously matching/mutated)
    try:
        out["mut_rates"] = out["mut_bases"] / out["info_bases"]
    except ZeroDivisionError:
        out["mut_rates"] = [m/i if i != 0 else None for m, i in zip(out["mut_bases"], out["info_bases"])]
    
    out['num_of_mutations'] =  query_muts(muts, SUB_N[0] | DELET[0] | INS_5[0] | INS_3[0], axis=1)

    """
    # Read in the bit vector
    bv = pa.concat_tables([po.read_table(os.path.join(bit_vector,b)) for b in os.listdir(bit_vector) if b.endswith(".orc")])
    muts = np.array(bv, dtype=np.uint8).T
    # Convert to a mutation profile
    out = dict()
    out['sequence'] = ''.join([c[0] for c in bv.column_names])
    out['num_aligned'] = muts.shape[0]
    out["match_bases"] = query_muts(muts, MATCH[0] | INS_5[0], set_type='subset')
    out["mod_bases_A"] = query_muts(muts, SUB_A[0], set_type='subset')
    out["mod_bases_C"] = query_muts(muts, SUB_C[0], set_type='subset')
    out["mod_bases_G"] = query_muts(muts, SUB_G[0], set_type='subset')
    out["mod_bases_T"] = query_muts(muts, SUB_T[0], set_type='subset')
    out["mod_bases_N"] = query_muts(muts, SUB_N[0], set_type='subset')
    out["del_bases"]   = query_muts(muts, DELET[0], set_type='subset')
    out["ins_bases"]   = query_muts(muts, INS_3[0], set_type='superset')
    # Can have any mutation, but not a match
    out["mut_bases"] = out["mod_bases_N"] #query_muts(muts, SUB_N[0] | DELET[0] | INS_3[0], set_type='superset')
    out["cov_bases"] = muts.astype(bool).sum(axis=0)  # i.e. not BLANK
    # Unambiguously matching or mutated (informative)
    out["info_bases"] = out["match_bases"] + out["mut_bases"]
    # Mutation rate (fraction mutated among all unambiguously matching/mutated)
    try:
        out["mut_rates"] = out["mut_bases"] / out["info_bases"]
    except ZeroDivisionError:
        out["mut_rates"] = [m/i if i != 0 else None for m, i in zip(out["mut_bases"], out["info_bases"])]
    
    out['num_of_mutations'] = np.histogram(query_muts(muts, SUB_N[0], axis=1), bins=range(0, muts.shape[1]))[0] # query_muts(muts, SUB_N[0] | DELET[0] | INS_3[0], axis=1)
    
    out['worst_cov_bases'] = min(out['cov_bases'])
    
    out.pop("match_bases")
    for k in out:
        if isinstance(out[k], np.ndarray):
            out[k] = list(out[k])
    return out
