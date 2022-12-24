import dreem
import os, sys
import click
import pandas as pd
import dreem.util as util
import numpy as np
import string
import json

from dreem.aggregate.library_samples import get_samples_info, get_library_info
from dreem.aggregate.rnastructure import add_rnastructure_predictions
from dreem.aggregate.poisson import compute_conf_interval
from dreem.util.cli import INPUT_DIR, LIBRARY, SAMPLES, SAMPLE, CLUSTERING_FILE, OUT_DIR, FASTA, RNASTRUCTURE_PATH, RNASTRUCTURE_TEMPERATURE, RNASTRUCTURE_FOLD_ARGS, RNASTRUCTURE_DMS, RNASTRUCTURE_DMS_MIN_UNPAIRED_VALUE, RNASTRUCTURE_DMS_MAX_PAIRED_VALUE, POISSON, VERBOSE, COORDS, PRIMERS, FILL, RNASTRUCTURE_PARTITION, RNASTRUCTURE_PROBABILITY
sys.path.append(os.path.dirname(__file__))
from mutation_count import generate_mut_profile_from_bit_vector
from dreem.util.files_sanity import check_library, check_samples
from rnastructure import RNAstructure
from dreem.util.fa import parse_fasta


def run(input_dir:str=INPUT_DIR, library:str=LIBRARY, samples:str=SAMPLES, sample:str=SAMPLE, clustering_file:str=CLUSTERING_FILE, out_dir:str=OUT_DIR, fasta:str = FASTA, rnastructure_path:str=RNASTRUCTURE_PATH, rnastructure_temperature:bool=RNASTRUCTURE_TEMPERATURE, rnastructure_fold_args:str=RNASTRUCTURE_FOLD_ARGS, rnastructure_dms:bool=RNASTRUCTURE_DMS, rnastructure_dms_min_unpaired_value:int=RNASTRUCTURE_DMS_MIN_UNPAIRED_VALUE, rnastructure_dms_max_paired_value:int=RNASTRUCTURE_DMS_MAX_PAIRED_VALUE, rnastructure_partition:bool=RNASTRUCTURE_PARTITION, rnastructure_probability:bool=RNASTRUCTURE_PROBABILITY, poisson:bool=POISSON, verbose:bool=VERBOSE, coords:str=COORDS, primers:str=PRIMERS, fill:bool=FILL):
    """Run the aggregate module.

    Reads in the bit vector files and aggregates them into a single file named [output]/output/aggregate/[name].csv.
    Adds on the library information, the sample information, RNAstructure predictions, and Poisson confidence intervals.

    Parameters from args:
    ---------------------

    input_dir: str
        Path to the bit vector file or list of paths to the bit vector files.
    library: str
        Path to a csv file with the library information.
    samples: str
        Path to a csv file with the sample information.
    sample: str
        Name to identify the row in samples.csv. Also the name for the output file. Default is the containing folder name.
    clustering_file: str
        Path to the clustering.json file.
    out_dir: str
        Path to the output folder (the sample).
    fasta: str
        Path to the fasta file.
    rnastructure_path: str
        Path to RNAstructure, to predict structure and free energy.
    rnastructure_temperature: bool
        Use sample.csv temperature values for RNAstructure.
    rnastructure_fold_args: str
        Arguments to pass to RNAstructure fold.
    rnastructure_dms: bool
        Use the DMS signal to amke predictions with RNAstructure.
    rnastructure_dms_min_unpaired_value: int
        Minimum unpaired value for using the dms signal as an input for RNAstructure.
    rnastructure_dms_max_paired_value: int
        Maximum paired value for using the dms signal as an input for RNAstructure.
    rnastructure_partition: bool
        Use RNAstructure partition function to predict free energy.
    rnastructure_probability: bool
        Use RNAstructure probability to predict free energy.
    poisson: bool
        Predict Poisson confidence intervals.
    verbose: bool
        Verbose output.
    coords: tuple
        coordinates for reference: '-c ref-name first last'
    primers: tuple
        primers for reference: '-p ref-name fwd rev'
    fill: bool
        Fill in coordinates of reference sequences for which neither coordinates nor primers were given (default: no).
        
    Returns:
    --------
    1 if successful, 0 otherwise.

    """

    # Extract the arguments
    library = check_library(pd.read_csv(library), fasta) if library is not None else None
    df_samples = check_samples(pd.read_csv(samples)) if samples is not None else None
    fasta = pd.DataFrame({k.decode("utf-8") :v.decode("utf-8")  for k,v in parse_fasta(fasta)}, index=[0]).T.reset_index().rename(columns={"index":"construct", 0:"sequence"})
    
    rnastructure = {
        'path': rnastructure_path,
        'temperature': rnastructure_temperature,
        'fold_args': rnastructure_fold_args,
        'dms': rnastructure_dms,
        'dms_min_unpaired_value': rnastructure_dms_min_unpaired_value,
        'dms_max_paired_value': rnastructure_dms_max_paired_value,
        'partition': rnastructure_partition,
        'probability': rnastructure_probability,
        'temp_folder': os.path.join(out_dir, 'temp', 'rnastructure')
    }
    
    bv_files = {}
    for construct in os.listdir(input_dir):
        if os.path.isfile(os.path.join(input_dir, construct)):
            continue
        bv_files[construct] = {}
        for file in os.listdir(os.path.join(input_dir, construct)):
            if file.endswith('.orc'):
                bv_files[construct][file.split('.')[0]] = os.path.join(input_dir, construct, file)
    
    # Find a name for the sample
    if sample is None:
        if df_samples is not None:
            raise ValueError('If samples is specified, sample must also be specified.')
        sample = os.path.basename(input_dir) 
    
    # Make folders
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.join(out_dir, 'temp'), exist_ok=True)

    # Read in the bit vectors
    if clustering_file is not None:
        with open(clustering_file, 'r') as f:
            clustering_file = json.load(f)    
    
    mut_profiles = {}
    for construct in bv_files:
        assert len(bv_files[construct]) > 0, 'No bit vectors found for construct {}'.format(construct)
        mut_profiles[construct] = {'sequence': fasta[fasta['construct'] == construct]['sequence'].values[0]}
        for section, path in bv_files[construct].items():
            assert library[(library['construct'] == construct)&(library['section'] == section)].shape[0] == 1, 'Library information not found for construct {} section {}'.format(construct, section)
            mut_profiles[construct][section] = {}
            mut_profiles[construct][section]['pop_avg'] = generate_mut_profile_from_bit_vector(path, clustering_file=clustering_file, verbose=verbose)
            for col in ['sequence']:
                mut_profiles[construct][section][col] = mut_profiles[construct][section]['pop_avg'].pop(col)
            for col in ['num_aligned']:
                mut_profiles[construct][col] = mut_profiles[construct][section]['pop_avg'].pop(col)
            mut_profiles[construct][section]['section_start'] = library[(library['construct'] == construct)&(library['section'] == section)]['section_start'].values[0]
            mut_profiles[construct][section]['section_end'] = library[(library['construct'] == construct)&(library['section'] == section)]['section_end'].values[0]
    
    if df_samples is not None:
        # Add the sample information
        mut_profiles = {**mut_profiles, **get_samples_info(df_samples, sample, verbose=verbose)}
    
    if rnastructure['path'] is not None:
        rna = RNAstructure(rnastructure)
    
    for construct in mut_profiles:
        if type(mut_profiles[construct]) is not dict:
            continue
        # Add the library information
        if library is not None:
            mut_profiles[construct] = {**mut_profiles[construct], **get_library_info(library, construct, verbose=verbose)}

        for section in mut_profiles[construct]:
            if type(mut_profiles[construct][section]) is not dict:
                continue
            # Add RNAstructure predictions

            if rnastructure['path'] is not None:
                mh = rna.run(mut_profiles[construct][section], sample, sequence_only=True)
                mut_profiles[construct][section] = {**mut_profiles[construct][section], **mh}

            continue
            for cluster in mut_profiles['constructs'][construct]['sections'][section]['clusters']:
                # Add Poisson confidence intervals
                if poisson:
                    mut_profiles['constructs'][construct]['sections'][section]['clusters'][cluster] = {**mut_profiles['constructs'][construct]['sections'][section]['clusters'][cluster], **compute_conf_interval(info_bases=cluster['info_bases'], mut_bases=cluster['mut_bases'], verbose=verbose)}
    # Write the output
    with open(os.path.join(out_dir, sample + '.json'), 'w') as f:
        json.dump(sort_dict(cast_dict(mut_profiles)), f, indent=4, cls=NpEncoder)

    return 1


def cast_to_json_compat(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj

def cast_dict(mp):
    for k, v in mp.items():
        mp[k] = cast_to_json_compat(v)
        if type(v) is dict:
            cast_dict(v)
    return mp
    
def sort_dict(mut_profiles):
    sorting_key = lambda item: 1000*{int:0, str:0, float:0, list:1, dict:2}[type(item[1])] + ord(item[0][0])
    mut_profiles = {k:mut_profiles[k] for k in sorted(mut_profiles)}
    mut_profiles = dict(sorted(mut_profiles.items(), key=sorting_key))
    for k, v in mut_profiles.items():
        if type(v) is not dict:
            continue
        mut_profiles[k] = dict(sorted(v.items(), key=sorting_key))
        for kk, vv, in v.items():
            if type(vv) is not dict:
                continue
            mut_profiles[k][kk] = dict(sorted(vv.items(), key=sorting_key))
    return mut_profiles


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)