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
from dreem.util.cli_args import INPUT_DIR, LIBRARY, SAMPLES, SAMPLE, CLUSTERING_FILE, OUT_DIR, RNASTRUCTURE_PATH, RNASTRUCTURE_TEMPERATURE, RNASTRUCTURE_FOLD_ARGS, RNASTRUCTURE_DMS, RNASTRUCTURE_DMS_MIN_UNPAIRED_VALUE, RNASTRUCTURE_DMS_MAX_PAIRED_VALUE, POISSON, VERBOSE, COORDS, PRIMERS, FILL, RNASTRUCTURE_PARTITION, RNASTRUCTURE_PROBABILITY
sys.path.append(os.path.dirname(__file__))
from mutation_count import generate_mut_profile_from_bit_vector



def run(input_dir:str=INPUT_DIR, library:str=LIBRARY, samples:str=SAMPLES, sample:str=SAMPLE, clustering_file:str=CLUSTERING_FILE, out_dir:str=OUT_DIR, rnastructure_path:str=RNASTRUCTURE_PATH, rnastructure_temperature:bool=RNASTRUCTURE_TEMPERATURE, rnastructure_fold_args:str=RNASTRUCTURE_FOLD_ARGS, rnastructure_dms:bool=RNASTRUCTURE_DMS, rnastructure_dms_min_unpaired_value:int=RNASTRUCTURE_DMS_MIN_UNPAIRED_VALUE, rnastructure_dms_max_paired_value:int=RNASTRUCTURE_DMS_MAX_PAIRED_VALUE, rnastructure_partition:bool=RNASTRUCTURE_PARTITION, rnastructure_probability:bool=RNASTRUCTURE_PROBABILITY, poisson:bool=POISSON, verbose:bool=VERBOSE, coords:str=COORDS, primers:str=PRIMERS, fill:bool=FILL):
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
    library = pd.read_csv(library) if library is not None else None
    df_samples = pd.read_csv(samples) if samples is not None else None
    rnastructure = {}
    rnastructure['path'] = rnastructure_path
    rnastructure['temperature'] = rnastructure_temperature
    rnastructure['fold_args'] = rnastructure_fold_args
    rnastructure['dms'] = rnastructure_dms
    rnastructure['dms_min_unpaired_value'] = rnastructure_dms_min_unpaired_value
    rnastructure['dms_max_paired_value'] = rnastructure_dms_max_paired_value
    rnastructure['partition'] = rnastructure_partition
    rnastructure['probability'] = rnastructure_probability
    rnastructure['temp_folder'] = os.makedirs(os.path.join(out_dir, 'temp', 'rnastructure'), exist_ok=True)
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
    
    mut_profiles = {'constructs': {}}
    for construct in bv_files:
        mut_profiles['constructs'][construct] = {}
        mut_profiles['constructs'][construct]['sections'] = {}
        for file, path in mut_profiles['constructs'][construct]['sections'].items():
            mut_profiles['constructs'][construct]['sections'][file] = generate_mut_profile_from_bit_vector(path, clustering_json=clustering_file, verbose=verbose)
                
    if df_samples is not None:
        # Add the sample information
        mut_profiles = {**mut_profiles, **get_samples_info(df_samples, sample, verbose=verbose)}
    
    for construct in mut_profiles['constructs']:
        if library is not None:
        # Add the library information
            mut_profiles['constructs'][construct] = {**mut_profiles['constructs'][construct], **get_library_info(library, construct, verbose=verbose)}

        if rnastructure['path'] is not None:
        # Add RNAstructure predictions
            mut_profiles['constructs'][construct] = {**mut_profiles['constructs'][construct], **add_rnastructure_predictions(rnastructure, sample, verbose=verbose)}

        if poisson:
        # Add Poisson confidence intervals
            for section in mut_profiles['constructs'][construct]['sections']:
                for cluster in mut_profiles['constructs'][construct]['sections'][section]['clusters']:
                    mut_profiles['constructs'][construct]['sections'][section]['clusters'][cluster] = {**mut_profiles['constructs'][construct]['sections'][section]['clusters'][cluster], **compute_conf_interval(info_bases=cluster['info_bases'], mut_bases=cluster['mut_bases'], verbose=verbose)}

    # Write the output
    with open(os.path.join(out_dir, sample + '.json'), 'w') as f:
        json.dump(mut_profiles, f, indent=4)

    return 1
