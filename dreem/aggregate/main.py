import dreem
import os
import click
import pandas as pd
import dreem.util as util
import numpy as np
import string

from dreem.aggregate.samples import add_samples_info
from dreem.aggregate.library import add_library_info
from dreem.aggregate.rnastructure import add_rnastructure_predictions
from dreem.aggregate.poisson import add_poisson_confidence_intervals
from dreem.util.cli_args import INPUT_DIR, LIBRARY, SAMPLES, SAMPLE, CLUSTERING_FILE, OUT_DIR, RNASTRUCTURE_PATH, RNASTRUCTURE_TEMPERATURE, RNASTRUCTURE_FOLD_ARGS, RNASTRUCTURE_DMS, RNASTRUCTURE_DMS_MIN_UNPAIRED_VALUE, RNASTRUCTURE_DMS_MAX_PAIRED_VALUE, POISSON, VERBOSE, COORDS, PRIMERS, FILL, RNASTRUCTURE_PARTITION, RNASTRUCTURE_PROBABILITY

def generate_mut_profile_from_bit_vector(bit_vector, clustering_json, verbose=False):
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
    df : pandas.DataFrame
        The mutation profile (one row per cluster).

    """
    # Read in the bit vector
    df = pd.read_orc(bit_vector)

    # Convert to a mutation profile

    ## TODO: This is a placeholder. Replace with the actual code.


    # Sanity check
    assert len(df) == 1, 'Mutation profile must have only one row.'

    return df




def run(input_dir:str=INPUT_DIR, library:str=LIBRARY, samples:str=SAMPLES, sample:str=SAMPLE, clustering_file:str=CLUSTERING_FILE, out_dir:str=OUT_DIR, rnastructure_path:str=RNASTRUCTURE_PATH, rnastructure_temperature:bool=RNASTRUCTURE_TEMPERATURE, rnastructure_fold_args:str=RNASTRUCTURE_FOLD_ARGS, rnastructure_dms:bool=RNASTRUCTURE_DMS, rnastructure_dms_min_unpaired_value:int=RNASTRUCTURE_DMS_MIN_UNPAIRED_VALUE, rnastructure_dms_max_paired_value:int=RNASTRUCTURE_DMS_MAX_PAIRED_VALUE, rnastructure_partition:bool=RNASTRUCTURE_PARTITION, rnastructure_probability:bool=RNASTRUCTURE_PROBABILITY, poisson:bool=POISSON, verbose:bool=VERBOSE, coords:str=COORDS, primers:str=PRIMERS, fill:bool=FILL):
    """Run the aggregate module.

    Reads in the bit vector files and aggregates them into a single file named [output]/output/aggregate/[name].csv.
    Adds on the library information, the sample information, RNAstructure predictions, and Poisson confidence intervals.

    Parameters from args:
    ---------------------

    input_dir: str
        Path to the bit vector file or list of paths to the bit vector files.
    library: str
        Csv file with the library information.
    samples: str
       Csv file with the sample information.
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
    bit_vector_names = [os.path.basename(f).split('.')[0][:-len('.orc')] for f in input_dir]
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
    
    # Remove this
    raise NotImplementedError('This module is not implemented yet')

    if sample is None:
        if df_samples is not None:
            raise ValueError('If samples is specified, sample must also be specified.')
        if 'bv_dir' in args.keys(): 
            sample = os.path.basename(args['bv_dir']) 
        else:
            sample = 'unnamed_sample_random_id_'+''.join(np.random.choice(string.ascii_lowercase) for _ in range(6)) 

    # Make folders
    output_folder = util.make_folder(os.path.join(root, 'output', 'aggregate') )
    temp_folder = util.make_folder(os.path.join(root, 'temp', 'aggregate') )

    # Read in the bit vectors
    df = {str: pd.DataFrame()}
    
    df_clustering = None if clustering is None else pd.read_json(clustering)
    for construct in bit_vector_names:
        df[construct] = generate_mut_profile_from_bit_vector(bit_vector[bit_vector_names.index(construct)], clustering_json=df_clustering, verbose=verbose)
    df = pd.concat(df, axis=1).reset_index(drop=True)
    
    if df_samples is not None:
        # Add the sample information
        df = add_samples_info(df, df_samples, sample, verbose=verbose)
    
    if df_library is not None:
        # Add the library information
        df = add_library_info(df, df_library, verbose=verbose)

    if rnastructure['path'] is not None:
        # Add RNAstructure predictions
        df = add_rnastructure_predictions(df, rnastructure, sample, verbose=verbose)

    if poisson:
        # Add Poisson confidence intervals
        df = add_poisson_confidence_intervals(df, sample, verbose=verbose)

    # Save the output
    if verbose:
        print('Saving the output to', os.path.join(output_folder, sample+'.csv'))
    df.to_csv(os.path.join(output_folder, sample), index=False)

    return 1
