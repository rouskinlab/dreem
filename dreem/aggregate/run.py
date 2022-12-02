import dreem
import os
import click
import pandas as pd
import dreem.util as util
import numpy as np
import string
from dreem.aggregate.aggregate import generate_mut_profile_from_bit_vector
from dreem.aggregate.samples import add_samples_info
from dreem.aggregate.library import add_library_info
from dreem.aggregate.rnastructure import add_rnastructure_predictions
from dreem.aggregate.poisson import add_poisson_confidence_intervals

@click.command()
@click.option('--bit_vector', '-bv', help='Path to the bit vector files', type=click.Path(exists=True), multiple=True)
@click.option('--bv_dir', '-bvd', help='Path to the folder containing the bit vector files', type=click.Path(exists=True))
@click.option('--samples', '-s', type=click.Path(exists=True), help='Path to the samples.csv file', default=None)
@click.option('--sample', type=str, help='Name to identify the row in samples.csv. Also the name for the output file.', default=None)
@click.option('--clustering', '-cl', type=click.Path(exists=True), help='Path to the clustering.json file', default=None)
@click.option('--library', '-l', type=click.Path(exists=True), help='Path to the library.csv file', default=None)
@click.option('--out_dir', '-o', default=os.getcwd(), type=click.Path(exists=True), help='Where to output the files')
@click.option('--rnastructure_path', '-rs', type=click.Path(exists=True), help='Path to RNAstructure, to predict structure and free energy', default=None)
@click.option('--rnastructure_temperature', '-rst', type=bool, help='Use sample.csv temperature values for RNAstructure', default=False)
@click.option('--rnastructure_fold_args', '-rsa', type=str, help='Arguments to pass to RNAstructure fold', default=None)
@click.option('--rnastructure_dms', '-rsd', type=bool, help='Use the DMS signal to amke predictions with RNAstructure', default=False)
@click.option('--rnastructure_dms_min_unpaired_value', '-rsdmin', type=int, help='Minimum unpaired value for using the dms signal as an input for RNAstructure', default=0.01)
@click.option('--rnastructure_dms_max_paired_value', '-rsdmax', type=int, help='Maximum paired value for using the dms signal as an input for RNAstructure', default=0.05)
@click.option('--rnastructure_partition', '-rspa', type=bool, help='Use RNAstructure partition function to predict free energy', default=False)
@click.option('--rnastructure_probability', '-rspr', type=bool, help='Use RNAstructure partition function to predict per-base mutation probability', default=False)
@click.option('--poisson', '-po', type=bool, help='Predict Poisson confidence intervals', default=True)
@click.option('--verbose', '-v', type=bool, help='Verbose output', default=False)


def run(**args):
    """Run the aggregate module.

    Reads in the bit vector files and aggregates them into a single file named [output]/output/aggregate/[name].csv.
    Adds on the library information, the sample information, RNAstructure predictions, and Poisson confidence intervals.

    Parameters from args:
    ---------------------

    bit_vector: str
        Path to the bit vector file or list of paths to the bit vector files.
    library: str
        Csv file with the library information.
    samples: str
       Csv file with the sample information.
    sample: str
        Name to identify the row in samples.csv. Also the name for the output file.
    clustering: str
        Path to the clustering.json file.
    out_dir: str
        Path to the output folder (the sample).
    name: str
        Name for the output file, for example the sample.
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
    poisson: bool
        Predict Poisson confidence intervals.
    verbose: bool
        Verbose output.
    

    Returns:
    --------
    1 if successful, 0 otherwise.

    """

    # Extract the arguments
    if 'bit_vector' in args.keys():
        bit_vector = args['bit_vector']
    elif 'bv_dir' in args.keys():
        bit_vector = [os.path.join(args['bv_dir'], f) for f in os.listdir(args['bv_dir']) if f.endswith('.orc')]
    else:
        raise ValueError('Either bit_vector or bv_dir must be specified.')
    bit_vector_names = [os.path.basename(f).split('.')[0][:-len('.orc')] for f in bit_vector]
    df_library = None if args['library'] is None else pd.read_csv(args['library'])
    df_samples = None if args['samples'] is None else pd.read_csv(args['samples'])
    root = args['out_dir']
    clustering = args['clustering']
    rnastructure = {}
    rnastructure['path'] = args['rnastructure_path']
    rnastructure['temperature'] = args['rnastructure_temperature']
    rnastructure['fold_args'] = args['rnastructure_fold_args']
    rnastructure['dms'] = args['rnastructure_dms']
    rnastructure['dms_min_unpaired_value'] = args['rnastructure_dms_min_unpaired_value']
    rnastructure['dms_max_paired_value'] = args['rnastructure_dms_max_paired_value']
    rnastructure['partition'] = args['rnastructure_partition']
    rnastructure['probability'] = args['rnastructure_probability']
    rnastructure['temp_folder'] = util.make_folder(os.path.join(root, 'temp','aggregate','rnastructure'))
    poisson = args['poisson']
    verbose = args['verbose']

    sample = args['sample']
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

if __name__ == '__main__':
    run()