import os, sys
import pandas as pd

from .library_samples import get_samples_info, get_library_info
#sys.path.append(os.path.dirname(__file__))
from .mutation_count import generate_mut_profile_from_bit_vector
from ..util.files_sanity import check_library, check_samples
from ..util.rnastructure import RNAstructure
from ..util.seq import parse_fasta
from ..util.dump import *
import logging
from ..util import docdef
from ..util import path

@docdef.auto()
def run( 
        fasta: str, 
        bv_files: list,
        *,
        library: str, 
        samples: str, 
        sample: str, 
        clustering_file: str, 
        out_dir: str, 
        rnastructure_path: str, 
        verbose: bool):
    """Run the aggregate module.

    Reads in the bit vector files and aggregates them into a single file named [output]/output/aggregate/[name].csv.
    Adds on the library information, the sample information, and RNAstructure predictions.

    Parameters from args:
    ---------------------

    bv_files: str
        Path to the bit vector files folders.
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
    verbose: bool
        Verbose output.

    Returns:
    --------
    1 if successful, 0 otherwise.

    """

    # Extract the arguments
    library = check_library(pd.read_csv(library), fasta, out_dir) if library is not None else None
    library['section_boundaries'] = library.apply(lambda x: str(x['section_start']) + '-' + str(x['section_end']), axis=1)
    df_samples = check_samples(pd.read_csv(samples)) if samples != '' else None
    fasta = pd.DataFrame({k :v.decode("utf-8")  for k,v in parse_fasta(fasta)}, index=[0]).T.reset_index().rename(columns={"index":"reference", 0:"sequence"})
    
    # Find a name for the sample
    if sample is None:
        if df_samples is not None:
            raise ValueError('If samples is specified, sample must also be specified.')
        sample = 'output'
    
    # Make folders
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.join(out_dir, 'temp'), exist_ok=True)

    # Read in the bit mut_vectors
    if clustering_file is not None:
        with open(clustering_file, 'r') as f:
            clustering_file = json.load(f)    
    
    mut_profiles = {}
    print('Reading in bit mut_vectors from {}...'.format(bv_files))
    
    bv_files = [b.replace('_report.json', '') for b in bv_files if b.endswith('_report.json')]
    
    for bv in bv_files:
        
        reference, section = bv.split('/')[-2], bv.split('/')[-1]
        
        assert len(library[(library['reference'] == reference)&(library['section_boundaries'] == section)]) < 2, 'Library information not unique for reference {} section {}'.format(reference, section)
        assert len(library[(library['reference'] == reference)&(library['section_boundaries'] == section)]) > 0, 'Library information not existing for reference {} section {}'.format(reference, section)
        section = library[(library['reference'] == reference)&(library['section_boundaries'] == section)]['section'].values[0]
        
        if not len(os.listdir(bv)) > 0:
            logging.warning('No bit vectors found for reference {}'.format(reference))
            continue
        
        if reference not in mut_profiles:
            mut_profiles[reference] = {'sequence': fasta[fasta['reference'] == reference]['sequence'].values[0]}
        # Add the library information
        mut_profiles[reference] = {**get_library_info(library, reference, verbose=verbose), **mut_profiles[reference]}
        mut_profiles[reference].pop('section_boundaries')
        
        assert library[(library['reference'] == reference)&(library['section'] == section)].shape[0] == 1, 'Library information not found for reference {} section {}'.format(reference, section)
        mut_profiles[reference][section] = {}
        mut_profiles[reference][section]['section_start'] = library[(library['reference'] == reference)&(library['section'] == section)]['section_start'].values[0]
        mut_profiles[reference][section]['section_end'] = library[(library['reference'] == reference)&(library['section'] == section)]['section_end'].values[0]
        mut_profiles[reference][section]['pop_avg'] = generate_mut_profile_from_bit_vector(bv, clustering_file=clustering_file, verbose=verbose)
        
        # TODO
        mut_profiles[reference][section]['pop_avg'].pop('sequence')
        mut_profiles[reference][section]['sequence'] = mut_profiles[reference]['sequence'][mut_profiles[reference][section]['section_start']-1:mut_profiles[reference][section]['section_end']]# mut_profiles[reference][section]['pop_avg'].pop('sequence')
        # assert mut_profiles[reference]['sequence'][mut_profiles[reference][section]['section_start']-1:mut_profiles[reference][section]['section_end']] == mut_profiles[reference][section]['sequence'], 'Sequence mismatch for reference {} section {}: {} vs {}'.format(reference, section, mut_profiles[reference]['sequence'][mut_profiles[reference][section]['section_start']-1:mut_profiles[reference][section]['section_end']], mut_profiles[reference][section]['sequence'])
        
        
        for col in ['num_aligned']:
            mut_profiles[reference][col] = mut_profiles[reference][section]['pop_avg'].pop(col)

    print('Done.')
    if df_samples is not None:
        # Add the sample information
        print('Adding sample information...')
        mut_profiles = {**mut_profiles, **get_samples_info(df_samples, sample, verbose=verbose)}
        print('Done.')
        
    
    print('Computing confidence intervals and RNAstructure predictions...')
    rna = RNAstructure(rnastructure_path=rnastructure_path)
    for reference in mut_profiles:
        if type(mut_profiles[reference]) is not dict:
            continue

        for section in mut_profiles[reference]:
            if type(mut_profiles[reference][section]) is not dict:
                continue
            # Add RNAstructure predictions

            mh = rna.run(mut_profiles[reference][section]['sequence'])
            mut_profiles[reference][section] = {**mut_profiles[reference][section], **mh}

    print('Done.')
    # Write the output
    print('Cast dictionary, size:', sys.getsizeof(mut_profiles))
    out = cast_dict(mut_profiles)
    print('Done.')
    print('Sort dictionary, size:', sys.getsizeof(mut_profiles))
    out = sort_dict(out)
    print('Done.')
    print('Dump the json, size', sys.getsizeof(json.dumps(out, cls=NpEncoder)))
    with open(os.path.join(out_dir, sample + '.json'), 'w') as f:
        f.write(json.dumps(out, cls=NpEncoder, indent=2))
    print('Done.')
    print('Done aggregating the data.')
    return 1
