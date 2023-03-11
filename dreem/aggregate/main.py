import os, sys
import pandas as pd
from click import command, pass_obj

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
from ..util.cli import (DreemCommandName, dreem_command,
                        opt_out_dir, opt_temp_dir, opt_save_temp,
                        opt_library, opt_samples, opt_fasta,
                        opt_bv_files,
                        opt_sample, opt_clustering_file,
                        opt_rnastructure_path, opt_rnastructure_use_temp,
                        opt_rnastructure_fold_args, opt_rnastructure_use_dms,
                        opt_rnastructure_dms_min_unpaired_value,
                        opt_rnastructure_dms_max_paired_value,
                        opt_rnastructure_deltag_ensemble,
                        opt_rnastructure_probability,
                        opt_verbose)
params = [
    opt_fasta,
    opt_bv_files,
    opt_library,
    opt_samples,
    opt_clustering_file,
    opt_sample,
    opt_rnastructure_path,
    opt_rnastructure_use_temp,
    opt_rnastructure_fold_args,
    opt_rnastructure_use_dms,
    opt_rnastructure_dms_min_unpaired_value,
    opt_rnastructure_dms_max_paired_value,
    opt_rnastructure_deltag_ensemble,
    opt_rnastructure_probability,
    opt_verbose,
    opt_out_dir,
    opt_temp_dir,
    opt_save_temp,
]


@command(DreemCommandName.AGGREGATE.value, params=params)
# Pass context object.
@pass_obj
# Turn into DREEM command.
@dreem_command(imports=("fasta", "bv_files"),
               result_key="dreem_output")
def cli(**kwargs):
    return run(**kwargs)


@docdef.auto()
def run( 
        fasta: str, 
        bv_files: list,
        *,
        out_dir: str, 
        temp_dir: str,
        save_temp: bool,
        library: str, 
        samples: str, 
        sample: str, 
        clustering_file: str, 
        rnastructure_path: str, 
        rnastructure_use_temp: bool,
        rnastructure_fold_args: str,
        rnastructure_use_dms: str,
        rnastructure_dms_min_unpaired_value: float,
        rnastructure_dms_max_paired_value: float,
        rnastructure_deltag_ensemble: bool,
        rnastructure_probability: bool,
        verbose: bool):
    """Run the aggregate module.

    """

    # Extract the arguments
    if library != '':
        library = check_library(pd.read_csv(library), fasta, out_dir) if library is not None else None
        library['section_boundaries'] = library.apply(lambda x: str(x['section_start']) + '-' + str(x['section_end']), axis=1)
    else:
        library = None
    if samples != '':
        df_samples = check_samples(pd.read_csv(samples)) if samples != '' else None
    else:
        df_samples = None
    fasta = pd.DataFrame({k :v.decode("utf-8")  for k,v in parse_fasta(fasta)}, index=[0]).T.reset_index().rename(columns={"index":"reference", 0:"sequence"})
    
    
    # Make folders
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(os.path.join(out_dir, 'temp'), exist_ok=True)

    # Read in the bit mut_vectors
    if clustering_file != '':
        with open(clustering_file, 'r') as f:
            clustering_file = json.load(f)    
    
    print('Reading in bit mut_vectors from {}...'.format(bv_files))
    
    reports_path = [b for b in bv_files if b.endswith('report.json')]
    bv_samples = [b for b in bv_files if not b.endswith('report.json')]
    
    all_samples = {os.path.dirname(os.path.dirname(s)):{} for s in reports_path}
    
    for sample in bv_samples:
                
        mut_profiles = {}

        for reference in os.listdir(os.path.join(sample)):
            
            mut_profiles[reference] = {'sequence': fasta[fasta['reference'] == reference]['sequence'].values[0]}
            
            # Add the library information 
            mut_profiles[reference] = {**get_library_info(library, reference, verbose=verbose), **mut_profiles[reference]}
        
            if not len(os.listdir(os.path.join(sample, reference))) > 0:
                logging.error('No bit vectors found for reference {}'.format(reference))
                continue
            
            for section in os.listdir(os.path.join(sample, reference)):
                bv = os.path.join(sample, reference, section)
                report = json.load(open(os.path.join(sample, reference, section, 'report.json'),'r'))
                
                assert report['reference'] == reference, 'Reference in report does not match reference in file path'
                assert report['sample'] == sample, 'Sample in report does not match sample in file path'
                
                mut_profiles[reference][section] = {
                    'section_start': report['section_start'],
                    'section_end': report['section_end'],
                    'sequence': report['sequence'],
                    'pop_avg': generate_mut_profile_from_bit_vector(bv, clustering_file=clustering_file, verbose=verbose)
                }
        
        for col in ['num_aligned']:
            mut_profiles[reference][col] = mut_profiles[reference][section]['pop_avg'].pop(col)
            
        all_samples[sample] = mut_profiles
    
    for report_path in reports_path:
        
        report = json.load(open(os.path.join(report_path),'r'))
        section_path = os.path.dirname(report_path)
        sample, reference, section = report['sample'], report['reference'], os.path.basename(section_path)
        
        all_samples[sample][reference][section] = {
            'section_start': report['section_start'],
            'section_end': report['section_end'],
            'sequence': report['sequence'],
            'pop_avg': generate_mut_profile_from_bit_vector(section_path, clustering_file=clustering_file, verbose=verbose)
        }
    

    print('Done.')
    if df_samples is not None:
        # Add the sample information
        print('Adding sample information...')
        for sample in all_samples:
            all_samples[sample] = {**all_samples[sample], **get_samples_info(df_samples, sample, verbose=verbose)}
        print('Done.')
        
    
    print('Computing confidence intervals and RNAstructure predictions...')
    print(all_samples)
    # rna = RNAstructure(rnastructure_path=rnastructure_path)
    # for sample, mut_profiles in all_samples.items():
    #     for reference in mut_profiles:
    #         if type(mut_profiles[reference]) is not dict:
    #             continue

    #         for section in mut_profiles[reference]:
    #             if type(mut_profiles[reference][section]) is not dict:
    #                 continue
    #             # Add RNAstructure predictions

    #             mh = rna.run(mut_profiles[reference][section]['sequence'])
    #             all_samples[sample][reference][section] = {**mut_profiles[reference][section], **mh}


    for sample, mut_profiles in all_samples.keys():
        # Write the output
        out = cast_dict(mut_profiles)
        out = sort_dict(out)
        with open(os.path.join(out_dir, sample + '.json'), 'w') as f:
            f.write(json.dumps(out, cls=NpEncoder, indent=2))
        print('Outputted to {}'.format(os.path.join(out_dir, sample + '.json')))
    return 1
