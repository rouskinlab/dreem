import os

import pandas as pd
from click import command

from .library_samples import get_samples_info, get_library_info
from .mutation_count import generate_mut_profile_from_bit_vector
from ..util import docdef
from ..util.cli import (opt_out_dir, opt_temp_dir, opt_save_temp,
                        opt_library, opt_samples, opt_fasta,
                        opt_bv_files, opt_clustering_file,
                        opt_rnastructure_path, opt_rnastructure_use_temp,
                        opt_rnastructure_fold_args, opt_rnastructure_use_dms,
                        opt_rnastructure_dms_min_unpaired_value,
                        opt_rnastructure_dms_max_paired_value,
                        opt_rnastructure_deltag_ensemble,
                        opt_rnastructure_probability,
                        opt_verbose)
from ..util.dump import *
from ..util.files_sanity import check_library, check_samples
from ..util.rnastructure import RNAstructure
from ..vector.profile import VectorReader
from ..util.dependencies import *

params = [
    opt_fasta,
    opt_bv_files,
    opt_library,
    opt_samples,
    opt_clustering_file,
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


@command("aggregate", params=params)
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
    
    check_rnastructure_exists(rnastructure_path)

    # Extract the arguments
    if library != '':
        library = check_library(pd.read_csv(library), fasta, out_dir) if library is not None else None
    else:
        library = None
    if samples != '':
        df_samples = check_samples(pd.read_csv(samples)) if samples != '' else None
    else:
        df_samples = None    
    
    # Make folders
    os.makedirs(out_dir, exist_ok=True)
    temp_dir = os.path.join(temp_dir, 'aggregate')
    #os.makedirs(temp_dir, exist_ok=True)

    print('Reading in bit vectors from {}...'.format(bv_files))

    reports_path = [b for b in bv_files if b.endswith('report.json')]
    
    
    for sample_path in [b for b in bv_files if not b.endswith('report.json')]:
        sample = os.path.basename(sample_path)
        for reference in os.listdir(os.path.join(sample_path)):
            for section in os.listdir(os.path.join(sample_path, reference)):
                if os.path.exists(reports_path[-1]):
                    reports_path.append(os.path.join(sample_path, reference, section, 'report.json'))
    
    all_samples = {}
    
    for report_path in reports_path:
        
        vectors = VectorReader.load(report_path)

        if vectors.sample not in all_samples:
            all_samples[vectors.sample] = {}
        if vectors.ref not in all_samples[vectors.sample]:
            all_samples[vectors.sample][vectors.ref] = {}

        all_samples[vectors.sample][vectors.ref][vectors.section] = {
            'section_start': vectors.end5,
            'section_end': vectors.end3,
            'sequence': vectors.seq.decode(),
            'pop_avg': generate_mut_profile_from_bit_vector(vectors, clustering_file=clustering_file, verbose=verbose)
        }
        
    for sample in all_samples:
        for reference in all_samples[sample]:
            # Add the library information 
            lib, section_tranlation = get_library_info(library, reference, verbose=verbose)
            all_samples[sample][reference] = {**lib, **all_samples[sample][reference]}
            
            for section in all_samples[sample][reference].copy().keys():
                if type(all_samples[sample][reference][section]) is not dict:
                    continue
                for col in ['num_aligned']:
                    all_samples[sample][reference][col] = all_samples[sample][reference][section]['pop_avg'].pop(col)
                if section in section_tranlation:
                    all_samples[sample][reference][section_tranlation[section]] = all_samples[sample][reference].pop(section)
            
    # Add the sample information
    print('Adding sample information...')
    for sample in all_samples:
        if df_samples is not None:
            all_samples[sample] = {**all_samples[sample], **get_samples_info(df_samples, sample, verbose=verbose)}
        else:
            all_samples[sample] = {**all_samples[sample], **{'sample': sample}}
    print('Done.')
    
    print('Computing confidence intervals and RNAstructure predictions...')
    
    rna = RNAstructure(rnastructure_path=rnastructure_path, temp = os.path.join(temp_dir,'rnastructure'))
    for sample, mut_profiles in all_samples.items():
        for reference in mut_profiles:
            if type(mut_profiles[reference]) is not dict:
                continue

            for section in mut_profiles[reference]:
                if type(mut_profiles[reference][section]) is not dict:
                    continue
                # Add RNAstructure predictions
                mh = rna.run(mut_profiles[reference][section]['sequence'])
                all_samples[sample][reference][section] = {**mut_profiles[reference][section], **mh}
    rna.dump_ledger()
    
    for sample, mut_profiles in all_samples.items():
        # Write the output
        out = cast_dict(mut_profiles)
        out = sort_dict(out)
        
        # Make lists in one line
        out = json.dumps(out, cls=NpEncoder, indent=2)
        backslashN = """
        """
        out = out.replace(']','[').split('[')
        out = [o.replace(backslashN+'  ','').replace(backslashN, '') if i%2 else o for i, o in enumerate(out)]
        out = out[0] + ''.join([('[',']')[i%2] + o for i, o in enumerate(out[1:])])
        
        # Write the output
        with open(os.path.join(out_dir, sample + '.json'), 'w') as f:
            f.write(out)
        print('Outputed to {}'.format(os.path.join(out_dir, sample + '.json')))
    return 1
