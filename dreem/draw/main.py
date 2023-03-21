

import os

import pandas as pd
from click import command, pass_obj

from ..util import docdef
from ..util.cli import (DreemCommandName, dreem_command,
                        opt_out_dir, opt_draw_input,
                        opt_flat, opt_coords, opt_primers,
                        opt_mutation_fraction, opt_mutation_fraction_identity,
                        opt_base_coverage, opt_mutations_in_barcodes,
                        opt_mutations_per_read_per_sample,  
                        )
from ..util.dump import *
from .study import Study

params = [
    opt_draw_input,
    opt_out_dir,
    opt_flat,
    opt_coords,
    opt_mutation_fraction,
    opt_mutation_fraction_identity,
    opt_base_coverage,
    opt_mutations_in_barcodes,
    opt_mutations_per_read_per_sample,
]


@command(DreemCommandName.DRAW.value, params=params)
# Pass context object.
@pass_obj
# Turn into DREEM command.
@dreem_command()
def cli(**kwargs):
    return run(**kwargs)


@docdef.auto()
def run( 
        inpt: tuple[str], 
        *,
        flat: list,
        out_dir: str,
        coords: tuple,
        mutation_fraction: bool,
        mutation_fraction_identity:bool,
        base_coverage: bool,
        mutations_in_barcodes: bool,
        mutations_per_read_per_sample: bool,
        ):
    """Run the draw command.

    """


    study = Study(
        data = inpt,
    )
    
    # Check output directory.
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)
 
    for sample in study.get_samples():
        
        path = os.path.join(out_dir, sample)
        if not os.path.isdir(path):
            os.makedirs(path)
        
        refs, starts, ends = select_coords(coords, sample, study)
        
        for ref, start, end in zip(refs, starts, ends):
            
            if flat:
                prefix = '__'.join([ref, str(start)+ '_' + str(end),''])
            else:
                prefix = '/'.join([sample, ref, str(start)+ '_' + str(end) , ''])
            
            section, base_index = find_section_for_these_coords(study, sample, ref, start, end)
            
            if mutation_fraction:
                study.mutation_fraction(
                    sample=sample,
                    reference=ref,
                    section=section,
                    cluster='pop_avg',
                    base_index = base_index,
                    to_html = os.path.join(path,prefix)+'mutation_fraction.html',
                )
            
            if mutation_fraction_identity:
                study.mutation_fraction_identity(
                    sample=sample,
                    reference=ref,
                    section=section,
                    cluster='pop_avg',
                    base_index = base_index,
                    to_html = os.path.join(path,prefix)+'mutation_fraction_identity.html',
                )
            
            if base_coverage:
                study.base_coverage(
                    sample=sample,
                    reference=ref,
                    section=section,
                    cluster='pop_avg',
                    base_index = base_index,
                    to_html = os.path.join(path,prefix)+'base_coverage.html',
                )
                
            if mutations_in_barcodes:
                study.mutations_in_barcodes(
                    sample=sample,
                    reference=ref,
                    section=section,
                    cluster='pop_avg',
                    base_index = base_index,
                    to_html = os.path.join(path,prefix)+'mutations_in_barcodes.html',
                )
                
            if mutations_per_read_per_sample:
                study.mutations_per_read_per_sample(
                    sample=sample,
                    reference=ref,
                    section=section,
                    cluster='pop_avg',
                    base_index = base_index,
                    to_html = os.path.join(path,prefix)+'mutations_per_read_per_sample.html',
                )               


def select_coords(coords, sample, study):
    
    refs, starts, ends = [], [], []
    
    if len(coords):
        for t in coords:
            refs.append(t[0])
            starts.append(t[1])
            ends.append(t[2])
    else:
        for _, row in study.df[study.df['sample']==sample].iterrows():
            refs.append(row['reference'])
            starts.append(row['section_start'])
            ends.append(row['section_end'])
    
    return refs, starts, ends
            
            
def find_section_for_these_coords(study, sample, ref, start, end):
    df = study.df[(study.df['sample']==sample)&(study.df['reference']==ref)]
    for _, row in df.iterrows():
        if row['section_start'] <= start and row['section_end'] >= end:
            return row['section'], np.arange(start-row['section_start'], end-row['section_start']+1)
    raise ValueError(f'No section found for {sample}, {ref}, {start}, {end}. Change your section in your library file or your coordinates.')    


