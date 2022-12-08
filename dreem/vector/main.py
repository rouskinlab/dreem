from dreem import util
import os
import pandas as pd
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import mprofile

def run(**args):
    """Run the vectoring pipeline.

    Turns each bam file into a vector file and outputs them in the directory `out_dir`.

    Parameters from args:
    -----------------------
    fasta: str
        Path to the reference FASTA file.
    input_dir: tuple
        Paths to the directory(ies) of bam files.
    out_dir: str
        Path to the output folder.
    parallel: str
        Parallelize the processing of mutational PROFILES or READS within each profile, turn parallelization OFF, or AUTO matically choose the parallelization method (default: auto).
    coords: tuple
        coordinates for reference: '-c ref-name first last'
    primers: tuple
        primers for reference: '-p ref-name fwd rev'
    fill: bool
        Fill in coordinates of reference sequences for which neither coordinates nor primers were given (default: no).
    
    Returns
    -------
    1 if successful, 0 otherwise.

    """
    # Extract the arguments
    input_dirs = args.pop('input_dir') if isinstance(args['input_dir'], tuple) else [args['input_dir']]
    temp_folder = os.path.join(args['out_dir'],'temp')
    
    # Make the arguments compatible with the vectoring module
    args['project_dir'] = args.pop('out_dir')
    args['ref_file'] = args.pop('fasta')   

    # Create the folders
    if not os.path.exists(args['project_dir']):
        os.makedirs(args['project_dir'])
    if not os.path.exists(temp_folder): 
        os.makedirs(temp_folder)
    
    # read library
    library = pd.read_csv(args['library'])

    for input_dir in input_dirs:
        # Extract the bam files in the input directories
        for bam in [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.bam')]:
            args['coords'] = [[r['section'], r['section_start'], r['section_end']] for _, r in library[library['construct']==bam.split('.')[0]].iterrows()]
            args['bam_files'] = bam
            mprofile.mp_gen(**args)
            print(f"{args=}")

    return 1
