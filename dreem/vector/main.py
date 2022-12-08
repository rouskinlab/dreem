from dreem import util
import os
import pandas as pd
import sys
sys.path.append(os.path.dirname(os.path.realpath(__file__)))

import mprofile

def run(fasta:str, input_dir, out_dir, library:str=None, parallel:str='auto', coords=[], primers=[], fill:bool=False):
    """Run the vectoring pipeline.

    Turns each bam file into a vector file and outputs them in the directory `out_dir`.

    Parameters from args:
    -----------------------
    fasta: str
        Path to the reference FASTA file.
    input_dir: tuple or str
        Paths to the directory(ies) of bam files.
    out_dir: str
        Path to the output folder.
    library: str
        Name of the library. Default is None.
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
    input_dirs = input_dir if isinstance(input_dir, tuple) else (input_dir,)
    temp_folder = os.path.join(out_dir, 'temp')
    
    # Make the arguments compatible with the vectoring module
    args = {
        'ref_file': fasta,
        'project_dir': out_dir,
        'parallel': parallel,
        'coords': coords,
        'primers': primers,
        'fill': fill
    }

    # Create the folders
    if not os.path.exists(args['project_dir']):
        os.makedirs(args['project_dir'])
    if not os.path.exists(temp_folder): 
        os.makedirs(temp_folder)
    
    # read library
    library = pd.read_csv(library)
    
    for input_dir in input_dirs:
        # Extract the bam files in the input directories
        for bam in os.listdir(input_dir):
            kwargs = args.copy()
            if not bam.endswith('.bam'):
                continue
            if library is not None:
                args['coords'] = [[r['section'], r['section_start'], r['section_end']] for _, r in library[library['construct']==bam.split('.')[0]].iterrows()]
            kwargs['bam_files'] = bam
            mprofile.mp_gen(**kwargs)
    return args

