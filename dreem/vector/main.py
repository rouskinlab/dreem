from dreem import util
import os


def generate_bitvectors(fasta, bam_dirs, output_folder, temp_folder, parallel, coords, primers, fill):
    """Generate the bitvectors for the given bam files.

    Parameters
    ----------
    fasta: str
        Path to the reference fasta file.
    bam_dirs: list of str
        List of paths to the bam files.
    output_folder: str
        Path to the output folder.
    temp_folder: str
        Path to the temp folder.
    parallel: int
        Number of parallel processes.
    coords: str
        Path to the coordinates file.
    primers: str
        Path to the primers file.
    fill: str
        Path to the fill file.

    Returns
    -------
    1 if successful, 0 otherwise.

    """
    try:
        __generate_bitvectors(fasta, bam_dirs, output_folder, temp_folder, parallel, coords, primers, fill)
    except Exception as e:
        print(e)
        return 0

def __generate_bitvectors(fasta, bam_dirs, output_folder, temp_folder, parallel, coords, primers, fill):

    for bam_dir in bam_dirs:
        path = util.make_folder(os.path.join(output_folder, bam_dir.split('/')[-1]))
        bam_files = [os.path.join(bam_dir, f) for f in os.listdir(bam_dir) if f.endswith('.bam')]
        bitvectors = [BAM2bitvector(bam, path) for bam in bam_files]


def BAM2bitvector(bam, output):
    """Convert a bam file to a bitvector file.

    Parameters
    ----------
    bam: str
        Path to the bam file.
    output: str
        Path to the output folder.

    Returns
    -------
    bitvector: str
        Path to the bitvector file.

    """
    
    # PLACEHOLDER TODO

    # Create the bitvector file
    bitvector = os.path.join(output, bam.split('/')[-1].replace('.bam', '.ord'))
    util.make_folder(output)

    # Create the bitvector
    with open(bitvector, 'w') as f:
        f.write('test')

    return bitvector



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
    fasta = args['fasta']
    input_dirs = args['input_dir']
    root = args['out_dir']
    temp_folder = os.path.join(root,'temp','vectoring')
    output_folder = os.path.join(root,'output','vectoring')

    # Create the folders
    util.make_folder(output_folder)
    util.make_folder(temp_folder)

    # Remove this
    raise NotImplementedError('This module is not implemented yet')

    assert generate_bitvectors(fasta, input_dirs, output_folder, temp_folder, parallel, coords, primers, fill), "Vectoring failed"

    return 1