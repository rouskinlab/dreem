import dreem
import click
import dreem.util as util
import os, sys
from dreem.vectoring.vectoring import generate_bitvectors

@click.command()
@click.option('--output', '-o', default=os.getcwd(), type=click.Path(exists=True), help='Where to output files')
@click.option('--fasta', '-fa', type=click.Path(exists=True), help='Path to the fasta file', required=True)
@click.option('--bam_dir', '-bd', help='Paths to the directory(ies) of bam files', type=click.Path(exists=True), required=True, multiple=True)
@click.option("-P", "--parallel",
              type=click.Choice(["profiles", "reads", "off", "auto"],
                                case_sensitive=False),
              default="auto",
              help="Parallelize the processing of mutational PROFILES or "
              "READS within each profile, turn parallelization OFF, or AUTO"
              "matically choose the parallelization method (default: auto).")
@click.option("-c", "--coords", type=(str, int, int), multiple=True, help="coordinates for reference: '-c ref-name first last'")
@click.option("-p", "--primers", type=(str, util.Primer, util.Primer), multiple=True, help="primers for reference: '-p ref-name fwd rev'")
@click.option("--fill/--no-fill", default=False,
              help="Fill in coordinates of reference sequences for which "
                   "neither coordinates nor primers were given (default: no).")

def run(**args):
    """Run the vectoring pipeline.

    Turns each bam file into a vector file. 

    Parameters from args:
    -----------------------
    fasta: str
        Path to the reference FASTA file.
    bam_dir: tuple
        Paths to the directory(ies) of bam files.
    output: str
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
    bam_dirs = args['bam_dir']
    root = args['output']
    parallel = args['parallel']
    coords = args['coords']
    primers = args['primers']
    fill = args['fill']
    temp_folder = os.path.join(root,'temp','vectoring')
    output_folder = os.path.join(root,'output','vectoring')

    # Create the folders
    util.make_folder(output_folder)
    util.make_folder(temp_folder)

    assert generate_bitvectors(fasta, bam_dirs, output_folder, temp_folder, parallel, coords, primers, fill), "Vectoring failed"

    return 1