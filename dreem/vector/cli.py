import dreem
import click
import dreem.util as util
import os, sys
from main import run

@click.command()
@click.option('--out_dir', '-o', default=os.getcwd(), type=click.Path(exists=True), help='Where to output files')
@click.option('--fasta', '-fa', type=click.Path(exists=True), help='Path to the fasta file', required=True)
@click.option('--input_dir', '-id', help='Paths to the directory(ies) of bam files', type=click.Path(exists=True), required=True, multiple=True)
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


def cli(**args):
    run(**args)
    
if __name__ == '__main__':
    run()