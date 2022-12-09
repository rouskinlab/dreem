import click
import os
from main import run

@click.command()
@click.option('-o','--out_dir',  default=os.getcwd(), type=click.Path(), help='Where to output files')
@click.option( '-fq1','--fastq1', help='Paths to the fastq1 file (forward primer). Enter multiple times for multiple files', type=click.Path(exists=True), required=True)
@click.option('-fq2', '--fastq2',  help='Paths to the fastq2 file (reverse primer). Enter multiple times for multiple files', type=click.Path(exists=True))
@click.option( '-l','--library', type=click.Path(exists=True), help='Path to the library.csv file')
@click.option('--barcode_start', '-bs', type=int, help='Start position of the barcode in the read')
@click.option('--barcode_end', '-be', type=int, help='End position of the barcode in the read')


def cli(**args):
    run(**args)

if __name__ == '__main__':
    cli()