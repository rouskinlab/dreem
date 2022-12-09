import dreem.util
import os
import click
from dreem.alignment import align_reads
from main import run


@click.command()
@click.option('--out_dir', '-o', default=os.getcwd(), type=click.Path(exists=True), help='Where to output files')
@click.option('--sample', type=click.Path(exists=True), help='Name of the sequence alignment map file(s) folder')
@click.option('--fasta', '-fa', type=click.Path(exists=True), help='Path to the fasta file', required=True)
@click.option('--fastq1', '-fq1', help='Paths to the fastq1 file (forward primer). Enter multiple times for multiple files', type=click.Path(exists=True), required=True)
@click.option('--fastq2', '-fq2', help='Paths to the fastq2 file (reverse primer). Enter multiple times for multiple files', type=click.Path(exists=True))

def cli(**args):
    run(**args)

if __name__ == '__main__':
    run()
