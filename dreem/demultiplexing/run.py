import pandas as pd
import dreem.util as util

from click_option_group import optgroup
import click
import os
from dreem.demultiplexing.demultiplexing import demultiplex

@click.command()
@optgroup.group('Files and folders paths')
@optgroup.option('--output', '-o', default=os.getcwd(), type=click.Path(exists=True), help='Where to output files')
@optgroup.option('--fastq1', '-fq1', help='Paths to the fastq1 file (forward primer). Enter multiple times for multiple files', multiple=True, type=click.Path(exists=True), required=True)
@optgroup.option('--fastq2', '-fq2', help='Paths to the fastq2 file (reverse primer). Enter multiple times for multiple files', multiple=True, type=click.Path(exists=True))
@optgroup.option('--library', '-l', type=click.Path(exists=True), help='Path to the library.csv file')

@optgroup.group('Demultiplexing parameters')
@optgroup.option('--barcode_start', '-bs', type=int, help='Start position of the barcode in the read')
@optgroup.option('--barcode_end', '-be', type=int, help='End position of the barcode in the read')


def run(**args):
    """Run the demultiplexing pipeline.

    Parameters from args:
    -----------------------
    library: str
        Path to the library file. Columns are (non-excusively): ['construct', 'barcode_start', 'barcode_end', 'barcode_sequence']
    fastq1: str
        Path to the FASTQ file or list of paths to the FASTQ files, forward primer.
    fastq2: str
        Path to the FASTQ file or list of paths to the FASTQ files, reverse primer.

    
    Returns
    -------
    1 if successful, 0 otherwise.

    """
    # Get the paths
    root = args['output']
    temp_folder = os.path.join(root,'temp','demultiplexing')
    output_folder = os.path.join(root,'output','demultiplexing')
    fastq1 = args['fastq1']
    fastq2 = args['fastq2']

    # Load the library
    library = pd.read_csv(args['library'])[["construct", "barcode_start", "barcode_end", "barcode_sequence"]].dropna()

    # Make the folders
    print(util.make_folder(output_folder))

    util.make_folder(temp_folder)

    # Demultiplex
    for f1 in args['fastq1']:
        for f2 in args['fastq2']:
            if f1[:-len('_R1.fastq')] == f2[:-len('_R2.fastq')]:
                sample = f1.split('/')[-1][:-len('_R1.fastq')]
                util.make_folder(os.path.join(output_folder, sample))

                ### TODO: Implement demultiplexing
                demultiplex(f1, f2, library, os.path.join(output_folder, sample), os.path.join(temp_folder, sample))

    return 1


if __name__ == '__main__':
    run()