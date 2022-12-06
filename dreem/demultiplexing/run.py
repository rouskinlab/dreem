import pandas as pd
import dreem.util as util

#from click_option_group import optgroup, RequiredMutuallyExclusiveOptionGroup
import click
import os
from dreem.demultiplexing.demultiplexing import demultiplex

@click.command()
#@optgroup.group('Files and folders paths')
@click.option('--out_dir', '-o', default=os.getcwd(), type=click.Path(), help='Where to output files')
@click.option('--fastq1', '-fq1', help='Paths to the fastq1 file (forward primer). Enter multiple times for multiple files', type=click.Path(exists=True), required=True)
@click.option('--fastq2', '-fq2', help='Paths to the fastq2 file (reverse primer). Enter multiple times for multiple files', type=click.Path(exists=True))
@click.option('--library', '-l', type=click.Path(exists=True), help='Path to the library.csv file')

#@optgroup.group('Demultiplexing parameters')
#@optgroup.option('--barcode_start', '-bs', type=int, help='Start position of the barcode in the read')
#@optgroup.option('--barcode_end', '-be', type=int, help='End position of the barcode in the read')


def run(**args):
    """Run the demultiplexing pipeline.

    Demultiplexes the reads and outputs one fastq file per construct in the directory `output_path`, using `temp_path` as a temp directory.

    Parameters from args:
    -----------------------
    library: str
        Path to the library file. Columns are (non-excusively): ['construct', 'barcode_start', 'barcode']
    fastq1: str
        Path to the FASTQ file or list of paths to the FASTQ files, forward primer.
    fastq2: str
        Path to the FASTQ file or list of paths to the FASTQ files, reverse primer.
    out_dir: str
        Name of the output directory.

    Returns
    -------
    1 if successful, 0 otherwise.

    """
    # Get the paths
    root = args['out_dir']
    temp_folder = os.path.join(root,'temp','demultiplexing')
    output_folder = os.path.join(root,'output','demultiplexing')
    fastq1 = args['fastq1']
    fastq2 = args['fastq2']

    # Load the library
    library = pd.read_csv(args['library'])[["construct", "barcode_start", "barcode_end", "barcode"]].dropna()

    # Make the folders
    util.make_folder(output_folder)
    util.make_folder(temp_folder)

    # Demultiplex
    for f1 in fastq1:
        for f2 in fastq2:
            if f1[:-len('_R1.fastq')] == f2[:-len('_R2.fastq')]:
                sample = f1.split('/')[-1][:-len('_R1.fastq')]
                util.make_folder(os.path.join(output_folder, sample))
                assert demultiplex(f1, f2, library, os.path.join(output_folder, sample), os.path.join(temp_folder, sample)), "Demultiplexing failed"
    return 1


if __name__ == '__main__':
    run()