import pandas as pd
import util

from click_option_group import optgroup
import click

@click.command()
@optgroup.group('Files and folders paths')
@optgroup.option('--root_dir', '-rd', default='', type=click.Path(exists=True), help='Where to output files and temp files')
@optgroup.option('--fasta', '-fa', type=click.Path(exists=True), help='Path to the fasta file', required=True)
@optgroup.option('--fastq1', '-fq1', help='Paths to the fastq1 file (forward primer). Enter multiple times for multiple files', multiple=True, type=click.Path(exists=True), required=True)
@optgroup.option('--fastq2', '-fq2', help='Paths to the fastq2 file (reverse primer). Enter multiple times for multiple files', multiple=True, type=click.Path(exists=True))
@optgroup.option('--samples', '-s', type=click.Path(exists=True), help='Path to the samples.csv file')
@optgroup.option('--library', '-l', type=click.Path(exists=True), help='Path to the library.csv file')

@optgroup.group('Demultiplexing parameters')
@optgroup.option('--demultiplexing', '-dx', type=bool, help='Use demultiplexing', default=False)
@optgroup.option('--barcode_start', '-bs', type=int, help='Start position of the barcode in the read')
@optgroup.option('--barcode_end', '-be', type=int, help='End position of the barcode in the read')


def run(args):
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
    paths = util.get_folders(args)
    temp_folder = paths['demultiplexing']['temp']
    output_folder = paths['demultiplexing']['output']
    fastq1 = args['fastq1']
    fastq2 = args['fastq2']

    # Load the library
    library = pd.read_csv(args['library'])[["construct", "barcode_start", "barcode_end", "barcode_sequence"]].dropna()

    # Demultiplex

    ### TODO: Implement demultiplexing


    ### TODO: Save the results

    # /output_folder/
    # —| /{sample_1}
    #     —| {construct_1}_1.fastq
    #     —| {construct_1}_2.fastq
    #     —| {construct_2.}_1.fastq
    #     …
    # —| /{sample_2}
    #     …

    return 1


if __name__ == '__main__':
    run()