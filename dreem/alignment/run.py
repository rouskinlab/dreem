import dreem.util
import os
import click
from dreem.alignment.alignment import align_reads

@click.command()
@click.option('--out_dir', '-o', default=os.getcwd(), type=click.Path(exists=True), help='Where to output files')
@click.option('--fasta', '-fa', type=click.Path(exists=True), help='Path to the fasta file', required=True)
@click.option('--fastq1', '-fq1', help='Paths to the fastq1 file (forward primer). Enter multiple times for multiple files', type=click.Path(exists=True), required=True)
@click.option('--fastq2', '-fq2', help='Paths to the fastq2 file (reverse primer). Enter multiple times for multiple files', type=click.Path(exists=True))

def run(**args):
    """Run the alignment module.

    Aligns the reads to the reference genome and outputs one bam file per construct in the directory `output_path`, using `temp_path` as a temp directory.

     /out_dir/
        —| {construct_1}.bam 
        —| {construct_2}.bam
        —| ...

    Parameters from args:
    -----------------------
    fasta: str
        Path to the reference FASTA file.
    fastq1: str
        Path to the FASTQ file or list of paths to the FASTQ files, forward primer.
    fastq2: str
        Path to the FASTQ file or list of paths to the FASTQ files, reverse primer.
    out_dir: str
        Path to the output folder (in general the sample).

    Returns
    -------
    1 if successful, 0 otherwise.

    """
    # Extract the arguments
    fasta = args['fasta']
    fastq1 = args['fastq1']
    fastq2 = args['fastq2']
    root = args['out_dir']
    temp_folder = os.path.join(root,'temp','alignment')
    output_folder = os.path.join(root,'output','alignment')

    # Make folders
    dreem.util.make_folder(output_folder)
    dreem.util.make_folder(temp_folder)

    # open the fasta file
    with open(fasta, 'r') as f:
        while line:= f.readline():
            if line.startswith('>'):
                construct = line[1:].strip()
                sequence = f.readline().strip()      
                assert align_reads(construct, sequence, fastq1, fastq2, output_folder, temp_folder), 'Alignment failed for construct {}'.format(construct)

if __name__ == '__main__':
    run()