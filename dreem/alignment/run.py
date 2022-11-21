import dreem.util
import os
import click

@click.command()
@click.option('--output', '-o', default=os.getcwd(), type=click.Path(exists=True), help='Where to output files')
@click.option('--sub_dir', '-sd', type=click.Path(), help='Name for a sub-directory for the output files, for example to group the constructs by sample', default=None)
@click.option('--fasta', '-fa', type=click.Path(exists=True), help='Path to the fasta file', required=True)
@click.option('--fastq1', '-fq1', help='Paths to the fastq1 file (forward primer). Enter multiple times for multiple files', type=click.Path(exists=True), required=True)
@click.option('--fastq2', '-fq2', help='Paths to the fastq2 file (reverse primer). Enter multiple times for multiple files', type=click.Path(exists=True))

def run(**args):
    """Run the alignment module.

    Parameters from args:
    -----------------------
    fasta: str
        Path to the reference FASTA file.
    fastq1: str
        Path to the FASTQ file or list of paths to the FASTQ files, forward primer.
    fastq2: str
        Path to the FASTQ file or list of paths to the FASTQ files, reverse primer.
    output: str
        Path to the output folder (the sample).
    sub_dir: str
        Name for a sub-directory for the output files, for example to group the constructs by sample.

    Returns
    -------
    1 if successful, 0 otherwise.

    """
    # Extract the arguments
    fasta = args['fasta']
    fastq1 = args['fastq1']
    fastq2 = args['fastq2']
    root = args['output']
    temp_folder = os.path.join(root,'temp','alignment', args['sub_dir'] if args['sub_dir'] else '')
    output_folder = os.path.join(root,'output','alignment', args['sub_dir'] if args['sub_dir'] else '')
    print(temp_folder)

    # Make folders
    dreem.util.make_folder(output_folder)
    dreem.util.make_folder(temp_folder)

    # TODO Implement alignment


    # TODO Save the results like this
    # /output_folder/
    #    —| {construct_1}.bam 
    #    —| {construct_2.}.bam
    #    —| ...
    # REPLACE THE CODE BELOW
    constructs = ['mttr-6-alt-h3']
    for construct in constructs:
        bam_content = f"""Just a placeholder."""
        dreem.util.run_cmd("echo {} > {}".format(bam_content, os.path.join(output_folder, f"{construct}.bam")))

    return 1

if __name__ == '__main__':
    run()