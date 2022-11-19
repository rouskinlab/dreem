import util



def run(fasta, fastq1, fastq2, output_folder, temp_folder):
    """Run the alignment pipeline.
    
    fasta: str
        Path to the reference FASTA file.
    fastq1: str
        Path to the FASTQ file or list of paths to the FASTQ files, forward primer.
    fastq2: str
        Path to the FASTQ file or list of paths to the FASTQ files, reverse primer.
    output_folder: str
        Path to the output folder.
    temp_folder: str
        Path to the temporary folder.

    Returns
    -------
    1 if successful, 0 otherwise.

    """
    
    # TODO Implement alignment

    # TODO Save the results
    # /output_folder/
    #    —| {construct_1}.bam
    #    —| {construct_2.}.bam
    #    —| ...

    return 1