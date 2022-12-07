import os
from typing import Optional, List

import dreem
from dreem.util import FastaParser


def align_reads(construct, sequence, fastq1, fastq2, output_folder, temp_folder):
    """Run the alignment module.

    Aligns the reads to the reference genome and outputs one bam file per construct in the directory `output_path`, using `temp_path` as a temp directory.

    Parameters from args:
    -----------------------
    construct: str
        Name of the construct.
    sequence: str
        Sequence of the construct.    
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
    try:
        __align_reads(construct, sequence, fastq1, fastq2, output_folder, temp_folder)
    except Exception as e:
        print(e)
        return 0
        

def __align_reads(construct, sequence, fastq1, fastq2, output_folder, temp_folder):     

    # Files
    reference_fasta = create_ref_fasta(construct, sequence, temp_folder)
    bam_file = os.path.join(output_folder, construct + '.bam')
    sam_file = os.path.join(temp_folder, construct + '.sam')
    
    ## Create sam file

    #TODO

    # Convert the sam file to bam
    __sam_to_bam(sam_file, bam_file)


def create_ref_fasta(construct, sequence, temp_folder):
    reference_fasta = os.path.join(temp_folder, construct + '.fasta')
    with open(reference_fasta, 'w') as f:
        f.write('>' + construct + '\n')
        f.write(sequence + '\n')
        f.close()
    return reference_fasta


def __sam_to_bam(sam_file, bam_file):
    dreem.util.run_cmd("samtools view -bS {} > {}".format(sam_file, bam_file))



def run(root: str, fasta: str, fastq1: List[str], fastq2: Optional[List[str]], **kwargs):
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
    temp_folder = os.path.join(root, "temp", "alignment")
    output_folder = os.path.join(root, "output", "alignment")

    # Make folders
    dreem.util.make_folder(output_folder)
    dreem.util.make_folder(temp_folder)

    # Align to each reference in the FASTA file.
    for ref, seq in FastaParser(fasta):
        align_reads(construct, sequence, fastq1, fastq2, output_folder, temp_folder), 'Alignment failed for construct {}'.format(construct)

