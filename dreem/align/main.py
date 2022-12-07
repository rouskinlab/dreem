import os
from typing import Optional, List

from multiprocessing import Pool

import dreem
from dreem.util import FastaParser, DEFAULT_PROCESSES


def align_demultiplexed(construct, sequence, fastq1, fastq2, output_folder, temp_folder):
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



def run(root: str, fasta: str, fastq1: List[str], fastq2: Optional[List[str]], demultiplexed: bool = False, **kwargs):
    """Run the alignment module.

    Aligns the reads to the reference genome and outputs one bam file per construct in the directory `output_path`, using `temp_path` as a temp directory.

     /out_dir/
        {sample_1}/
            —| {ref_1}.bam
            —| {ref_2}.bam
            —| ...
        {sample_2}/
        ...

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
    demultiplexed: bool
        Whether the FASTQ files were demultiplexed (default: False).
        If True:
            Assume that each FASTQ file contains reads from ONE sample and ONE reference.
            This happens after the (optional) demultiplexing step, whose output follows this structure:
            Every FASTQ file given must be named with the reference name, and if paired-end, followed by the mate (1 or 2).
            Every FASTQ file must be inside a directory with the sample name.
                {sample_1}/
                    |- {ref_1}_R1.fastq
                    |- {ref_1}_R2.fastq
                    |- {ref_2}_R1.fastq
                    |- {ref_2}_R2.fastq
                    |- ...
                {sample_2}/
                ...
            Since each FASTQ is known to correspond to only one reference (and may misalign if aligned to all references in the FASTA file),
            for each reference in the FASTA file, create a new FASTA file with only that reference and align all FASTQ files with the corresponding name to that reference.
        If False (default):
            Assume each FASTQ file contains reads from ONE sample and ONE OR MORE references.
            This is the typical structure of the FASTQ files directly from a sequencer.
            Every FASTQ file given must be named with the sample name, and if paired-end, followed by the mate (1 or 2).
            The directory of the FASTQ files does not matter.
                {sample_1}_R1.fastq
                {sample_1}_R2.fastq
                {sample_2}_R1.fastq
                {sample_2}_R2.fastq
                ...

    Returns
    -------
    1 if successful, 0 otherwise.

    """

    temp_folder = os.path.join(root, "temp", "alignment")
    output_folder = os.path.join(root, "output", "alignment")

    # Make folders
    dreem.util.make_folder(output_folder)
    dreem.util.make_folder(temp_folder)

    if demultiplexed:
        # Align to each reference in the FASTA file.
        args = ((ref, seq, fastq1, fastq2, output_folder, temp_folder)
                for ref, seq in FastaParser(fasta).parse())
        with Pool(DEFAULT_PROCESSES) as pool:
            pool.starmap(align_demultiplexed, args)
    else:
        align_combined()
