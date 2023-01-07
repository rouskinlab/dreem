import pathlib
import sys,os
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from dreem.util.cli import FASTQ2, DEFAULT_DEMULTIPLEXED


from dreem.util.fq import get_fastq_name
from dreem.align.align import align_pipeline, demultiplexed_fun


def run(out_dir: str, fasta: str, fastq: str, fastq2: str=FASTQ2,
        demultiplexed: bool=DEFAULT_DEMULTIPLEXED, **kwargs):
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
    out_dir: str
        Path to the output folder (in general the sample). 
    fasta: str
        Path to the reference FASTA file.
    fastq1: str
        Path to the FASTQ file or list of paths to the FASTQ files, forward primer.
    fastq2: str
        Path to the FASTQ file or list of paths to the FASTQ files, reverse primer.
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
    interleaved: bool
        Whether the FASTQ files are interleaved (default: False).
    
    """
    if demultiplexed:
        sample = fastq.split('/')[-2] # pathlib.PosixPath(fastq).stem
        demultiplexed_fun(out_dir, fasta, sample, fastq, fastq2=fastq2,
                            **kwargs)
    else:
        sample = get_fastq_name(fastq, fastq2)
        align_pipeline(out_dir, fasta, sample, fastq, fastq2=fastq2,
                       **kwargs)
