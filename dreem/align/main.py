from collections import Counter
import os, sys
from typing import Optional, List

from multiprocessing import Pool
from dreem.util.cli_args import OUT_DIR, FASTA, FASTQ1, FASTQ2, INTERLEAVED, DEMULTIPLEXING, COORDS, PRIMERS, FILL

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from dreem.util.util import FastaParser, FastaWriter, DEFAULT_PROCESSES, name_temp_file, try_remove,  TEMP_DIR, DNA, run_cmd
from fastq import FastqTrimmer, FastqAligner, SamSorter, get_fastq_name, get_fastq_pairs, SamRemoveEqualMappers, SamOutputter, SamSplitter


def _align(root_dir: str, ref_file: str, sample: str, fastq: str,
           fastq2: str = "", interleaved: bool = False):
    paired = interleaved or fastq2
    # Trim the FASTQ file.
    fq_trim = FastqTrimmer(root_dir, ref_file, sample, paired, fastq, fastq2)
    # Align the FASTQ to the reference.
    sam_aln = FastqAligner(root_dir, ref_file, sample, paired, *fq_trim.run())
    # Remove equally mapping reads.
    sam_rem = SamRemoveEqualMappers(root_dir, ref_file, sample, sam_aln.run(),
                                    paired)
    fq_trim.clean()
    # Sort the SAM file and output a BAM file.
    bam_sort = SamSorter(root_dir, ref_file, sample, sam_rem.run())
    sam_aln.clean()
    # Split the BAM file into one file for each reference.
    bams_split = SamSplitter(root_dir, ref_file, sample, bam_sort.run())
    sam_rem.clean()
    # Move the BAM files to the final output directory.
    bams_out = [SamOutputter(root_dir, ref_file, sample, bam)
                for bam in bams_split.run()]
    bam_sort.clean()
    bams = [bam_out.run() for bam_out in bams_out]
    bams_split.clean()
    return bams


def _align_demultiplexed(root_dir: str, ref: bytes, seq: DNA, sample: str, fastq1: str, fastq2: str = "", interleaved: bool = False):
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
    
    """

    # Write a temporary FASTA file for the reference.
    temp_dir = os.path.join(root_dir, TEMP_DIR, "alignment")
    temp_fasta_dir = os.path.join(temp_dir, "fasta")
    temp_fasta = os.path.join(temp_fasta_dir, f"{ref.decode()}.fasta")
    try:
        FastaWriter(temp_fasta, {ref: seq}).write()
        _align(root_dir, temp_fasta, sample, fastq1, fastq2, interleaved)
    finally:
        try_remove(temp_fasta)
        try:
            os.rmdir(temp_fasta_dir)
        except OSError:
            pass


def run(out_dir: str= OUT_DIR, fasta: str= FASTA, fastq1: str= FASTQ1, fastq2: str = FASTQ2, demultiplexed: bool= DEMULTIPLEXING, interleaved: bool= INTERLEAVED, verbose: bool= False):
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
    
    Returns
    -------
    1 if successful, 0 otherwise.

    """
    #fastq_names = list(map(get_fastq_name, fastq, fq2s_list))
    if demultiplexed:
        args = list()
        sample = os.path.basename(fastq1)
        if fastq2:
            if fastq2 != fastq1:
                raise ValueError("fastq1 and fastq2 must be equal")
            pairs = get_fastq_pairs(fastq1)
        else:
            pairs = {get_fastq_name(fq): (os.path.join(fastq1, fq), FASTQ2)
                     for fq in os.listdir(fastq1)}
        args = list()
        for ref, seq in FastaParser(fasta).parse():
            try:
                fq1, fq2 = pairs[ref]
            except KeyError:
                pass
            else:
                arg = (out_dir, ref, seq, sample, fq1, fq2)
                args.append(arg)
        if args:
            with Pool(DEFAULT_PROCESSES) as pool:
                pool.starmap(_align_demultiplexed, args)
    else:
        sample = get_fastq_name(fastq1, fastq2)
        _align(out_dir, fasta, sample, fastq1, fastq2)
