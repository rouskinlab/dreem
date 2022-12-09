from collections import Counter
import os, sys
from typing import Optional, List

from multiprocessing import Pool

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import dreem
from dreem.util.util import FastaParser, FastaWriter, DEFAULT_PROCESSES, name_temp_file, try_remove,  TEMP_DIR, DNA, run_cmd
from fastq import FastqInterleaver, FastqTrimmer, FastqMasker, FastqAligner, SamSorter, get_fastq_name, get_fastq_pairs, SamRemoveEqualMappers, SamOutputter, SamSplitter


def _align(root_dir: str, ref_file: str, sample: str, fastq1: str,
           fastq2: Optional[str] = None, interleaved: bool = False):
    if interleaved:
        if fastq2:
            raise ValueError(
                "fastq2 cannot be given if fastq1 is interleaved.")
        paired = True
    else:
        if fastq2:
            paired = True
            fastq1 = FastqInterleaver(root_dir, ref_file, sample, fastq1,
                                     fastq2).run()
        else:
            paired = False
    # Trim the FASTQ file.
    fastq1 = FastqTrimmer(root_dir, ref_file, sample, fastq1, paired).run()
    # Mask any remaining low-quality bases with N.
    fastq1 = FastqMasker(root_dir, ref_file, sample, fastq1, paired).run()
    # Align the FASTQ to the reference.
    sam = FastqAligner(root_dir, ref_file, sample, fastq1, paired).run()
    # Remove equally mapping reads.
    sam = SamRemoveEqualMappers(root_dir, ref_file, sample, sam, paired).run()
    # Sort the SAM file and output a BAM file.
    bam = SamSorter(root_dir, ref_file, sample, sam, paired).run()
    # Split the BAM file into one file for each reference.
    bam_dir = SamSplitter(root_dir, ref_file, sample, bam, paired).run()
    # Move the BAM files to the final output directory.
    bams = SamOutputter(root_dir, ref_file, sample, bam_dir, paired).run()


def _align_demultiplexed(root_dir: str, ref: bytes, seq: DNA, sample: str, fastq1: str, fastq2: Optional[str] = None, interleaved: bool = False):
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

    """

    # Write a temporary FASTA file for the reference.
    temp_dir = os.path.join(root_dir, TEMP_DIR, "alignment")
    temp_fasta = os.path.join(temp_dir, f"{ref.decode()}.fasta")
    try:
        FastaWriter(temp_fasta, {ref: seq}).write()
        _align(root_dir, temp_fasta, sample, fastq1, fastq2, interleaved)
    finally:
        try_remove(temp_fasta)


def run(out_dir: str, fasta: str, fastq1: List[str], fastq2: Optional[List[str]], demultiplexed: bool = False, **kwargs):
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
    coords: tuple
        coordinates for reference: '-c ref-name first last'
    primers: tuple
        primers for reference: '-p ref-name fwd rev'
    fill: bool
        Fill in coordinates of reference sequences for which neither coordinates nor primers were given (default: no).
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
    #fastq_names = list(map(get_fastq_name, fastq, fq2s_list))
    if demultiplexed:
        args = list()
        sample = os.path.basename(fastq1)
        if fastq2:
            if fastq2 != fastq1:
                raise ValueError("fastq1 and fastq2 must be equal")
            pairs = get_fastq_pairs(fastq1)
        else:
            pairs = {get_fastq_name(fq): (os.path.join(fastq1, fq), None)
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
        _align(out_dir, fasta, sample, fastq1, fastq2, **kwargs)
