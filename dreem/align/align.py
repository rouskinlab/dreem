import os
from multiprocessing import Pool

from dreem.base.dflt import NUM_PROCESSES
from dreem.base.fa import FastaParser, FastaWriter
from dreem.base.fq import FastqAligner, FastqTrimmer, get_fastq_pairs, get_fastq_name
from dreem.base.path import TEMP_DIR, try_remove
from dreem.base.seq import DNA
from dreem.base.xam import SamSorter, SamRemoveEqualMappers, SamOutputter, SamSplitter


def align_pipeline(out_dir: str, ref_file: str, sample: str, fastq: str,
                   fastq2: str = "", interleaved: bool = False):
    paired = interleaved or fastq2
    # Trim the FASTQ file.
    trimmer = FastqTrimmer(out_dir, ref_file, sample, paired, fastq, fastq2)
    trimmed = trimmer.run()
    # Align the FASTQ to the reference.
    aligner = FastqAligner(out_dir, ref_file, sample, paired, *trimmed)
    aligned = aligner.run()
    trimmer.clean()
    # Remove equally mapping reads.
    remover = SamRemoveEqualMappers(out_dir, ref_file, sample, paired, aligned)
    removed = remover.run()
    aligner.clean()
    # Sort the SAM file and output a BAM file.
    sorter = SamSorter(out_dir, ref_file, sample, removed)
    sorted = sorter.run()
    remover.clean()
    # Split the BAM file into one file for each reference.
    splitter = SamSplitter(out_dir, ref_file, sample, sorted)
    splits = splitter.run()
    sorter.clean()
    # Move the BAM files to the final output directory.
    bams = list()
    for split in splits:
        outputter = SamOutputter(out_dir, ref_file, sample, split)
        bams.append(outputter.run())
    splitter.clean()
    return bams


def align_single_ref(out_dir: str, ref: bytes, seq: DNA, sample: str, fastq1: str, fastq2: str = "", interleaved: bool = False):
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
    temp_dir = os.path.join(out_dir, TEMP_DIR, "alignment")
    temp_fasta_dir = os.path.join(temp_dir, "fasta")
    temp_fasta = os.path.join(temp_fasta_dir, f"{ref.decode()}.fasta")
    try:
        FastaWriter(temp_fasta, {ref: seq}).write()
        align_pipeline(out_dir, temp_fasta, sample, fastq1, fastq2, interleaved)
    finally:
        try_remove(temp_fasta)
        try:
            os.rmdir(temp_fasta_dir)
        except OSError:
            pass


def align_demultiplexed(out_dir: str, fasta: str, sample: str, fastq: str, fastq2: str = ""):
    fq_dir = os.path.dirname(fastq)
    if fastq2:
        if fastq2 != fastq:
            raise ValueError("fastq1 and fastq2 must be equal")
        pairs = get_fastq_pairs(fq_dir)
    else:
        pairs = {get_fastq_name(fq): (os.path.join(fq_dir, fq), "")
                 for fq in os.listdir(fq_dir)}
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
        n_procs = min(len(args), NUM_PROCESSES)
        with Pool(n_procs) as pool:
            pool.starmap(align_single_ref, args)
