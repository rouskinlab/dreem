import os
from multiprocessing import Pool

from dreem.util.cli import DEFAULT_INTERLEAVED_INPUT, DEFAULT_INTERLEAVE_OUTPUT, DEFAULT_TRIM, DEFAULT_NEXTSEQ_TRIM
from dreem.util.dflt import NUM_PROCESSES
from dreem.util.fa import FastaParser, FastaWriter
from dreem.util.fq import FastqAligner, FastqTrimmer, get_fastq_pairs, get_fastq_name
from dreem.util.star import starstarmap
from dreem.util.path import TEMP_DIR, try_remove
from dreem.util.seq import DNA
from dreem.util.xam import SamSorter, SamRemoveEqualMappers, SamOutputter, SamSplitter


def align_pipeline(out_dir: str, ref_file: str, sample: str, fastq: str,
             fastq2: str, trim: bool=DEFAULT_TRIM,
             interleaved_in: bool=DEFAULT_INTERLEAVED_INPUT,
             interleave_out: bool=DEFAULT_INTERLEAVE_OUTPUT,
             nextseq_trim: bool=DEFAULT_NEXTSEQ_TRIM):
    paired = interleaved_in or fastq2
    # Trim the FASTQ file.
    if trim:
        trimmer = FastqTrimmer(out_dir, ref_file, sample, paired, fastq, fastq2,
                               interleave_out=interleave_out)
        trimmer.run(nextseq_trim=nextseq_trim)
        fqs = trimmer.outputs
    else:
        fqs = [fastq]
        if fastq2:
            fqs.append(fastq2)

    # Align the FASTQ to the reference.
    aligner = FastqAligner(out_dir, ref_file, sample, paired, *fqs)
    aligner.run()
    trimmer.clean()
    # Remove equally mapping reads.
    rmequal = SamRemoveEqualMappers(out_dir, ref_file, sample, paired,
                                    aligner.sam_out)
    rmequal.run()
    aligner.clean()
    # Sort the SAM file and output a BAM file.
    sorter = SamSorter(out_dir, ref_file, sample, rmequal.sam_out)
    sorter.run()
    rmequal.clean()
    # Split the BAM file into one file for each reference.
    splitter = SamSplitter(out_dir, ref_file, sample, sorter.bam_out)
    splitter.run()
    sorter.clean()
    # Move the BAM files to the final output directory.
    bams = list()
    for bam in splitter.bams_out:
        outputter = SamOutputter(out_dir, ref_file, sample, bam)
        outputter.run()
        bams.append(outputter.bam_out)
    splitter.clean()
    return bams


def _dmplex_ref(out_dir: str, ref: bytes, seq: DNA, sample: str, fastq1: str,
                fastq2: str, **kwargs):
    """Run the alignment module.

    Aligns the reads to the reference genome and outputs one bam file perconstruct in the directory `output_path`, using `temp_path` as a temp directory.

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
    os.makedirs(temp_fasta_dir, exist_ok=True)
    
   # try:
    FastaWriter(temp_fasta, {ref: seq}).write()
    align_pipeline(out_dir, temp_fasta, sample, fastq1, fastq2, **kwargs)

   # finally:
   #     try_remove(temp_fasta)
   #     try:
   #         os.rmdir(temp_fasta_dir)
   #     except OSError:
   #         pass


def demultiplexed_fun(out_dir: str, fasta: str, sample: str, fastq: str,
                  fastq2: str, **kwargs):
    fq_dir = os.path.dirname(fastq)
    if fastq2:
        if fastq2.replace('_R2.','_R1.') != fastq:
            raise ValueError("fastq1 and fastq2 must be equal")
        #pairs = get_fastq_pairs(fq_dir)
        pairs = {get_fastq_name(fastq)[:-len('_R1')].encode('ascii'): (fastq, fastq2)}

    else:
        pairs = {get_fastq_name(fq): (os.path.join(fq_dir, fq), "")
                 for fq in os.listdir(fq_dir)}
    align_args = list()
    align_kwargs = list()
    for ref, seq in FastaParser(fasta).parse():
        try:
            fq1, fq2 = pairs[ref]
        except KeyError:
            pass
        else:
            arg = (out_dir, ref, seq, sample, fq1, fq2)
            align_args.append(arg)
            align_kwargs.append(kwargs)
            
    assert len(align_args) == len(align_kwargs)
    if align_args:
        n_procs = min(len(align_args), NUM_PROCESSES)
        with Pool(n_procs) as pool:
            starstarmap(pool.starmap, _dmplex_ref, align_args, align_kwargs)
            