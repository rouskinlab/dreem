import logging
from multiprocessing import Pool

from dreem.util.cli import DEFAULT_NEXTSEQ_TRIM
from dreem.util.dflt import NUM_PROCESSES
from dreem.util.seq import FastaParser, FastaWriter
from dreem.util.reads import (FastqAligner, FastqTrimmer, FastqUnit,
                              BamAlignSorter, BamSplitter,
                              SamRemoveEqualMappers,
                              get_demultiplexed_fastq_files,
                              get_demultiplexed_fastq_pairs)
from dreem.util.stargs import starstarmap
from dreem.util import path
from dreem.util.seq import DNA


def _align(top_dir: path.TopDirPath,
           fasta: path.RefsetSeqInFilePath | path.OneRefSeqTempFilePath,
           fastq: FastqUnit,
           nextseq_trim: bool = DEFAULT_NEXTSEQ_TRIM):
    # Trim the FASTQ file(s).
    trimmer = FastqTrimmer(top_dir, fastq)
    fastq = trimmer.run(nextseq_trim=nextseq_trim)
    # Align the FASTQ to the reference.
    aligner = FastqAligner(top_dir, fastq, fasta)
    xam_path = aligner.run()
    trimmer.clean()
    # Remove equally mapping reads.
    remover = SamRemoveEqualMappers(top_dir, xam_path)
    xam_path = remover.run()
    aligner.clean()
    # Sort the SAM file and output a BAM file.
    sorter = BamAlignSorter(top_dir, xam_path)
    xam_path = sorter.run()
    remover.clean()
    # Split the BAM file into one file for each reference.
    splitter = BamSplitter(top_dir, xam_path, fasta)
    bams = splitter.run()
    sorter.clean()
    return bams


def all_refs(top_dir: str, fasta: str,
             fastqs: str, fastqi: str, fastq1: str, fastq2: str,
             phred_enc: int, **kwargs):
    fqs = path.SampleReadsInFilePath.parse_path(fastqs) if fastqs else None
    fqi = path.SampleReadsInFilePath.parse_path(fastqi) if fastqi else None
    fq1 = path.SampleReads1InFilePath.parse_path(fastq1) if fastq1 else None
    fq2 = path.SampleReads2InFilePath.parse_path(fastq2) if fastq2 else None
    return _align(path.TopDirPath.parse_path(top_dir),
                  path.RefsetSeqInFilePath.parse_path(fasta),
                  FastqUnit.wrap(fastqs=fqs, fastqi=fqi,
                                 fastq1=fq1, fastq2=fq2,
                                 phred_enc=phred_enc),
                  **kwargs)


def _get_fq_units_from_dir(fastqs_dir: str, fastqi_dir: str, fastq12_dir: str,
                           phred_enc: int):
    if (count := bool(fastqs_dir) + bool(fastqi_dir) + bool(fastq12_dir)) > 1:
        raise TypeError(f"Got {count} arguments for FASTQ dirs (expected 1)")
    if fastqs_dir:
        return get_demultiplexed_fastq_files(
            path.SampleInDirPath.parse_path(fastqs_dir), False, phred_enc)
    if fastqi_dir:
        return get_demultiplexed_fastq_files(
            path.SampleInDirPath.parse_path(fastqi_dir), True, phred_enc)
    if fastq12_dir:
        return get_demultiplexed_fastq_pairs(
            path.SampleInDirPath.parse_path(fastq12_dir), phred_enc)
    raise TypeError("Got no arguments for FASTQ dirs (expected 1)")


def _one_ref(top_dir: str, ref: str, seq: DNA, fq_unit: FastqUnit, **kwargs):
    # Write a temporary FASTA file for the reference.
    fasta = path.OneRefSeqTempFilePath(top=top_dir,
                                       partition=path.Partition.TEMP,
                                       module=path.Module.ALIGN,
                                       step=path.TempStep.ALIGN_ALIGN,
                                       ref=ref,
                                       ext=path.FASTA_EXTS[0])
    fasta.path.parent.mkdir(parents=True, exist_ok=True)
    try:
        FastaWriter(fasta.path, {ref: seq}).write()
        bams = _align(path.TopDirPath.parse_path(top_dir),
                      fasta, fq_unit, **kwargs)
        if len(bams) != 1:
            raise ValueError(f"Expected 1 BAM file, got {len(bams)}")
        return bams[0]
    finally:
        fasta.path.unlink(missing_ok=True)


def each_ref(top_dir: str, refs_file: str,
             fastqs_dir: str, fastqi_dir: str, fastq12_dir: str,
             phred_enc: int, **kwargs):
    fq_units = _get_fq_units_from_dir(fastqs_dir, fastqi_dir, fastq12_dir,
                                      phred_enc)
    align_args = list()
    align_kwargs = list()
    for ref, seq in FastaParser(refs_file).parse():
        try:
            fq_unit = fq_units[ref]
        except KeyError:
            logging.warning(f"No FASTQ files for reference '{ref}'")
        else:
            align_args.append((top_dir, ref, seq, fq_unit))
            align_kwargs.append(kwargs)
    if align_args:
        n_procs = min(len(align_args), NUM_PROCESSES)
        with Pool(n_procs) as pool:
            bams = list(starstarmap(pool.starmap, _one_ref,
                                    align_args, align_kwargs))
    else:
        bams = list()
    return bams
