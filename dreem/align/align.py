import logging
import os
from multiprocessing import Pool
from typing import Optional

from dreem.util.cli import DEFAULT_INTERLEAVED_INPUT, DEFAULT_INTERLEAVE_OUTPUT, DEFAULT_TRIM, DEFAULT_NEXTSEQ_TRIM
from dreem.util.dflt import NUM_PROCESSES
from dreem.util.seq import FastaParser, FastaWriter
from dreem.util.reads import FastqAligner, FastqTrimmer, get_demultiplexed_fastq_pairs, FastqUnit, BamAlignSorter, SamRemoveEqualMappers, BamSplitter
from dreem.util.stargs import starstarmap
from dreem.util.path import TEMP_DIR, try_remove, BasePath, FastaSampleOutPath, MOD_ALN, FastaSegment, FastaInPath
from dreem.util.seq import DNA


def _align(base_path: BasePath,
           ref_file: FastaInPath | FastaSampleOutPath,
           fq_unit: FastqUnit,
           trim: bool=DEFAULT_TRIM,
           nextseq_trim: bool=DEFAULT_NEXTSEQ_TRIM):
    # Trim the FASTQ file(s).
    if trim:
        trimmer = FastqTrimmer(base_path, fq_unit)
        fq_unit = trimmer.run(nextseq_trim=nextseq_trim)
    # Align the FASTQ to the reference.
    aligner = FastqAligner(base_path, fq_unit, ref_file)
    xam_path = aligner.run()
    trimmer.clean()
    # Remove equally mapping reads.
    remover = SamRemoveEqualMappers(base_path, xam_path)
    xam_path = remover.run()
    aligner.clean()
    # Sort the SAM file and output a BAM file.
    sorter = BamAlignSorter(base_path, xam_path)
    xam_path = sorter.run()
    remover.clean()
    # Split the BAM file into one file for each reference.
    splitter = BamSplitter(base_path, xam_path, ref_file)
    bams = splitter.run()
    sorter.clean()
    return bams


def all_refs(base_dir: str, ref_file: str,
             fastqu: str, fastqi: str, fastq1: str, fastq2: str, **kwargs):
    base_path = BasePath.parse(base_dir)
    ref_path = FastaInPath.parse(ref_file)
    fq_unit = FastqUnit(fastqu=fastqu, fastqi=fastqi, fastq1=fastq1, fastq2=fastq2)
    return _align(base_path, ref_path, fq_unit, **kwargs)


def _get_fq_inputs(fastqu_dir: str, fastqi_dir: str, fastq12_dir: str):
    def assert_inputs_exist(exist: bool):
        if bool(fq_inputs) != exist:
            error = "no" if exist else "more than one"
            raise ValueError(f"Received {error} argument for FASTQ directory.")
    fq_inputs = dict()
    if fastqu_dir:
        assert_inputs_exist(False)
        # FIXME
    if fastqi_dir:
        assert_inputs_exist(False)
        # FIXME
    if fastq12_dir:
        assert_inputs_exist(False)
        fq_inputs = get_demultiplexed_fastq_pairs(fastq12_dir)
    assert_inputs_exist(True)
    return fq_inputs


def _one_ref(base_dir: str, ref: str, seq: DNA, fq_unit: FastqUnit, **kwargs):
    # Write a temporary FASTA file for the reference.
    base_path = BasePath.parse(base_dir)
    fasta_name = FastaSegment.format(ref)
    fasta = FastaSampleOutPath(base_dir, TEMP_DIR, MOD_ALN, fq_unit.sample,
                               fasta_name)
    fasta.path.parent.mkdir(parents=True, exist_ok=True)
    try:
        FastaWriter(fasta.path, {ref: seq}).write()
        bams = _align(base_path, fasta, fq_unit, **kwargs)
        assert len(bams) == 1
        return bams[0]
    finally:
        try_remove(fasta.path)


def each_ref(base_dir: str, refs_file: str,
             fastqu_dir: str, fastqi_dir: str, fastq12_dir: str, **kwargs):
    base_path = BasePath.parse(base_dir)
    fq_inputs = _get_fq_inputs(fastqu_dir, fastqi_dir, fastq12_dir)
    align_args = list()
    align_kwargs = list()
    for ref, seq in FastaParser(refs_file).parse():
        try:
            fq_unit = fq_inputs[ref]
        except KeyError:
            logging.warning(f"No FASTQ files for reference '{ref}'")
        else:
            align_args.append((base_path, ref, seq, fq_unit))
            align_kwargs.append(kwargs)
            
    assert len(align_args) == len(align_kwargs)
    if align_args:
        n_procs = min(len(align_args), NUM_PROCESSES)
        with Pool(n_procs) as pool:
            bams = starstarmap(pool.starmap, _one_ref, align_args, align_kwargs)
    return bams
