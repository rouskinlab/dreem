from collections import Counter, defaultdict
from multiprocessing import Pool

from dreem.util.cli import DEFAULT_NEXTSEQ_TRIM
from dreem.util.dflt import NUM_PROCESSES
from dreem.util.seq import FastaParser, FastaWriter
from dreem.util.reads import (FastqAligner, FastqTrimmer, FastqUnit,
                              BamAlignSorter, BamSplitter,
                              SamRemoveEqualMappers)
from dreem.util.stargs import starstarmap
from dreem.util import path


def confirm_no_duplicate_samples(fq_units: list[FastqUnit]):
    # Count the number of times each sample and reference occurs.
    samples = defaultdict(int)
    sample_refs = defaultdict(lambda: defaultdict(int))
    for fq_unit in fq_units:
        if fq_unit.demult:
            sample_refs[fq_unit.sample][fq_unit.ref] += 1
        else:
            samples[fq_unit.sample] += 1
    # Find duplicates.
    dups = set()
    # Duplicate whole-sample FASTQs
    dups = dups | {sample for sample, count in samples.items() if count > 1}
    # Duplicate demultiplexed FASTQs
    dups = dups | {(sample, ref) for sample, refs in sample_refs.items()
                   for ref, count in refs.items() if count > 1}
    # Duplicate samples between whole-sample and demultiplexed FASTQs
    dups = dups | (set(samples) & set(sample_refs))
    if dups:
        # If there are any duplicate samples, raise an error.
        raise ValueError(f"Got duplicate samples/refs: {dups}")


def write_temp_ref_files(top_dir: path.TopDirPath,
                         refset_file: path.RefsetSeqInFilePath,
                         fq_units: list[FastqUnit]):
    ref_files: dict[str, path.OneRefSeqTempFilePath] = dict()
    # Determine which reference sequences need to be written.
    refs = {fq_unit.ref for fq_unit in fq_units if fq_unit.demult}
    if refs:
        # Only parse the FASTA if there are any references to write.
        for ref, seq in FastaParser(refset_file.path).parse():
            if ref in refs:
                ref_file = path.OneRefSeqTempFilePath(
                    top=top_dir.top,
                    partition=path.Partition.TEMP,
                    module=path.Module.ALIGN,
                    step=path.TempStep.ALIGN_ALIGN,
                    ref=ref,
                    ext=path.FASTA_EXTS[0])
                ref_file.path.parent.mkdir(parents=True, exist_ok=True)
                FastaWriter(ref_file.path, {ref: seq}).write()
                ref_files[ref] = ref_file
    return ref_files


def run_steps(top_dir: path.TopDirPath,
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


def run_steps_parallel(top_path: str,
                       refset_path: str,
                       fq_units: list[FastqUnit],
                       **kwargs):
    # Confirm that there are no duplicate samples.
    confirm_no_duplicate_samples(fq_units)
    # Generate the paths.
    top_dir = path.TopDirPath.parse_path(top_path)
    refset_file = path.RefsetSeqInFilePath.parse_path(refset_path)
    # Write the temporary FASTA files for demultiplexed FASTQs.
    ref_files = write_temp_ref_files(top_dir, refset_file, fq_units)
    try:
        align_args = list()
        align_kwargs = list()
        for fq_unit in fq_units:
            fasta = ref_files[fq_unit.ref] if fq_unit.demult else refset_file
            align_args.append((top_dir, fasta, fq_unit))
            align_kwargs.append(kwargs)
        if align_args:
            n_procs = min(len(align_args), NUM_PROCESSES)
            with Pool(n_procs) as pool:
                bams = tuple(starstarmap(pool.starmap, run_steps,
                                         align_args, align_kwargs))
        else:
            bams = tuple()
        return bams
    finally:
        # Always delete the temporary files before exiting.
        for ref_file in ref_files.values():
            ref_file.path.unlink(missing_ok=True)
