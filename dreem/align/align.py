import itertools
import logging
from collections import defaultdict
from multiprocessing import Pool

from dreem.util.cli import DEFAULT_NEXTSEQ_TRIM, ParallelChoice
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
                ref_path = ref_file.path
                ref_path.parent.mkdir(parents=True, exist_ok=True)
                logging.info(f"Writing temporary reference file: {ref_path}")
                FastaWriter(ref_path, {ref: seq}).write()
                ref_files[ref] = ref_file
    return ref_files


def run_steps_fq(top_dir: path.TopDirPath,
                 max_cpus: int,
                 fasta: path.RefsetSeqInFilePath | path.OneRefSeqTempFilePath,
                 fastq: FastqUnit,
                 nextseq_trim: bool = DEFAULT_NEXTSEQ_TRIM):
    # Trim the FASTQ file(s).
    trimmer = FastqTrimmer(top_dir, max_cpus, fastq)
    fastq = trimmer.run(nextseq_trim=nextseq_trim)
    # Align the FASTQ to the reference.
    aligner = FastqAligner(top_dir, max_cpus, fastq, fasta)
    xam_path = aligner.run()
    trimmer.clean()
    # Remove equally mapping reads.
    remover = SamRemoveEqualMappers(top_dir, max_cpus, xam_path)
    xam_path = remover.run()
    aligner.clean()
    # Sort the SAM file and output a BAM file.
    sorter = BamAlignSorter(top_dir, max_cpus, xam_path)
    xam_path = sorter.run()
    remover.clean()
    # Split the BAM file into one file for each reference.
    splitter = BamSplitter(top_dir, max_cpus, xam_path, fasta)
    bams = splitter.run()
    sorter.clean()
    return bams


def run_steps_fqs(top_path: str,
                  refset_path: str,
                  fq_units: list[FastqUnit],
                  parallel: str,
                  max_cpus: int,
                  **kwargs):
    if not fq_units:
        raise ValueError(f"No FASTQ files given")
    n_fqs = len(fq_units)  # > 0
    # Confirm that there are no duplicate samples.
    confirm_no_duplicate_samples(fq_units)
    # Generate the paths.
    top_dir = path.TopDirPath.parse_path(top_path)
    refset_file = path.RefsetSeqInFilePath.parse_path(refset_path)
    # Write the temporary FASTA files for demultiplexed FASTQs.
    ref_files = write_temp_ref_files(top_dir, refset_file, fq_units)
    try:
        # Determine if the computation should be parallelized
        # - broadly (i.e. run all FASTQs simultaneously in parallel and limit
        #   Cutadapt and Bowtie2 to single-threaded mode),
        # - deeply (i.e. run the FASTQs sequentially with Cutadapt and Bowtie2
        #   running in multithreaded mode),
        # - or with parallelization off (i.e. run the FASTQs sequentially and
        #   limit Cutadapt and Bowtie2 to single-threaded mode).
        parallel_broad = (parallel == ParallelChoice.BROAD
                          or (parallel == ParallelChoice.AUTO and n_fqs > 1))
        parallel_deep = (parallel == ParallelChoice.DEEP
                         or (parallel == ParallelChoice.AUTO and n_fqs == 1))
        # Compute the maximum number of processes allowed per FASTQ.
        procs_per_fq = max((max_cpus // n_fqs if parallel_broad
                            else max_cpus if parallel_deep
                            else 1), 1)
        # Next determine the positional and keyword arguments for each FASTQ.
        align_args = list()
        align_kwargs = list()
        for fq_unit in fq_units:
            # If the FASTQ is demultiplexed, then align to the FASTA containing
            # the one reference corresponding to the FASTQ; otherwise, align to
            # the oringal FASTA file that may contain more than one sequence.
            if fq_unit.demult:
                try:
                    fasta = ref_files[fq_unit.ref]
                except KeyError:
                    raise ValueError(f"Reference '{fq_unit.ref}' not found in "
                                     f"FASTA file {refset_file.path}")
            else:
                fasta = refset_file
            align_args.append((top_dir, procs_per_fq, fasta, fq_unit))
            align_kwargs.append(kwargs)
            print("FASTQ Unit", fq_unit, align_args[-1], align_kwargs[-1])
        if parallel_broad:
            # Process all FASTQs simultaneously in parallel.
            n_procs = max(min(len(fq_units), max_cpus), 1)
            with Pool(n_procs) as pool:
                bams = tuple(itertools.chain(*starstarmap(pool.starmap,
                                                          run_steps_fq,
                                                          align_args,
                                                          align_kwargs)))
        else:
            # Process each FASTQ file sequentially.
            bams = tuple(itertools.chain(*starstarmap(itertools.starmap,
                                                      run_steps_fq,
                                                      align_args,
                                                      align_kwargs)))
    finally:
        # Always delete the temporary files before exiting.
        for ref_file in ref_files.values():
            logging.info(f"Deleting temporary reference file: {ref_file.path}")
            ref_file.path.unlink(missing_ok=True)
    # Return a tuple of the final alignment map files.
    return bams
