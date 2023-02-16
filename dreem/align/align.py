import itertools
import logging
from collections import defaultdict
from multiprocessing import Pool

from dreem.util.seq import FastaParser, FastaWriter
from dreem.align.reads import (FastqAligner, FastqTrimmer, FastqUnit,
                               BamAlignSorter, BamSplitter,
                               SamRemoveEqualMappers)
from dreem.util.stargs import starstarmap
from dreem.util import path
from dreem.util.util import get_num_parallel


def check_for_duplicates(fq_units: list[FastqUnit]):
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
    return dups


def write_temp_ref_files(temp_path: path.TopDirPath,
                         refset_file: path.RefsetSeqInFilePath,
                         fq_units: list[FastqUnit]):
    """ Write temporary FASTA files, each containing one reference that
    corresponds to a FASTQ file from demultiplexing. """
    ref_paths: dict[str, path.OneRefSeqStepFilePath] = dict()
    # Only the reference sequences of FASTQ files that have come from
    # demultiplexing need to be written.
    refs = {fq_unit.ref for fq_unit in fq_units if fq_unit.demult}
    if refs:
        # Parse the FASTA only if there are any references to write.
        for ref, seq in FastaParser(refset_file.path).parse():
            if ref in refs:
                # Write the reference sequence to a temporary FASTA file
                # only if at least one demultiplexed FASTQ file uses it.
                ref_path = path.OneRefSeqStepFilePath(
                    top=temp_path.top,
                    module=path.Module.ALIGN,
                    step=path.Step.ALIGN_ALIGN,
                    ref=ref,
                    ext=path.FASTA_EXTS[0])
                ref_path.path.parent.mkdir(parents=True, exist_ok=True)
                logging.info(f"Writing temporary FASTA file: {ref_path.path}")
                FastaWriter(ref_path.path, {ref: seq}).write()
                ref_paths[ref] = ref_path
    return ref_paths


def run_steps_fq(out_path: path.TopDirPath,
                 temp_path: path.TopDirPath,
                 num_cpus: int,
                 fasta: path.RefsetSeqInFilePath | path.OneRefSeqStepFilePath,
                 fastq: FastqUnit,
                 *,
                 trim: bool,
                 trim_minq1: int,
                 trim_minq2: int,
                 trim_adapt15: str,
                 trim_adapt13: str,
                 trim_adapt25: str,
                 trim_adapt23: str,
                 trim_minover: int,
                 trim_maxerr: float,
                 trim_indels: bool,
                 trim_nextseq: bool,
                 trim_discard_trimmed: bool,
                 trim_discard_untrimmed: bool,
                 trim_minlen: int,
                 align_local: bool,
                 align_unal: bool,
                 align_disc: bool,
                 align_mixed: bool,
                 align_dove: bool,
                 align_cont: bool,
                 align_score: str,
                 align_minl: int,
                 align_maxl: int,
                 align_gbar: int,
                 align_slen: int,
                 align_sint: str,
                 align_exten: int,
                 align_reseed: int,
                 align_pad: int,
                 align_orient: str):
    # Trim the FASTQ file(s).
    trimmer = FastqTrimmer(temp_path, num_cpus, fastq)
    if trim:
        fastq = trimmer.run(trim_minq1=trim_minq1,
                            trim_minq2=trim_minq2,
                            trim_adapt15=trim_adapt15,
                            trim_adapt13=trim_adapt13,
                            trim_adapt25=trim_adapt25,
                            trim_adapt23=trim_adapt23,
                            trim_minover=trim_minover,
                            trim_maxerr=trim_maxerr,
                            trim_indels=trim_indels,
                            trim_nextseq=trim_nextseq,
                            trim_discard_trimmed=trim_discard_trimmed,
                            trim_discard_untrimmed=trim_discard_untrimmed,
                            trim_minlen=trim_minlen)
    # Align the FASTQ to the reference.
    aligner = FastqAligner(temp_path, num_cpus, fastq, fasta)
    xam_path = aligner.run(align_local=align_local,
                           align_unal=align_unal,
                           align_disc=align_disc,
                           align_mixed=align_mixed,
                           align_dove=align_dove,
                           align_cont=align_cont,
                           align_score=align_score,
                           align_minl=align_minl,
                           align_maxl=align_maxl,
                           align_gbar=align_gbar,
                           align_slen=align_slen,
                           align_sint=align_sint,
                           align_exten=align_exten,
                           align_reseed=align_reseed,
                           align_pad=align_pad,
                           align_orient=align_orient)
    trimmer.clean()
    # Remove equally mapping reads.
    remover = SamRemoveEqualMappers(temp_path, num_cpus, xam_path)
    xam_path = remover.run()
    aligner.clean()
    # Sort the SAM file and output a BAM file.
    sorter = BamAlignSorter(temp_path, num_cpus, xam_path)
    xam_path = sorter.run()
    remover.clean()
    # Split the BAM file into one file for each reference.
    splitter = BamSplitter(out_path, num_cpus, xam_path, fasta)
    bam_paths = splitter.run()
    sorter.clean()
    return bam_paths


def run_steps_fqs(out_dir: str,
                  temp_dir: str,
                  refset_file: str,
                  fq_units: list[FastqUnit],
                  parallel: bool,
                  max_procs: int,
                  **kwargs) -> tuple[path.OneRefAlignmentInFilePath, ...]:
    n_fqs = len(fq_units)
    if n_fqs == 0:
        logging.critical("No FASTQ files were given for alignment.")
        return ()
    # Confirm that there are no duplicate samples and references.
    if dups := check_for_duplicates(fq_units):
        logging.critical(f"Got duplicate samples/refs: {dups}")
        return ()
    if max_procs < 1:
        logging.warning("Maximum CPUs must be ≥ 1: setting to 1")
        max_procs = 1
    # Generate the paths.
    out_path = path.TopDirPath.parse_path(out_dir)
    temp_path = path.TopDirPath.parse_path(temp_dir)
    refset_path = path.RefsetSeqInFilePath.parse_path(refset_file)
    # Write the temporary FASTA files for demultiplexed FASTQs.
    ref_paths = write_temp_ref_files(temp_path, refset_path, fq_units)
    try:
        # Determine how to parallelize the tasks.
        n_tasks_parallel, n_procs_per_task = get_num_parallel(n_fqs,
                                                              max_procs,
                                                              parallel,
                                                              hybrid=True)
        # List the arguments for each task, including n_procs_per_task.
        align_args = list()
        align_kwargs = list()
        for fq in fq_units:
            if fq.demult:
                # If the FASTQ came from demultiplexing (so contains
                # reads from only one reference), then align to the
                # FASTA of only that reference.
                try:
                    fasta = ref_paths[fq.ref]
                except KeyError:
                    # If the FASTA with that reference does not exist,
                    # then log an error and skip this FASTQ.
                    logging.error(f"Skipping {', '.join(fq.pathstrs)} "
                                  f"because reference '{fq.ref}' not found "
                                  f"in FASTA file '{refset_path.path}'")
                    continue
            else:
                # If the FASTQ may contain reads from ≥ 2 references,
                # then align to the FASTA file with all references.
                fasta = refset_path
            align_args.append((out_path, temp_path, n_procs_per_task, fasta, fq))
            align_kwargs.append(kwargs)
        if n_tasks_parallel > 1:
            # Process multiple FASTQs simultaneously.
            with Pool(n_tasks_parallel) as pool:
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
        for ref_file in ref_paths.values():
            logging.info(f"Deleting temporary reference file: {ref_file.path}")
            ref_file.path.unlink(missing_ok=True)
    # Return a tuple of the final alignment map files.
    return bams
