from collections import defaultdict
from functools import partial
import itertools
import logging
from multiprocessing import Pool

from ..align.reads import (BamAlignSorter, BamOutputter, BamSplitter,
                           FastqAligner, FastqTrimmer, FastqUnit,
                           SamRemoveEqualMappers)
from ..util import path
from ..util.docdef import autodef, autodoc
from ..util.seq import FastaParser, FastaWriter
from ..util.util import get_num_parallel


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
                         refset_path: path.RefsetSeqInFilePath,
                         fq_units: list[FastqUnit]):
    """ Write temporary FASTA files, each containing one reference that
    corresponds to a FASTQ file from demultiplexing. """
    ref_paths: dict[str, path.OneRefSeqStepFilePath] = dict()
    # Only the reference sequences of FASTQ files that have come from
    # demultiplexing need to be written.
    refs = {fq_unit.ref for fq_unit in fq_units if fq_unit.demult}
    if refs:
        # Parse the FASTA only if there are any references to write.
        for ref, seq in FastaParser(refset_path.path).parse():
            if ref in refs:
                # Write the reference sequence to a temporary FASTA file
                # only if at least one demultiplexed FASTQ file uses it.
                ref_path = path.OneRefSeqStepFilePath(
                    top=temp_path.top,
                    module=path.Module.ALIGN.value,
                    step=path.Step.ALIGN_ALIGN.value,
                    ref=ref,
                    ext=path.FASTA_EXTS[0])
                ref_path.path.parent.mkdir(parents=True, exist_ok=True)
                logging.info(f"Writing temporary FASTA file: {ref_path}")
                FastaWriter(ref_path.path, {ref: seq}).write()
                ref_paths[ref] = ref_path
    return ref_paths


def infer_outputs(out_path: path.TopDirPath,
                  fasta: path.RefsetSeqInFilePath | path.OneRefSeqStepFilePath,
                  fastq: FastqUnit):
    translator = path.AlignmentInToAlignmentOut.inverse()
    return [
        translator.trans_inst(path.OneRefAlignmentOutFilePath(
            top=out_path.top,
            module=path.Module.ALIGN,
            sample=fastq.sample,
            ref=ref,
            ext=path.BAM_EXT))
        for ref, _ in FastaParser(fasta.path).parse()
    ]


@autodoc(ret_doc="List of paths to binary alignment map (BAM) files")
@autodef()
def run_steps_fq(out_dir: path.TopDirPath,
                 temp_dir: path.TopDirPath,
                 save_temp: bool,
                 num_cpus: int,
                 fasta: path.RefsetSeqInFilePath | path.OneRefSeqStepFilePath,
                 fastq: FastqUnit,
                 /, *,
                 rerun: bool,
                 resume: bool,
                 fastqc: bool,
                 fastqc_extract: bool,
                 cut: bool,
                 cut_q1: int,
                 cut_q2: int,
                 cut_g1: str,
                 cut_a1: str,
                 cut_g2: str,
                 cut_a2: str,
                 cut_o: int,
                 cut_e: float,
                 cut_indels: bool,
                 cut_nextseq: bool,
                 cut_discard_trimmed: bool,
                 cut_discard_untrimmed: bool,
                 cut_m: int,
                 bt2_local: bool,
                 bt2_discordant: bool,
                 bt2_mixed: bool,
                 bt2_dovetail: bool,
                 bt2_contain: bool,
                 bt2_unal: bool,
                 bt2_score_min: str,
                 bt2_i: int,
                 bt2_x: int,
                 bt2_gbar: int,
                 bt2_l: int,
                 bt2_s: str,
                 bt2_d: int,
                 bt2_r: int,
                 bt2_dpad: int,
                 bt2_orient: str,
                 rem_buffer: int) -> list[path.OneRefAlignmentInFilePath]:
    """ Run all steps of alignment for one FASTQ file or one pair of
    mated FASTQ files. """
    # Determine whether alignment needs to be run.
    if not rerun:
        # If alignment is not required to be rerun, then check if all
        # the expected output files already exist.
        outputs = infer_outputs(out_dir, fasta, fastq)
        if all(output.path.is_file() for output in outputs):
            # If all the output files already exist, just return them.
            logging.warning(
                f"Skipping alignment for {' and '.join(fastq.str_paths)} "
                "because all expected output files already exist. "
                "Add --rerun to rerun.")
            return outputs
    # Trim the FASTQ file(s).
    trimmer = FastqTrimmer(top_dir=temp_dir, num_cpus=num_cpus, fastq=fastq,
                           save_temp=save_temp, resume=resume)
    if fastqc:
        trimmer.qc_input(fastqc_extract)
    if cut:
        if resume and all(p.is_file() for p in trimmer.output.paths):
            fastq = trimmer.output
        else:
            fastq = trimmer.run(cut_q1=cut_q1,
                                cut_q2=cut_q2,
                                cut_g1=cut_g1,
                                cut_a1=cut_a1,
                                cut_g2=cut_g2,
                                cut_a2=cut_a2,
                                cut_o=cut_o,
                                cut_e=cut_e,
                                cut_indels=cut_indels,
                                cut_nextseq=cut_nextseq,
                                cut_discard_trimmed=cut_discard_trimmed,
                                cut_discard_untrimmed=cut_discard_untrimmed,
                                cut_m=cut_m)
            if fastqc:
                trimmer.qc_output(fastqc_extract)
    # Align the FASTQ to the reference.
    aligner = FastqAligner(top_dir=temp_dir, num_cpus=num_cpus, fastq=fastq,
                           fasta=fasta, save_temp=save_temp, resume=resume)
    xam_path = aligner.run(bt2_local=bt2_local,
                           bt2_discordant=bt2_discordant,
                           bt2_mixed=bt2_mixed,
                           bt2_dovetail=bt2_dovetail,
                           bt2_contain=bt2_contain,
                           bt2_unal=bt2_unal,
                           bt2_score_min=bt2_score_min,
                           bt2_i=bt2_i,
                           bt2_x=bt2_x,
                           bt2_gbar=bt2_gbar,
                           bt2_l=bt2_l,
                           bt2_s=bt2_s,
                           bt2_d=bt2_d,
                           bt2_r=bt2_r,
                           bt2_dpad=bt2_dpad,
                           bt2_orient=bt2_orient)
    trimmer.clean()
    # Remove equally mapping reads.
    remover = SamRemoveEqualMappers(top_dir=temp_dir,
                                    num_cpus=num_cpus,
                                    xam=xam_path,
                                    save_temp=save_temp,
                                    resume=resume)
    xam_path = remover.run(rem_buffer=rem_buffer)
    aligner.clean()
    # Sort the SAM file and output a BAM file.
    sorter = BamAlignSorter(top_dir=temp_dir,
                            num_cpus=num_cpus,
                            xam=xam_path,
                            save_temp=save_temp,
                            resume=resume)
    xam_path = sorter.run()
    remover.clean()
    # Split the BAM file into one file for each reference.
    splitter = BamSplitter(top_dir=temp_dir,
                           num_cpus=num_cpus,
                           xam=xam_path,
                           fasta=fasta,
                           save_temp=save_temp,
                           resume=resume)
    bams = splitter.run()
    sorter.clean()
    # Output the BAM files and generate an index for each.
    bams = [BamOutputter(top_dir=out_dir,
                         num_cpus=num_cpus,
                         xam=bam,
                         save_temp=save_temp,
                         resume=resume).run()
            for bam in bams]
    splitter.clean()
    return bams


@autodoc()
@autodef()
def run_steps_fqs(fasta: str,
                  fq_units: list[FastqUnit],
                  /, *,
                  out_dir: str,
                  temp_dir: str,
                  save_temp: bool,
                  parallel: bool,
                  max_procs: int,
                  **kwargs) -> tuple[path.OneRefAlignmentInFilePath, ...]:
    """ Run all steps of alignment for one or more FASTQ files or pairs
    of mated FASTQ files. """
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
    # Generate the paths for the output, temporary files, and the input
    # FASTA file.
    out_path = path.TopDirPath.parse(out_dir)
    temp_path = path.TopDirPath.parse(temp_dir)
    refset_path = path.RefsetSeqInFilePath.parse(fasta)
    # Write the temporary FASTA files for demultiplexed FASTQs.
    temp_ref_paths = write_temp_ref_files(temp_path, refset_path, fq_units)
    try:
        # Determine how to parallelize each alignment task.
        n_tasks_parallel, procs_per_task = get_num_parallel(n_fqs,
                                                            max_procs,
                                                            parallel,
                                                            hybrid=True)
        # One alignment task will be created for each FASTQ unit.
        # Get the arguments for each task, including procs_per_task.
        iter_args = list()
        for fq in fq_units:
            if fq.demult:
                # If the FASTQ came from demultiplexing (so contains
                # reads from only one reference), then align to the
                # FASTA of only that reference.
                try:
                    fasta_path = temp_ref_paths[fq.ref]
                except KeyError:
                    # If the FASTA with that reference does not exist,
                    # then log an error and skip this FASTQ.
                    logging.error(f"Skipping {', '.join(fq.str_paths)} because "
                                  f"reference '{fq.ref}' was not found in "
                                  f"FASTA file: {refset_path}")
                    continue
            else:
                # If the FASTQ may contain reads from ≥ 1 references,
                # then align to the FASTA file with all references.
                fasta_path = refset_path
            # Add these arguments to the lists of arguments that will be
            # passed to run_steps_fq.
            iter_args.append((out_path, temp_path, save_temp, procs_per_task,
                              fasta_path, fq))
        partial_run_steps_fq = partial(run_steps_fq, **kwargs)
        if n_tasks_parallel > 1:
            # Process multiple FASTQ files simultaneously.
            with Pool(n_tasks_parallel) as pool:
                bams = pool.starmap(partial_run_steps_fq, iter_args)
        else:
            # Process the FASTQ files sequentially.
            bams = itertools.starmap(partial_run_steps_fq, iter_args)
    finally:
        if not save_temp:
            # Delete the temporary files before exiting.
            for ref_file in temp_ref_paths.values():
                logging.debug(
                    f"Deleting temporary reference file: {ref_file}")
                ref_file.path.unlink(missing_ok=True)
    # Return a tuple of the final alignment map files.
    return tuple(itertools.chain(*bams))
