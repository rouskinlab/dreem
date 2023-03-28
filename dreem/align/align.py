from collections import defaultdict
from functools import partial
import itertools
import logging
from multiprocessing import Pool
import pathlib

from ..align.reads import (FastqUnit, run_fastqc, run_cutadapt, run_bowtie2,
                           dedup_sam, sort_xam, view_xam, index_bam_file,
                           index_fasta_file, get_fasta_index_paths)
from ..util import path
from ..util.seq import FastaParser, FastaWriter
from ..util.util import get_num_parallel

logger = logging.getLogger(__name__)


def deduplicate_fqs(fq_units: list[FastqUnit]):
    """ Return the FASTQ units with no duplicate samples or
    sample-reference combinations. """
    # Count the number of times each sample and reference occurs.
    sample_fqs: dict[str, FastqUnit] = dict()
    oneref_fqs: dict[str, dict[str, FastqUnit]] = defaultdict(lambda: dict())
    dedup_fqs: dict[str | tuple[str, str], FastqUnit] = dict()
    for fq_unit in fq_units:
        if fq_unit.sample in sample_fqs:
            # Another FASTQ unit comes from the same whole sample.
            logger.warning(f"Got >1 FASTQ from sample '{fq_unit.sample}'")
            try:
                # Remove the duplicate FASTQ unit.
                dedup_fqs.pop(fq_unit.sample)
            except KeyError:
                # It may have been removed already.
                pass
            # Skip this FASTQ unit.
            continue
        if fq_unit.ref:
            # The FASTQ unit comes from one reference in one sample.
            sample_ref = fq_unit.sample, fq_unit.ref
            if fq_unit.ref in oneref_fqs[fq_unit.sample]:
                # Another FASTQ unit comes from the same reference and
                # sample.
                logger.warning(
                    f"Got >1 FASTQ from sample and reference {sample_ref}")
                try:
                    # Remove the duplicate FASTQ unit.
                    dedup_fqs.pop(sample_ref)
                except KeyError:
                    # It may have been removed already.
                    pass
                # Skip this FASTQ unit.
                continue
            # Add the FASTQ to the full and deduplicated collections.
            oneref_fqs[fq_unit.sample][fq_unit.ref] = fq_unit
            dedup_fqs[sample_ref] = fq_unit
        else:
            # The FASTQ unit comes from an entire sample.
            if fq_unit.sample in oneref_fqs:
                # A demultiplexed FASTQ unit comes from the same sample.
                logger.warning(f"Got >1 FASTQ from sample '{fq_unit.sample}'")
                for ref in oneref_fqs[fq_unit.sample]:
                    try:
                        # Remove the duplicate FASTQ unit.
                        dedup_fqs.pop((fq_unit.sample, ref))
                    except KeyError:
                        # It may have been removed already.
                        pass
                    # Skip this FASTQ unit.
                    continue
            # Add the FASTQ to the full and deduplicated collections.
            sample_fqs[fq_unit.sample] = fq_unit
            dedup_fqs[fq_unit.sample] = fq_unit
    # Return the list of FASTQ units without duplicates.
    return list(dedup_fqs.values())


def write_temp_ref_files(temp_dir: pathlib.Path,
                         fasta: pathlib.Path,
                         refs: set[str]):
    """ Write temporary FASTA files, each containing one reference that
    corresponds to a FASTQ file from demultiplexing. """
    ref_paths: dict[str, pathlib.Path] = dict()
    if refs:
        # Parse the FASTA only if there are any references to write.
        for ref, seq in FastaParser(fasta).parse():
            if ref in refs:
                # Write the reference sequence to a temporary FASTA file
                # only if at least one demultiplexed FASTQ file uses it.
                ref_path = path.OneRefSeqTempFilePath(top=temp_dir,
                                                      module=path.Module.ALIGN.value,
                                                      step=path.Step.ALIGN_REFS.value,
                                                      ref=ref,
                                                      ext=path.FASTA_EXTS[0]).path
                # Create the parent directory.
                logger.debug(f"Creating directory: {ref_path.parent}")
                ref_path.parent.mkdir(parents=True, exist_ok=True)
                # Write the temporary FASTA file.
                logger.debug(f"Writing temporary FASTA file: {ref_path}")
                FastaWriter(ref_path, {ref: seq}).write()
                ref_paths[ref] = ref_path
    if missing := sorted(refs - set(ref_paths.keys())):
        logger.warning(f"Missing references in {fasta}: {', '.join(missing)}")
    return ref_paths


def index_temp_ref_file(fasta: str | pathlib.Path,
                        temp_dir: str | pathlib.Path,
                        n_procs: int):
    """ Build a temporary Bowtie2 index for a FASTA file. """
    prefix = path.RefsetBowtie2IndexTempPrefix(top=temp_dir,
                                               module=path.Module.ALIGN,
                                               step=path.Step.ALIGN_REFS).path
    index_fasta_file(fasta, prefix, n_procs)
    return prefix


def infer_outputs(out_dir: str | pathlib.Path,
                  sample: str,
                  refs: list[str]):
    """ Infer the paths of the BAM files that are expected to be output
    from the alignment. """
    return [path.OneRefAlignmentOutFilePath(top=str(out_dir),
                                            module=path.Module.ALIGN,
                                            sample=sample,
                                            ref=ref,
                                            ext=path.BAM_EXT).path
            for ref in refs]


def run_steps_fq(fq_inp: FastqUnit,
                 fasta: pathlib.Path,
                 bowtie2_index: pathlib.Path,
                 *,
                 n_procs: int,
                 out_dir: str,
                 temp_dir: str,
                 save_temp: bool,
                 rerun: bool,
                 fastqc: bool,
                 qc_extract: bool,
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
                 bt2_orient: str) -> list[str]:
    """ Run all steps of alignment for one FASTQ file or one pair of
    mated FASTQ files. """
    # Get the name of the set of references.
    refset = path.RefsetSeqInFilePath.parse(fasta).refset
    # List the reference names.
    refs = [ref for ref, _ in FastaParser(fasta).parse()]
    if not refs:
        logger.critical(f"No reference sequences for {fq_inp} in {fasta}")
        return list()
    # Determine whether alignment needs to be run.
    if not rerun:
        # If alignment is not required to be rerun, then check if all
        # the expected output files already exist.
        expected_outputs = infer_outputs(out_dir, fq_inp.sample, refs)
        if all(out_path.is_file() for out_path in expected_outputs):
            # If all the output files already exist, just return them.
            logger.warning(f"Skipping alignment for {fq_inp} because all "
                           f"expected output files already exist.")
            return list(map(str, expected_outputs))
    # FASTQC of the input
    if fastqc:
        fastqc_out = path.FastqcOutDirPath(top=str(out_dir),
                                           module=path.Module.ALIGN,
                                           sample=fq_inp.sample,
                                           fastqc=path.Fastqc.QC_INPUT).path
        run_fastqc(fq_inp, fastqc_out, qc_extract)
    # Trimming
    if cut:
        fq_cut = fq_inp.edit(top=str(temp_dir),
                             module=path.Module.ALIGN,
                             step=path.Step.ALIGN_TRIM)
        run_cutadapt(fq_inp, fq_cut,
                     n_procs=n_procs,
                     cut_q1=cut_q1,
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
        # FASTQC after trimming
        if fastqc:
            fastqc_out = path.FastqcOutDirPath(top=str(out_dir),
                                               module=path.Module.ALIGN,
                                               sample=fq_inp.sample,
                                               fastqc=path.Fastqc.QC_TRIM).path
            run_fastqc(fq_cut, fastqc_out, qc_extract)
    else:
        fq_cut = None
    # Alignment
    sam_aligned = path.RefsetAlignmentTempFilePath(top=str(temp_dir),
                                                   module=path.Module.ALIGN,
                                                   step=path.Step.ALIGN_ALIGN,
                                                   sample=fq_inp.sample,
                                                   refset=refset,
                                                   ext=path.SAM_EXT).path
    run_bowtie2(fq_inp if fq_cut is None else fq_cut,
                bowtie2_index,
                sam_aligned,
                n_procs=n_procs,
                bt2_local=bt2_local,
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
    if not save_temp and fq_cut is not None:
        for fq_file in fq_cut.paths:
            fq_file.path.unlink(missing_ok=True)
    # Deduplication
    sam_deduped = path.RefsetAlignmentTempFilePath(top=str(temp_dir),
                                                   module=path.Module.ALIGN,
                                                   step=path.Step.ALIGN_DEDUP,
                                                   sample=fq_inp.sample,
                                                   refset=refset,
                                                   ext=path.SAM_EXT).path
    dedup_sam(sam_aligned, sam_deduped)
    if not save_temp:
        sam_aligned.unlink(missing_ok=True)
    # Sorting by coordinate
    bam_sorted = path.RefsetAlignmentTempFilePath(top=str(temp_dir),
                                                  module=path.Module.ALIGN,
                                                  step=path.Step.ALIGN_SORT,
                                                  sample=fq_inp.sample,
                                                  refset=refset,
                                                  ext=path.BAM_EXT).path
    sort_xam(sam_deduped, bam_sorted, name=False)
    if not save_temp:
        sam_deduped.unlink(missing_ok=True)
    # Splitting into one BAM file for each reference
    index_bam_file(bam_sorted, n_procs=n_procs)  # index is required
    bams_out: list[str] = list()
    for ref in refs:
        bam_split = path.OneRefAlignmentOutFilePath(top=str(out_dir),
                                                    module=path.Module.ALIGN,
                                                    sample=fq_inp.sample,
                                                    ref=ref,
                                                    ext=path.BAM_EXT).path
        view_xam(bam_sorted, bam_split, ref=ref)
        bams_out.append(str(bam_split))
    if not save_temp:
        bam_sorted.unlink(missing_ok=True)
    # Return a list of the split BAM files.
    return bams_out


def run_steps_fqs(fq_units: list[FastqUnit],
                  fasta: pathlib.Path,
                  *,
                  max_procs: int,
                  parallel: bool,
                  out_dir: pathlib.Path,
                  temp_dir: pathlib.Path,
                  save_temp: bool,
                  **kwargs) -> tuple[str, ...]:
    """ Run all steps of alignment for one or more FASTQ files or pairs
    of mated FASTQ files. """
    # Remove duplicate FASTQ units.
    fq_units = deduplicate_fqs(fq_units)
    # Ensure that at least one FASTQ unit was given.
    if not fq_units:
        logger.critical("No FASTQ files were given for alignment.")
        return ()
    # Validate the maximum number of processes.
    if max_procs < 1:
        logger.warning("Maximum CPUs must be ≥ 1: setting to 1")
        max_procs = 1
    # Write the temporary FASTA file for each demultiplexed FASTQ.
    temp_refs = set(filter(None, (fq_unit.ref for fq_unit in fq_units)))
    temp_ref_paths = write_temp_ref_files(temp_dir, fasta, temp_refs)
    # Temporary index files for each fasta without a pre-built index.
    temp_bowtie2_indexes: dict[str, pathlib.Path] = dict()
    if all(index.is_file() for index in get_fasta_index_paths(fasta)):
        # Bowtie2 index for the main FASTA exists.
        fasta_bowtie2_index = fasta.with_suffix("")
    else:
        # Bowtie2 index does not already exist.
        fasta_bowtie2_index = None
    try:
        # One alignment task will be created for each FASTQ unit.
        # Get the arguments for each task, including n_procs_per_task.
        iter_args = list()
        for fq in fq_units:
            if fq.ref:
                # If the FASTQ came from demultiplexing (so contains
                # reads from only one reference), then align to the
                # FASTA of only that reference.
                try:
                    temp_fasta = temp_ref_paths[fq.ref]
                except KeyError:
                    # If the FASTA with that reference does not exist,
                    # then log an error and skip this FASTQ.
                    logger.error(
                        f"Skipped FASTQ {fq} because reference '{fq.ref}' "
                        f"was not found in FASTA file {fasta}")
                    continue
                # Build a Bowtie2 index for the single FASTA file.
                temp_bowtie2_indexes[fq.ref] = temp_fasta.with_suffix("")
                index_fasta_file(temp_fasta,
                                 temp_bowtie2_indexes[fq.ref],
                                 max_procs)
                # Add these arguments to the lists of arguments that
                # will be passed to run_steps_fq.
                iter_args.append((fq, temp_fasta, temp_bowtie2_indexes[fq.ref]))
            else:
                # If the FASTQ may contain reads from ≥ 1 references,
                # then align to the FASTA file with all references.
                if fasta_bowtie2_index is None:
                    # Build a Bowtie2 index for the main FASTA file.
                    refset = path.RefsetSeqInFilePath.parse(fasta).refset
                    fasta_bowtie2_index = path.RefsetBowtie2IndexTempPrefix(
                        top=str(temp_dir), module=path.Module.ALIGN.value,
                        step=path.Step.ALIGN_REFS.value, refset=refset).path
                    fasta_bowtie2_index.parent.mkdir(parents=True,
                                                     exist_ok=True)
                    index_fasta_file(fasta, fasta_bowtie2_index, max_procs)
                # Add these arguments to the lists of arguments that
                # will be passed to run_steps_fq.
                iter_args.append((fq, fasta, fasta_bowtie2_index))
        # Determine how to parallelize each alignment task.
        n_tasks_parallel, n_procs_per_task = get_num_parallel(len(fq_units),
                                                              max_procs,
                                                              parallel,
                                                              hybrid=True)
        # Pass the keyword arguments to every call of run_steps_fq.
        partial_run_steps_fq = partial(run_steps_fq,
                                       **{**dict(n_procs=n_procs_per_task,
                                                 out_dir=out_dir,
                                                 temp_dir=temp_dir,
                                                 save_temp=save_temp),
                                          **kwargs})
        if n_tasks_parallel > 1:
            # Process multiple FASTQ files simultaneously.
            with Pool(n_tasks_parallel) as pool:
                tasks = pool.starmap(partial_run_steps_fq, iter_args)
                bams = tuple(itertools.chain(*tasks))
        else:
            # Process the FASTQ files sequentially.
            tasks = itertools.starmap(partial_run_steps_fq, iter_args)
            bams = tuple(itertools.chain(*tasks))
    finally:
        if not save_temp:
            # Delete the temporary files before exiting.
            for ref_file in temp_ref_paths.values():
                logger.debug(
                    f"Deleting temporary reference file: {ref_file}")
                ref_file.unlink(missing_ok=True)
            # Delete the temporary index files before exiting.
            for ref, prefix in temp_bowtie2_indexes.items():
                for index_file in get_fasta_index_paths(prefix):
                    logger.debug(
                        f"Deleting temporary index file: {index_file}")
                    index_file.unlink(missing_ok=True)
    # Return a tuple of the final alignment map files.
    return bams
