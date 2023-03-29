from functools import partial
import itertools
import logging
from multiprocessing import Pool
import pathlib
from typing import Iterable

from ..align.reads import (FastqUnit, run_fastqc, run_cutadapt, run_bowtie2,
                           dedup_sam, sort_xam, view_xam, index_bam_file,
                           index_fasta_file, get_fasta_index_paths)
from ..util import path
from ..util.seq import FastaParser, FastaWriter
from ..util.util import get_num_parallel

logger = logging.getLogger(__name__)


def write_temp_ref_files(temp_dir: pathlib.Path,
                         refset_path: pathlib.Path,
                         refs: set[str],
                         n_procs: int):
    """ Write temporary FASTA files, each containing one reference that
    corresponds to a FASTQ file from demultiplexing. """
    ref_paths: dict[str, tuple[pathlib.Path, pathlib.Path]] = dict()
    if refs:
        # Parse the FASTA only if there are any references to write.
        for ref, seq in FastaParser(refset_path).parse():
            if ref in refs:
                # Write the reference sequence to a temporary FASTA file
                # only if at least one demultiplexed FASTQ file uses it.
                ref_path = path.OneRefSeqTempFilePath(top=str(temp_dir),
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
                # Build a Bowtie2 index of the temporary FASTA file.
                index_prefix = ref_path.with_suffix("")
                index_fasta_file(ref_path, index_prefix, n_procs)
                # Record the temporary FASTA and index prefix.
                ref_paths[ref] = ref_path, index_prefix
    if missing := sorted(refs - set(ref_paths.keys())):
        logger.warning(f"Missing references in {refset_path}: "
                       + ", ".join(missing))
    return ref_paths


def fq_pipeline(fq_inp: FastqUnit,
                fasta: pathlib.Path,
                bowtie2_index: pathlib.Path,
                *,
                n_procs: int,
                out_dir: str,
                temp_dir: str,
                save_temp: bool,
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
                bt2_orient: str) -> list[pathlib.Path]:
    """ Run all steps of the alignment pipeline for one FASTQ file or
    one pair of mated FASTQ files. """
    # Get the name of the set of references.
    refset = path.RefsetSeqInFilePath.parse(fasta).refset
    # List the reference names.
    refs = [ref for ref, _ in FastaParser(fasta).parse()]
    if not refs:
        logger.critical(f"No reference sequences for {fq_inp} in {fasta}")
        return list()
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
    bams_out: list[pathlib.Path] = list()
    for ref in refs:
        bam_split = path.OneRefAlignmentOutFilePath(top=str(out_dir),
                                                    module=path.Module.ALIGN,
                                                    sample=fq_inp.sample,
                                                    ref=ref,
                                                    ext=path.BAM_EXT).path
        view_xam(bam_sorted, bam_split, ref=ref)
        bams_out.append(bam_split)
    if not save_temp:
        bam_sorted.unlink(missing_ok=True)
    # Return a list of the split BAM files.
    return bams_out


def fqs_pipeline(fq_units: list[FastqUnit],
                 fasta: pathlib.Path, *,
                 max_procs: int,
                 parallel: bool,
                 out_dir: pathlib.Path,
                 temp_dir: pathlib.Path,
                 save_temp: bool,
                 **kwargs) -> list[pathlib.Path]:
    """ Run all steps of alignment for one or more FASTQ files or pairs
    of mated FASTQ files. """
    # Validate the maximum number of processes.
    if max_procs < 1:
        logger.warning("Maximum CPUs must be ≥ 1: setting to 1")
        max_procs = 1
    # Get the name of the reference for every demultiplexed FASTQ.
    temp_refs = set(filter(None, (fq_unit.ref for fq_unit in fq_units)))
    # Write a temporary FASTA file and Bowtie2 index for each
    # demultiplexed FASTQ.
    temp_ref_paths = write_temp_ref_files(temp_dir, fasta, temp_refs, max_procs)
    # Check if the main FASTA file already has a Bowtie2 index.
    if all(index.is_file() for index in get_fasta_index_paths(fasta)):
        # Bowtie2 index for the main FASTA exists.
        main_bowtie2_index = fasta.with_suffix("")
    else:
        # Bowtie2 index does not already exist.
        main_bowtie2_index = None
    try:
        # Get the arguments for each alignment task.
        iter_args: list[tuple[FastqUnit, pathlib.Path, pathlib.Path]] = list()
        # One alignment task will be created for each FASTQ unit.
        for fq in fq_units:
            if fq.ref:
                # If the FASTQ came from demultiplexing (so contains
                # reads from only one reference), then align to the
                # FASTA of only that reference.
                try:
                    temp_fasta, temp_index = temp_ref_paths[fq.ref]
                except KeyError:
                    # If the FASTA with that reference does not exist,
                    # then log an error and skip this FASTQ.
                    logger.critical(
                        f"Skipped FASTQ {fq} because reference '{fq.ref}' "
                        f"was not found in FASTA file {fasta}")
                    continue

                # Add these arguments to the lists of arguments that
                # will be passed to fq_pipeline.
                iter_args.append((fq, temp_fasta, temp_index))
            else:
                # If the FASTQ may contain reads from ≥ 1 references,
                # then align to the FASTA file with all references.
                if main_bowtie2_index is None:
                    # Determine the name of the set of references.
                    refset = path.RefsetSeqInFilePath.parse(fasta).refset
                    # Determine the path of the temporary Bowtie2 index
                    # of the main FASTA file.
                    main_bowtie2_index = path.RefsetBowtie2IndexTempPrefix(
                        top=str(temp_dir), module=path.Module.ALIGN.value,
                        step=path.Step.ALIGN_REFS.value, refset=refset).path
                    # Make its parent directory if it does not exist.
                    main_bowtie2_index.parent.mkdir(parents=True, exist_ok=True)
                    # Build the Bowtie2 index.
                    index_fasta_file(fasta, main_bowtie2_index, max_procs)
                # Add these arguments to the lists of arguments that
                # will be passed to fq_pipeline.
                iter_args.append((fq, fasta, main_bowtie2_index))
        # Determine how to parallelize each alignment task.
        n_tasks_parallel, n_procs_per_task = get_num_parallel(len(fq_units),
                                                              max_procs,
                                                              parallel,
                                                              hybrid=True)
        # Pass the keyword arguments to every call of fq_pipeline.
        partial_fq_pipeline = partial(fq_pipeline,
                                      **{**dict(n_procs=n_procs_per_task,
                                                out_dir=out_dir,
                                                temp_dir=temp_dir,
                                                save_temp=save_temp),
                                         **kwargs})
        if n_tasks_parallel > 1:
            # Process the FASTQ files simultaneously.
            with Pool(n_tasks_parallel) as pool:
                tasks = pool.starmap(partial_fq_pipeline, iter_args)
                bams = list(itertools.chain(*tasks))
        else:
            # Process the FASTQ files sequentially.
            tasks = itertools.starmap(partial_fq_pipeline, iter_args)
            bams = list(itertools.chain(*tasks))
    finally:
        if not save_temp:
            # Delete the temporary files before exiting.
            for ref_file, index_prefix in temp_ref_paths.values():
                # Reference file
                logger.debug(
                    f"Deleting temporary reference file: {ref_file}")
                ref_file.unlink(missing_ok=True)
                # Index files
                for index_file in get_fasta_index_paths(index_prefix):
                    logger.debug(
                        f"Deleting temporary index file: {index_file}")
                    index_file.unlink(missing_ok=True)
    # Return the final alignment map files.
    return bams


def figure_alignments(fq_units: list[FastqUnit], refs: set[str]):
    """ Return a ```dict``` of every expected alignment of a sample to a
    reference sequence. Check for and remove duplicates. """
    # Map each combination of a sample and reference to a FASTQ unit.
    alignments: dict[tuple[str, str], FastqUnit] = dict()
    # Keep track of any duplicate sample-reference pairs.
    duplicates: set[tuple[str, str]] = set()
    for fq_unit in fq_units:
        # Determine which references the FASTQ reads could come from.
        if fq_unit.ref is None:
            # The FASTQ contains reads from potentially all references.
            fq_refs = refs
        else:
            # The FASTQ contains reads from only one reference.
            # Confirm that the reference actually exists.
            if fq_unit.ref not in refs:
                logger.critical(
                    f"Reference '{fq_unit.ref}' of {fq_unit} not found")
                continue
            fq_refs = fq_unit.ref,
        # Add each sample-reference pair to the expected alignments.
        for ref in fq_refs:
            sample_ref = fq_unit.sample, ref
            if sample_ref in duplicates:
                # Skip the sample-reference pair if it is a duplicate.
                continue
            try:
                # Test if the sample-reference pair is already in the
                # dict of alignments. If so, then remove it.
                alignments.pop(sample_ref)
            except KeyError:
                # If not, then add the FASTQ to the dict of alignments,
                # keyed by its sample-reference pair.
                alignments[sample_ref] = fq_unit
            else:
                # If so, then flag it as a duplicate.
                logger.warning(f"Duplicate sample and reference: {sample_ref}")
                duplicates.add(sample_ref)
    # Return a duplicate-free dict of alignments.
    return alignments


def check_fqs_bams(alignments: dict[tuple[str, str], FastqUnit],
                   out_dir: pathlib.Path):
    """ Return every FASTQ unit on which alignment must be run and every
    expected BAM file that already exists. """
    alignments_missing: dict[tuple[str, str], FastqUnit] = dict()
    bams_existing: list[pathlib.Path] = list()
    for (sample, ref), fq_unit in alignments.items():
        # Determine the path of the BAM file expected to result from the
        # alignment of the sample to the reference.
        bam_expect = path.OneRefAlignmentOutFilePath(top=str(out_dir),
                                                     module=path.Module.ALIGN,
                                                     sample=sample,
                                                     ref=ref,
                                                     ext=path.BAM_EXT).path
        if bam_expect.is_file():
            # If the BAM file already exists, then add it to the dict of
            # BAM files that have already been aligned.
            bams_existing.append(bam_expect)
        else:
            # If at least one BAM file for a FASTQ unit does not exist,
            # then align the FASTQ.
            alignments_missing[sample, ref] = fq_unit
    return alignments_missing, bams_existing


def merge_nondemult_fqs(fq_units: Iterable[FastqUnit]):
    """ For every FASTQ that is not demultiplexed, merge all the keys
    that map to the FASTQ into one key: (sample, None). Merging ensures
    that every non-demultiplexed FASTQ is aligned only once to the whole
    set of references, not once for every reference in the set. This
    function is essentially the inverse of ```figure_alignments```. """
    merged: dict[tuple[str, str | None], FastqUnit] = dict()
    for fq_unit in fq_units:
        merged[fq_unit.sample, fq_unit.ref] = fq_unit
    return list(merged.values())


def list_fqs_bams(fq_units: list[FastqUnit],
                  refs: set[str],
                  out_dir: pathlib.Path):
    """ List every FASTQ that needs to be aligned and every expected BAM
    file that already exists. """
    # Determine all possible alignments of a sample and reference.
    alignments = figure_alignments(fq_units, refs)
    # Determine which alignments need to be / have already been run.
    alignments_missing, bams_existing = check_fqs_bams(alignments, out_dir)
    # Merge entries for each non-demultiplexed FASTQ.
    fqs_to_align = merge_nondemult_fqs(alignments_missing.values())
    return fqs_to_align, set(bams_existing)


def get_bam_files(fq_units: list[FastqUnit],
                  fasta: pathlib.Path,
                  *,
                  out_dir: pathlib.Path,
                  rerun: bool,
                  **kwargs) -> tuple[str, ...]:
    """ Run the alignment pipeline and return a tuple of all BAM files
    from the pipeline. """
    # Get the names of all reference sequences.
    refs = {ref for ref, _ in FastaParser(fasta).parse()}
    if not refs:
        raise ValueError(f"Got no reference sequences in {fasta}")
    if rerun:
        # Rerun all alignments.
        fqs_to_align = fq_units
        bams: set[pathlib.Path] = set()
    else:
        # Run only the alignments whose outputs do not yet exist.
        fqs_to_align, bams = list_fqs_bams(fq_units, refs, out_dir)
    if fqs_to_align:
        # Align all FASTQs that need to be aligned.
        bams_new = set(fqs_pipeline(fq_units, fasta, out_dir=out_dir, **kwargs))
    else:
        logger.warning("All given FASTQ files have already been aligned")
        bams_new = set()
    # Merge the existing and new BAM paths into a tuple of strings.
    return tuple(map(str, bams | bams_new))
