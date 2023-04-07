from logging import getLogger
from pathlib import Path
import re
from typing import BinaryIO

from ..util import path
from ..util.shell import run_cmd, SAMTOOLS_CMD

logger = getLogger(__name__)

# SAM file format specifications
SAM_HEADER = b"@"
SAM_DELIMITER = b"\t"
SAM_ALIGN_SCORE = b"AS:i:"
SAM_EXTRA_SCORE = b"XS:i:"


def index_bam(bam: Path, n_padd: int = 0):
    """ Build an index of a BAM file using ```samtools index```. """
    logger.info(f"Began building BAM index of {bam}")
    cmd = [SAMTOOLS_CMD, "index", "-@", n_padd - 1, bam]
    run_cmd(cmd)
    logger.info(f"Ended building BAM index of {bam}: "
                f"{bam.with_suffix(path.BAI_EXT)}")


def sort_xam(xam_inp: Path, xam_out: Path, *,
             name: bool = False, n_padd: int = 0):
    """ Sort a SAM or BAM file using ```samtools sort```. """
    logger.info(f"Began sorting {xam_inp}")
    cmd = [SAMTOOLS_CMD, "sort", "-@", n_padd]
    if name:
        # Sort by name instead of coordinate.
        cmd.append("-n")
    cmd.extend(["-o", xam_out, xam_inp])
    # Make the output directory.
    xam_out.parent.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Ensured directory: {xam_out.parent}")
    run_cmd(cmd)
    logger.info(f"Ended sorting {xam_inp} to {xam_out}")


def view_xam(xam_inp: Path,
             xam_out: Path, *,
             ref: str | None = None,
             end5: int | None = None,
             end3: int | None = None,
             n_padd: int = 0):
    """ Convert between SAM and BAM formats, or extract reads aligning
    to a specific reference/section using ```samtools view```. """
    logger.info(f"Began viewing {xam_inp}")
    cmd = [SAMTOOLS_CMD, "view", "-@", n_padd, "-h"]
    if xam_out.suffix == path.BAM_EXT:
        # Write a binary (BAM) file.
        cmd.append("-b")
    # Input and output files
    cmd.extend(("-o", xam_out, xam_inp))
    # Reference and section specification
    if ref is not None:
        if end5 is not None and end3 is not None:
            # View only reads aligning to a section of this reference.
            cmd.append(f"{ref}:{end5}-{end3}")
        else:
            # View only reads aligning to this reference.
            cmd.append(ref)
            if end5 is not None:
                logger.warning(f"Got end5 = {end5} but not end3")
            if end3 is not None:
                logger.warning(f"Got end3 = {end3} but not end5")
    elif end5 is not None or end5 is not None:
        logger.warning(f"Options end5 and end3 require a reference name")
    # Make the output directory.
    xam_out.parent.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Ensured directory: {xam_out.parent}")
    # Run the command.
    run_cmd(cmd)
    logger.info(f"Ended viewing {xam_inp} as {xam_out}")


def dedup_sam(sam_inp: Path, sam_out: Path):
    """ Remove SAM reads that map equally to multiple locations. """

    logger.info(f"Began deduplicating {sam_inp}")

    pattern_a = re.compile(SAM_ALIGN_SCORE + rb"(\d+)")
    pattern_x = re.compile(SAM_EXTRA_SCORE + rb"(\d+)")

    min_fields = 11
    max_flag = 4095  # 2^12 - 1

    def get_score(line: bytes, ptn: re.Pattern[bytes]):
        return (float(match.groups()[0])
                if (match := ptn.search(line)) else None)

    def is_best_align(line: bytes):
        return ((score_x := get_score(line, pattern_x)) is None
                or score_x < get_score(line, pattern_a))

    def read_is_paired(line: bytes):
        info = line.split()
        if len(info) < min_fields:
            raise ValueError(f"Invalid SAM line:\n{line.decode()}")
        flag = int(info[1])
        if flag < 0 or flag > max_flag:
            raise ValueError(f"Invalid SAM flag: {flag}")
        return bool(flag % 2)

    def write_summary_single(written: int, skipped: int):
        total = written + skipped
        f_written = (round(100 * written / total, 3) if total
                     else float('nan'))
        return (f"Summary of removing multiply-mapped reads from {sam_inp}\n"
                f"Total number of single-end reads: {total}\n"
                f"Uniquely-mapped reads (written):  {written}\n"
                f"Multiply-mapped reads (skipped):  {skipped}\n"
                f"Percent uniquely mapped (written): {f_written} %")

    def write_summary_paired(pairs_written: int, pairs_skipped: int,
                             mates_written: int, mates_skipped: int):
        mates = mates_written + mates_skipped
        pairs = pairs_written + pairs_skipped
        pairs_one_mate = 2 * pairs - mates
        pairs_both_mates = pairs - pairs_one_mate
        f_written = (round(100 * pairs_written / pairs, 3) if pairs
                     else float('nan'))
        return (f"Summary of removing multiply-mapped reads from {sam_inp}\n"
                f"Total number of mates aligning: {mates}\n"
                f"Total number of pairs of mates: {pairs}\n"
                f"Pairs with both mates aligning: {pairs_both_mates}"
                f"Pairs with a mate not aligning: {pairs_one_mate}"
                f"Uniquely-mapped pairs (written):  {pairs_written}\n"
                f"Multiply-mapped pairs (skipped):  {pairs_skipped}\n"
                f"Percent uniquely mapped (written): {f_written} %")

    def iter_single(sam: BinaryIO, line: bytes):
        """ For each read, yield the best-scoring alignment, excluding
        reads that aligned equally well to multiple locations. """
        n_copy = 0
        n_skip = 0
        while line:
            if is_best_align(line):
                n_copy += 1
                yield line
            else:
                n_skip += 1
            line = sam.readline()
        logger.debug(write_summary_single(n_copy, n_skip))

    def iter_paired(sam: BinaryIO, line1: bytes):
        """ For each pair of reads, yield the pair of alignments for
        which both the forward alignment and the reverse alignment in
        the pair scored best among all alignments for the forward and
        reverse reads, respectively. Exclude pairs for which the forward
        and/or reverse read aligned equally well to multiple locations,
        or for which the best alignments for the forward and reverse
        reads individually are not part of the same alignment pair. """
        n_pairs_written = 0
        n_pairs_skipped = 0
        n_mates_written = 0
        n_mates_skipped = 0
        while line1:
            if SAM_DELIMITER not in line1:
                # If the line is blank, skip it.
                line1 = sam.readline()
                continue
            # Check if there is at least one more line in the file.
            if (line2 := sam.readline()) and SAM_DELIMITER in line2:
                # Check if the reads on lines 1 and 2 are mates.
                if (line1.split(SAM_DELIMITER, 1)[0]
                        == line2.split(SAM_DELIMITER, 1)[0]):
                    # If they are mates of each other, check if the
                    # lines represent the best alignment for both.
                    if is_best_align(line1) and is_best_align(line2):
                        # If so, then yield both mates.
                        n_pairs_written += 1
                        n_mates_written += 2
                        yield line1
                        yield line2
                    else:
                        # Otherwise, skip both mates.
                        n_pairs_skipped += 1
                        n_mates_skipped += 2
                    # Read the next line from the file.
                    line1 = sam.readline()
                    continue
            # If line2 does not have the mate of line1, then check only
            # line1 now and consider line2 on the next iteration.
            if is_best_align(line1):
                # Yield only line1.
                n_pairs_written += 1
                n_mates_written += 1
                yield line1
            else:
                # Skip line1.
                n_pairs_skipped += 1
                n_mates_skipped += 1
            # Since line2 was not used yet, consider it on the next
            # iteration instead of reading another line from the file.
            line1 = line2
        logger.debug(write_summary_paired(n_pairs_written, n_pairs_skipped,
                                          n_mates_written, n_mates_skipped))

    # Make the output directory.
    sam_out.parent.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Ensuring directory: {sam_out.parent}")

    # Deduplicate the alignments.
    with (open(sam_inp, "rb") as sami, open(sam_out, "xb") as samo):
        # Copy the entire header from the input to the output SAM file.
        n_header = 0
        while (line_ := sami.readline()).startswith(SAM_HEADER):
            samo.write(line_)
            n_header += 1
        logger.debug(
            f"Copied {n_header} header lines from {sam_inp} to {sam_out}")
        # Copy only the reads that mapped best to one location from the
        # input to the output SAM file.
        if line_:
            is_paired = read_is_paired(line_)
            # Determine whether the reads in the SAM file are paired.
            if is_paired:
                endedness = "paired"
                lines = iter_paired(sami, line_)
            else:
                endedness = "single"
                lines = iter_single(sami, line_)
            logger.debug(f"Deduplicating {endedness}-end SAM file: {sam_inp}")
            # Write the selected lines to the output SAM file.
            for line_ in lines:
                samo.write(line_)
        else:
            logger.critical(f"SAM file {sam_inp} contained no reads")

    logger.info(f"Ended deduplicating {sam_inp} to {sam_out}")


'''
primer1 = "CAGCACTCAGAGCTAATACGACTCACTATA"
primer1rc = "TATAGTGAGTCGTATTAGCTCTGAGTGCTG"
primer2 = "TGAAGAGCTGGAACGCTTCACTGA"
primer2rc = "TCAGTGAAGCGTTCCAGCTCTTCA"
adapters5 = (primer1, primer2rc)
adapters3 = (primer2, primer1rc)
'''
