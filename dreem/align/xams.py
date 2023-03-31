from logging import getLogger
from pathlib import Path
import re
from typing import BinaryIO

from ..util import path
from ..util.logs import log_process
from ..util.shell import run_cmd, SAMTOOLS_CMD


logger = getLogger(__name__)


# SAM file format specifications
SAM_HEADER = b"@"
SAM_DELIMITER = b"\t"
SAM_ALIGN_SCORE = b"AS:i:"
SAM_EXTRA_SCORE = b"XS:i:"


def index_bam(bam: Path, n_padd: int = 0):
    """ Build an index of a BAM file using ```samtools index```. """
    logger.info(f"Building BAM index of {bam}: {bam.with_suffix(path.BAI_EXT)}")
    cmd = [SAMTOOLS_CMD, "index", "-@", n_padd - 1, bam]
    process = run_cmd(cmd, capture_output=True)
    log_process(logger, process)


def sort_xam(xam_inp: Path, xam_out: Path, *,
             name: bool = False, n_padd: int = 0):
    """ Sort a SAM or BAM file using ```samtools sort```. """
    logger.info(f"Sorting {xam_inp} to {xam_out}")
    cmd = [SAMTOOLS_CMD, "sort", "-@", n_padd]
    if name:
        # Sort by name instead of coordinate.
        cmd.append("-n")
    cmd.extend(["-o", xam_out, xam_inp])
    # Make the output directory.
    xam_out.parent.mkdir(parents=True, exist_ok=True)
    process = run_cmd(cmd, capture_output=True)
    log_process(logger, process)


def view_xam(xam_inp: Path,
             xam_out: Path, *,
             ref: str | None = None,
             end5: int | None = None,
             end3: int | None = None,
             n_padd: int = 0):
    """ Convert between SAM and BAM formats, or extract reads aligning
    to a specific reference/section using ```samtools view```. """
    logger.info(f"Viewing {xam_inp} as {xam_out}")
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
    # Run the command.
    process = run_cmd(cmd, capture_output=True)
    log_process(logger, process)


def dedup_sam(sam_inp: Path, sam_out: Path):
    """ Remove SAM reads that map equally to multiple locations. """

    logger.info(f"Deduplicating {sam_inp} to {sam_out}")

    pattern_a = re.compile(SAM_ALIGN_SCORE + rb"(\d+)")
    pattern_x = re.compile(SAM_EXTRA_SCORE + rb"(\d+)")

    min_fields = 11
    max_flag = 4095  # 2^12 - 1

    def get_score(line: bytes, ptn: re.Pattern[bytes]):
        return (float(match.groups()[0])
                if (match := ptn.search(line)) else None)

    def is_best_alignment(line: bytes):
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

    def iter_paired(sam: BinaryIO, line1: bytes):
        """ For each pair of reads, yield the pair of alignments for
        which both the forward alignment and the reverse alignment in
        the pair scored best among all alignments for the forward and
        reverse reads, respectively. Exclude pairs for which the forward
        and/or reverse read aligned equally well to multiple locations,
        or for which the best alignments for the forward and reverse
        reads individually are not part of the same alignment pair. """
        while line1:
            if not (line2 := sam.readline()):
                raise ValueError(f"SAM file {sam_inp} is paired-end but has a "
                                 f"mate 1 with no mate 2: '{line1.decode()}'")
            if ((name1 := line1.split(SAM_DELIMITER, 1)[0])
                    != (name2 := line2.split(SAM_DELIMITER, 1)[0])):
                raise ValueError(f"SAM file {sam_inp} has a mate 1 ({name1}) "
                                 f"and mate 2 ({name2}) with different names")
            if is_best_alignment(line1) and is_best_alignment(line2):
                yield line1
                yield line2
            line1 = sam.readline()

    def iter_single(sam: BinaryIO, line: bytes):
        """ For each read, yield the best-scoring alignment, excluding
        reads that aligned equally well to multiple locations. """
        while line:
            if is_best_alignment(line):
                yield line
            line = sam.readline()

    # Make the output directory.
    sam_out.parent.mkdir(parents=True, exist_ok=True)

    # Deduplicate the alignments.
    with (open(sam_inp, "rb") as sami, open(sam_out, "xb") as samo):
        # Copy the entire header from the input to the output SAM file.
        while (line_ := sami.readline()).startswith(SAM_HEADER):
            samo.write(line_)
        # Copy only the reads that mapped best to one location from the
        # input to the output SAM file.
        if line_:
            # Determine whether the reads in the SAM file are paired.
            if read_is_paired(line_):
                lines = iter_paired(sami, line_)
            else:
                lines = iter_single(sami, line_)
            # Write the selected lines to the output SAM file.
            for line_ in lines:
                samo.write(line_)
        else:
            logger.warning(f"SAM file {sam_inp} contains no reads")


'''
primer1 = "CAGCACTCAGAGCTAATACGACTCACTATA"
primer1rc = "TATAGTGAGTCGTATTAGCTCTGAGTGCTG"
primer2 = "TGAAGAGCTGGAACGCTTCACTGA"
primer2rc = "TCAGTGAAGCGTTCCAGCTCTTCA"
adapters5 = (primer1, primer2rc)
adapters3 = (primer2, primer1rc)
'''
