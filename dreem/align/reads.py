import itertools
import logging
import pathlib
import re
from functools import cached_property
from typing import BinaryIO, Iterable

from ..util import cli, path
from ..util.logs import log_process
from ..util.shell import (run_cmd, BOWTIE2_CMD, BOWTIE2_BUILD_CMD,
                          CUTADAPT_CMD, FASTQC_CMD, SAMTOOLS_CMD)

logger = logging.getLogger(__name__)

# SAM file format specifications
SAM_HEADER = b"@"
SAM_DELIMITER = b"\t"
SAM_ALIGN_SCORE = b"AS:i:"
SAM_EXTRA_SCORE = b"XS:i:"

# Bowtie2 parameters
MATCH_BONUS = "1"
MISMATCH_PENALTY = "1,1"
N_PENALTY = "0"
REF_GAP_PENALTY = "0,1"
READ_GAP_PENALTY = "0,1"
METRICS_INTERVAL = 60  # Write metrics once every 60 seconds.


def get_fasta_index_paths(fasta: pathlib.Path):
    """ Return the Bowtie 2 index paths for a FASTA file. """
    return [fasta.with_suffix(ext) for ext in path.BOWTIE2_INDEX_EXTS]


def index_fasta_file(fasta: pathlib.Path,
                     prefix: pathlib.Path,
                     n_procs: int = 1):
    """ Build a Bowtie2 index of a FASTA file. """
    logger.info(f"Building Bowtie2 index of FASTA {fasta}: {prefix}")
    # Generate and run the command.
    cmd = [BOWTIE2_BUILD_CMD, "--threads", n_procs, fasta, prefix]
    process = run_cmd(cmd, capture_output=True)
    log_process(logger, process)


def index_bam_file(bam: pathlib.Path, n_procs: int = 1):
    """ Build an index of a BAM file using ```samtools index```. """
    logger.info(f"Building BAM index of {bam}: {bam.with_suffix(path.BAI_EXT)}")
    cmd = [SAMTOOLS_CMD, "index", "-@", n_procs - 1, bam]
    process = run_cmd(cmd, capture_output=True)
    log_process(logger, process)


def sort_xam(xam_inp: pathlib.Path, xam_out: pathlib.Path, name: bool):
    """ Sort a SAM or BAM file using ```samtools sort```. """
    logger.info(f"Sorting {xam_inp} to {xam_out}")
    cmd = [SAMTOOLS_CMD, "sort"]
    if name:
        cmd.append("-n")
    cmd.extend(["-o", xam_out, xam_inp])
    # Make the output directory.
    xam_out.parent.mkdir(parents=True, exist_ok=True)
    process = run_cmd(cmd, capture_output=True)
    log_process(logger, process)


def view_xam(xam_inp: pathlib.Path,
             xam_out: pathlib.Path,
             ref: str | None = None,
             end5: int | None = None,
             end3: int | None = None):
    """ Convert between SAM and BAM formats, or extract reads aligning
    to a specific reference/section using ```samtools view```. """
    logger.info(f"Viewing {xam_inp} as {xam_out}")
    cmd = [SAMTOOLS_CMD, "view", "-h"]
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


class FastqUnit(object):
    """
    Unified interface for handling the following sets of sequencing reads:

    - One FASTQ file of single-end reads from one sample
    - One FASTQ file of interleaved, paired-end reads from one sample
    - Two FASTQ files of mate 1 and mate 2 paired-end reads from one sample
    - One FASTQ file of single-end reads originating from one reference sequence
      in one sample
    - One FASTQ file of interleaved, paired-end reads originating from one
      reference sequence in one sample
    - Two FASTQ files of mate 1 and mate 2 paired-end reads originating from one
      reference sequence in one sample

    """

    MAX_PHRED_ENC = 127  # 2^7 - 1

    KEYF_SINGLE = "fastqs"
    KEYF_INTER = "fastqi"
    KEYF_MATE1 = "fastq1"
    KEYF_MATE2 = "fastq2"

    KEYD_SINGLE = "fastqs_dir"
    KEYD_INTER = "fastqi_dir"
    KEYD_MATE12 = "fastq12_dir"

    BOWTIE2_FLAGS = {KEYF_SINGLE: "-U",
                     KEYF_INTER: "--interleaved",
                     KEYF_MATE1: "-1",
                     KEYF_MATE2: "-2"}

    def __init__(self, /, *,
                 fastqs: pathlib.Path | None = None,
                 fastqi: pathlib.Path | None = None,
                 fastq1: pathlib.Path | None = None,
                 fastq2: pathlib.Path | None = None,
                 phred_enc: int,
                 one_ref: bool):
        if fastqs:
            if fastqi or fastq1 or fastq2:
                raise TypeError("Got too many FASTQ files")
            self.inputs: dict[str, pathlib.Path] = {self.KEYF_SINGLE: fastqs}
            self.paired = False
            self.interleaved = False
        elif fastqi:
            if fastq1 or fastq2:
                raise TypeError("Got too many FASTQ files")
            self.inputs: dict[str, pathlib.Path] = {self.KEYF_INTER: fastqi}
            self.paired = True
            self.interleaved = True
        elif fastq1:
            if not fastq2:
                raise TypeError("Got fastq1 but not fastq2")
            self.inputs: dict[str, pathlib.Path] = {self.KEYF_MATE1: fastq1,
                                                    self.KEYF_MATE2: fastq2}
            self.paired = True
            self.interleaved = False
        elif fastq2:
            raise TypeError("Got fastq2 but not fastq1")
        if phred_enc < 0 or phred_enc > self.MAX_PHRED_ENC:
            raise ValueError(f"Invalid Phred encoding: {phred_enc}")
        self.phred_enc = phred_enc
        self.one_ref = one_ref
        logger.debug(f"Instantiated a {self.__class__.__name__} with "
                     + ", ".join(f"{k} = {v} ({type(v).__name__})"
                                 for k, v in self.inputs.items())
                     + f", phred_enc = {phred_enc}, one_ref = {one_ref}")

    @property
    def phred_arg(self):
        return f"--phred{self.phred_enc}"

    @property
    def kind(self):
        if self.paired:
            if self.interleaved:
                return "paired-end interleaved FASTQ file"
            return "paired-end separated FASTQ files"
        return "single-end FASTQ file"

    @cached_property
    def parent(self) -> pathlib.Path:
        """ Return the parent directory of the FASTQ file(s). """
        parents = [inp.parent for inp in self.inputs.values()]
        if not parents:
            raise TypeError("Not parent directory")
        if any(parent != parents[0] for parent in parents[1:]):
            raise ValueError("More than one parent directory")
        return parents[0]

    @cached_property
    def paths(self) -> list[(path.OneRefReadsInFilePath |
                             path.OneRefReads1InFilePath |
                             path.OneRefReads2InFilePath |
                             path.SampleReadsInFilePath |
                             path.SampleReads1InFilePath |
                             path.SampleReads2InFilePath)]:
        """ Return a path instance for each input FASTQ file. """
        if self.one_ref:
            dtypes = {self.KEYF_SINGLE: path.OneRefReadsInFilePath,
                      self.KEYF_INTER: path.OneRefReadsInFilePath,
                      self.KEYF_MATE1: path.OneRefReads1InFilePath,
                      self.KEYF_MATE2: path.OneRefReads2InFilePath}
        else:
            dtypes = {self.KEYF_SINGLE: path.SampleReadsInFilePath,
                      self.KEYF_INTER: path.SampleReadsInFilePath,
                      self.KEYF_MATE1: path.SampleReads1InFilePath,
                      self.KEYF_MATE2: path.SampleReads2InFilePath}
        return [dtypes[key].parse(fq_path)
                for key, fq_path in self.inputs.items()]

    @cached_property
    def sample(self):
        """ Return the name of the sample of the FASTQ file(s). """
        samples: list[str] = list({fq.sample for fq in self.paths})
        if len(samples) == 0:
            raise ValueError("No sample names")
        if len(samples) > 1:
            raise ValueError(f"Sample names of input files {self} disagree: "
                             + " ≠ ".join(samples))
        return samples[0]

    @cached_property
    def ref(self):
        """ Return the name of the reference of the FASTQ file(s). """
        refs: list[str] = list()
        # Get the reference name from each FASTQ that names a reference.
        for fq in self.paths:
            try:
                refs.append(fq.ref)
            except AttributeError:
                pass
        if len(refs) == 0:
            # No FASTQ files named a reference.
            return None
        if len(refs) > 1:
            raise ValueError(f"Reference names of input files {self} disagree: "
                             + " ≠ ".join(refs))
        return refs[0]

    def is_compatible_fasta(self, fasta: pathlib.Path, one_ref_fasta: bool):
        """
        Return whether a given FASTA file is compatible with the FASTQ,
        which means one of the following is true:
        - The FASTA has a set of references (one_ref_fasta=False) and
          the FASTQ has reads from an entire sample.
        - The FASTA has one reference (one_ref_fasta=True) and the FASTQ
          has reads from only the same reference.
        """
        if one_ref_fasta:
            return self.ref == path.OneRefSeqInFilePath.parse(fasta).ref
        return not self.one_ref

    @property
    def cutadapt_input_args(self):
        """ Return input file arguments for Cutadapt. """
        return tuple(self.inputs.values())

    @property
    def bowtie2_inputs(self):
        """ Return input file arguments for Bowtie2. """
        return tuple(itertools.chain(*[(self.BOWTIE2_FLAGS[key], fq)
                                       for key, fq in self.inputs.items()]))

    def edit(self, **fields):
        """
        Return a new FastqUnit by modifying all the FASTQ paths, using
        the given fields.

        Parameters
        ----------
        **fields: Any
            Keyword arguments to determine the new paths

        Returns
        -------
        FastqUnit
            New FastqUnit with the same keys and modified paths
        """
        return FastqUnit(**{k: fq.edit(**{**fields, path.EXT_KEY: fq.ext}).path
                            for k, fq in zip(self.inputs, self.paths,
                                             strict=True)},
                         phred_enc=self.phred_enc,
                         one_ref=self.one_ref)

    @classmethod
    def _from_files(cls, /, *,
                    key: str,
                    fq_files: Iterable[pathlib.Path],
                    one_ref: bool,
                    phred_enc: int):
        """ Yield a FastqUnit for each given path to a FASTQ file of
        single-end or interleaved paired-end reads. """
        for fq in fq_files:
            if fq.suffix not in path.FQ_EXTS:
                logger.debug(f"Skipping non-FASTQ file: '{fq}'")
                continue
            try:
                yield FastqUnit(**{key: fq}, phred_enc=phred_enc, one_ref=one_ref)
            except path.PathError as error:
                logger.critical(
                    f"Failed to parse path '{fq}' as a FASTQ: {error}")

    @classmethod
    def _from_demult_files(cls, /, *,
                           key: str,
                           fq_dirs: Iterable[str],
                           phred_enc: int):
        """ Yield FastqUnits for each given directory of demultiplexed
        FASTQ files of single-end or interleaved paired-end reads. """
        for fq_dir in fq_dirs:
            yield from cls._from_files(key=key,
                                       fq_files=pathlib.Path(fq_dir).iterdir(),
                                       phred_enc=phred_enc,
                                       one_ref=True)

    @classmethod
    def _from_demult_pairs(cls, /, *, fq_dirs: Iterable[str], phred_enc: int):
        """
        Yield a FastqUnit for each pair of mated, demultiplexed FASTQ files
        in a directory.

        Parameters
        ----------
        phred_enc: int
            Phred score encoding
        fq_dirs: tuple[pathlib.Path]
            Directories containing the FASTQ files

        Return
        ------
        Iterable[FastqUnit]
            One for each FASTQ pair in the directory
        """
        for fq_dir in map(pathlib.Path, fq_dirs):
            # Create empty dictionaries to store the FASTQ files for
            # mates 1 and 2, keyed by the name of the reference sequence
            # for each FASTQ file.
            mates1: dict[str, str] = dict()
            mates2: dict[str, str] = dict()
            # Process every FASTQ file in the directory.
            for fq_file in map(str, fq_dir.iterdir()):
                # Determine if the file is a mate 1 or 2 FASTQ file.
                # Note that using is1 = fq_file.suffix in path.FQ1_EXTS
                # would not work because each item in FQ1_EXTS includes
                # the part of the path that indicates whether the file
                # is mate 1 or 2, which comesbefore the true extension.
                # For example, the file 'mysample_R1.fq' would match the
                # item '_R1.fq' in FQ1_EXTS, but fq_file.suffix == '.fq'
                is1 = any(fq_file.endswith(ext) for ext in path.FQ1_EXTS)
                is2 = any(fq_file.endswith(ext) for ext in path.FQ2_EXTS)
                if is1 and is2:
                    # There should be no way for this error to happen,
                    # but catching it just in case.
                    logger.critical(f"FASTQ path matched both mates: {fq_file}")
                    continue
                if is1:
                    # The file name matched an extension for a mate 1
                    # FASTQ file.
                    fq_path = path.OneRefReads1InFilePath.parse(fq_file)
                    if fq_path.ref in mates1:
                        logger.critical(
                            f"Got >1 mate 1 FASTQ file for ref '{fq_path.ref}'")
                        continue
                    # Add the path to the dict of mate 1 files, keyed
                    # by reference.
                    mates1[fq_path.ref] = fq_file
                elif is2:
                    # The file name matched an extension for a mate 2
                    # FASTQ file.
                    fq_path = path.OneRefReads2InFilePath.parse(fq_file)
                    if fq_path.ref in mates2:
                        logger.critical(
                            f"Got >1 mate 2 FASTQ file for ref '{fq_path.ref}'")
                        continue
                    # Add the path to the dict of mate 2 files, keyed
                    # by reference.
                    mates2[fq_path.ref] = fq_file
                else:
                    # If a file name does not match the expected FASTQ
                    # name format, log a message but keep going, since
                    # the presence or absence of one FASTQ file will
                    # not affect the others, and this file might be an
                    # extraneous file, such as a FASTQC report or a
                    # .DS_Store file on macOS.
                    logger.debug(f"Skipping non-FASTQ file: '{fq_file}'")
            # Iterate through all references from mate 1 and mate 2 files.
            for ref in set(mates1) | set(mates2):
                fastq1 = mates1.get(ref)
                fastq2 = mates2.get(ref)
                if fastq1 is not None and fastq2 is not None:
                    # Yield a new FastqUnit from the paired FASTQ files.
                    fq_paths = {cls.KEYF_MATE1: fastq1, cls.KEYF_MATE2: fastq2}
                    yield cls(**fq_paths, phred_enc=phred_enc, one_ref=True)
                elif fastq1 is None:
                    logger.critical(f"Missing mate 1 for reference '{ref}'")
                else:
                    logger.critical(f"Missing mate 2 for reference '{ref}'")

    @classmethod
    def _from_sample_files(cls, /, *,
                           key: str,
                           fq_files: Iterable[str],
                           phred_enc: int):
        """ Yield a FastqUnit for each given path to a FASTQ file of
        single-end or interleaved paired-end reads from a sample. """
        yield from cls._from_files(key=key,
                                   fq_files=map(pathlib.Path, fq_files),
                                   phred_enc=phred_enc,
                                   one_ref=False)

    @classmethod
    def _from_sample_pairs(cls, /, *,
                           fq1_files: Iterable[str],
                           fq2_files: Iterable[str],
                           phred_enc: int):
        for fq1_file, fq2_file in zip(fq1_files, fq2_files, strict=True):
            fq_paths = {cls.KEYF_MATE1: pathlib.Path(fq1_file),
                        cls.KEYF_MATE2: pathlib.Path(fq2_file)}
            try:
                yield FastqUnit(**fq_paths, phred_enc=phred_enc, one_ref=False)
            except (TypeError, ValueError) as error:
                logger.critical("Failed to pair up FASTQ files of mate 1 reads "
                                f"({fq1_file}) and mate 2 reads ({fq2_file}) "
                                f"due to the following error: {error}")

    @classmethod
    def from_strs(cls, /, *, phred_enc: int, **fastq_args: tuple[str, ...]):
        """
        Yield a FastqUnit for each FASTQ file (or each pair of mate 1
        and mate 2 FASTQ files) whose paths are given as strings.

        Parameters
        ----------
        phred_enc: int
            ASCII offset for encoding Phred scores
        fastq_args: tuple[str]
            FASTQ files, given as tuples of file path string:
            - fastqs: FASTQ files of single-end reads
            - fastqi: FASTQ files of interleaved paired-end reads
            - fastq1: FASTQ files of mate 1 paired-end reads; must
                      correspond 1-for-1 (in order) with fastq2
            - fastq2: FASTQ files of mate 2 paired-end reads; must
                      correspond 1-for-1 (in order) with fastq1
            - fastqs_dir: Directory of FASTQ files of single-end reads
            - fastqi_dir: Directory of FASTQ files of interleaved
                          paired-end reads
            - fastq12_dir: Directory of FASTQ files of separate mate 1
                           and mate 2 paired-end reads; for every FASTQ
                           file of mate 1 reads, there must be a FASTQ
                           file of mate 2 reads with the same sample
                           name, and vice versa.

        Yield
        -----
        FastqUnit
            FastqUnit representing the FASTQ or pair of FASTQ files.
            The order is determined primarily by the order of keyword
            arguments; within each keyword argument, by the order of
            file or directory paths; and for directories, by the order
            in which ```os.path.listdir``` returns file paths.
        """
        keys = set()
        # Directories of single-end and interleaved paired-end FASTQs
        for key in (cls.KEYD_SINGLE, cls.KEYD_INTER):
            keys.add(key)
            fq_dirs = fastq_args.get(key, ())
            yield from cls._from_demult_files(fq_dirs=fq_dirs,
                                              phred_enc=phred_enc,
                                              key=key)
        # Directories of separate mate 1 and mate 2 FASTQs
        keys.add(key := cls.KEYD_MATE12)
        fq_dirs = fastq_args.get(key, ())
        yield from cls._from_demult_pairs(fq_dirs=fq_dirs,
                                          phred_enc=phred_enc)
        # FASTQ files of single-end and interleaved paired-end reads
        for key in (cls.KEYF_SINGLE, cls.KEYF_INTER):
            keys.add(key)
            fq_files = fastq_args.get(key, ())
            yield from cls._from_sample_files(fq_files=fq_files,
                                              phred_enc=phred_enc,
                                              key=key)
        # FASTQ files of separate mate 1 and mate 2 paired-end reads
        keys.add(key := cls.KEYF_MATE1)
        fq1_files = fastq_args.get(key, ())
        keys.add(key := cls.KEYF_MATE2)
        fq2_files = fastq_args.get(key, ())
        yield from cls._from_sample_pairs(fq1_files=fq1_files,
                                          fq2_files=fq2_files,
                                          phred_enc=phred_enc)
        if extras := set(fastq_args) - keys:
            raise TypeError(f"Got extra keyword arguments: {extras}")

    def __str__(self):
        return f"{self.kind} {' and '.join(map(str, self.inputs.values()))}"


def run_fastqc(fq_unit: FastqUnit, out_dir: pathlib.Path, extract: bool):
    """ Run FASTQC on the given FASTQ unit. """
    # FASTQC command, including whether to extract automatically
    cmd = [FASTQC_CMD, "--extract" if extract else "--noextract"]
    # FASTQC output directory
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd.extend(["-o", out_dir])
    # Input FASTQ files
    cmd.extend(fq_unit.inputs.values())
    # Run FASTQC
    logger.info(f"FASTQC of {fq_unit}")
    process = run_cmd(cmd, capture_output=True)
    log_process(logger, process)


def run_cutadapt(fq_inp: FastqUnit,
                 fq_out: FastqUnit, *,
                 n_procs: int,
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
                 cut_m: int):
    """ Trim adapters and low-quality bases with Cutadapt. """
    logger.info(f"Trimming {fq_inp} to {fq_out}")
    # Cutadapt command
    cmd = [CUTADAPT_CMD, "--cores", n_procs]
    # Low-quality trimming
    if cut_nextseq:
        if cut_q1 > 0:
            cmd.extend(["--nextseq-trim", cut_q1])
    else:
        if cut_q1 > 0:
            cmd.extend(["-q", cut_q1])
        if cut_q2 > 0:
            cmd.extend(["-Q", cut_q2])
    # Adapter trimming
    adapters = {"g": cut_g1, "a": cut_a1,
                "G": cut_g2, "A": cut_a2}
    for arg, adapter in adapters.items():
        if adapter and (fq_inp.paired or arg.islower()):
            for adapt in adapter:
                cmd.extend([f"-{arg}", adapt])
    cmd.extend(["-O", cut_o])
    cmd.extend(["-e", cut_e])
    cmd.extend(["-m", cut_m])
    if not cut_indels:
        cmd.append("--no-indels")
    if cut_discard_trimmed:
        cmd.append("--discard-trimmed")
    if cut_discard_untrimmed:
        cmd.append("--discard-untrimmed")
    # FASTQ format
    if fq_inp.interleaved:
        cmd.append("--interleaved")
    # Output files
    for flag, value in zip(("-o", "-p"), fq_out.paths, strict=False):
        cmd.extend([flag, value])
    # Input files
    cmd.extend(fq_inp.cutadapt_input_args)
    # Make the output directory.
    fq_out.parent.mkdir(parents=True, exist_ok=True)
    # Run Cutadapt.
    process = run_cmd(cmd, capture_output=True)
    log_process(logger, process)


def run_bowtie2(fq_inp: FastqUnit,
                index_pfx: pathlib.Path,
                sam_out: pathlib.Path, *,
                n_procs: int,
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
                bt2_orient: str):
    """ Align reads to the reference with Bowtie 2. """
    logger.info(f"Aligning {fq_inp} with reference {index_pfx} to {sam_out}")
    # Bowtie2 command
    cmd = [BOWTIE2_CMD]
    # Resources
    cmd.extend(["--threads", n_procs])
    # Alignment
    if bt2_local:
        cmd.append("--local")
    cmd.extend(["--gbar", bt2_gbar])
    cmd.extend(["--dpad", bt2_dpad])
    cmd.extend(["-L", bt2_l])
    cmd.extend(["-i", bt2_s])
    cmd.extend(["-D", bt2_d])
    cmd.extend(["-R", bt2_r])
    # Scoring
    cmd.append(fq_inp.phred_arg)
    cmd.append("--ignore-quals")
    cmd.extend(["--ma", MATCH_BONUS])
    cmd.extend(["--mp", MISMATCH_PENALTY])
    cmd.extend(["--np", N_PENALTY])
    cmd.extend(["--rfg", REF_GAP_PENALTY])
    cmd.extend(["--rdg", READ_GAP_PENALTY])
    # Filtering
    if not bt2_unal:
        cmd.append("--no-unal")
    cmd.extend(["--score-min", bt2_score_min])
    cmd.extend(["-I", bt2_i])
    cmd.extend(["-X", bt2_x])
    # Mate pair orientation
    orientations = tuple(op.value for op in cli.MateOrientationOption)
    if bt2_orient in orientations:
        cmd.append(f"--{bt2_orient}")
    else:
        cmd.append(f"--{orientations[0]}")
        logger.warning(f"Invalid mate orientation: '{bt2_orient}'; "
                       f"defaulting to '{orientations[0]}'")
    if not bt2_discordant:
        cmd.append("--no-discordant")
    if not bt2_contain:
        cmd.append("--no-contain")
    if bt2_dovetail:
        cmd.append("--dovetail")
    if not bt2_mixed:
        cmd.append("--no-mixed")
    # Formatting
    cmd.append("--xeq")
    # Metrics
    cmd.extend(["--met-stderr", "--met", METRICS_INTERVAL])
    # Input and output files
    cmd.extend(["-S", sam_out])
    cmd.extend(["-x", index_pfx])
    cmd.extend(fq_inp.bowtie2_inputs)
    # Make the output directory.
    sam_out.parent.mkdir(parents=True, exist_ok=True)
    # Run alignment.
    process = run_cmd(cmd, capture_output=True)
    log_process(logger, process)


def dedup_sam(sam_inp: pathlib.Path, sam_out: pathlib.Path):
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
