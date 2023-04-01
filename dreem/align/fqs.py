from functools import cached_property
from itertools import chain
from logging import getLogger
from pathlib import Path
from typing import Iterable

from ..util import path
from ..util.shell import run_cmd, BOWTIE2_CMD, CUTADAPT_CMD, FASTQC_CMD


logger = getLogger(__name__)


# Bowtie2 parameters
# Consider this example: Ref = ACGT, Read = AG
# Assume that we want to minimize the number of edits needed to convert
# the reference into the read sequence. The smallest number of edits is
# two, specifically these two deletions (/) from the reference: [A/G/]
# which gets a score of (2 * match - 2 * gap_open - 2 * gap_extend).
# But there are two alternative alignment, each with 3 edits:
# [Ag//] and [A//g] (substitutions marked in lowercase). Each gets the
# score (match - substitution - gap_open - 2 * gap_extend).
# In order to favor the simpler alignment with two edits,
# (2 * match - 2 * gap_open - 2 * gap_extend) must be greater than
# (match - substitution - gap_open - 2 * gap_extend); this simplifies to
# (substitution > gap_open - match). Thus, the substitution penalty and
# match bonus must be relatively large, and the gap open penalty small.
# We want to avoid introducing too many gaps, especially to prevent the
# introduction of an insertion and a deletion from scoring better than
# one substitution. Consider this example: Ref = ATAT, Read = ACTT
# The simplest alignment (the smallest number of mutations) is ActT,
# which gets a score of (2 * match - 2 * substitution).
# Another alignment with indels is A{C}T/T, where {C} means a C was
# inserted into the read and the / denotes an A deleted from the read.
# This alignment scores (3 * match - 2 * gap_open - 2 * gap_extend).
# Thus, (2 * match - 2 * substitution) must be greater than
# (3 * match - 2 * gap_open - 2 * gap_extend), which simplifies to
# (2 * gap_open + 2 * gap_extend > match + 2 * substitution).
# There are two easy solutions to these inequalities:
# - Bowtie v2.5 defaults: 6 > 5 - 2 and 2*5 + 2*3 > 2 + 2*6
# - Set every value to 1: 1 > 1 - 1 and 2*1 + 2*1 > 1 + 2*1
MATCH_BONUS = "1"
MISMATCH_PENALTY = "1,1"
N_PENALTY = "0"
REF_GAP_PENALTY = "1,1"
READ_GAP_PENALTY = "1,1"
METRICS_INTERVAL = 60  # Write metrics once every 60 seconds.


class FastqUnit(object):
    """
    Unified interface for the following sets of sequencing reads:

    - One FASTQ file of single-end reads from one sample
    - One FASTQ file of interleaved, paired-end reads from one sample
    - Two FASTQ files of mate 1 and 2 paired-end reads from one sample
    - One FASTQ file of single-end reads originating from one reference
      sequence in one sample
    - One FASTQ file of interleaved, paired-end reads originating from
      one reference sequence in one sample
    - Two FASTQ files of mate 1 and mate 2 paired-end reads originating
      from one reference sequence in one sample

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
                 fastqs: Path | None = None,
                 fastqi: Path | None = None,
                 fastq1: Path | None = None,
                 fastq2: Path | None = None,
                 phred_enc: int,
                 one_ref: bool):
        if fastqs:
            if fastqi or fastq1 or fastq2:
                raise TypeError("Got too many FASTQ files")
            self.inputs: dict[str, Path] = {self.KEYF_SINGLE: fastqs}
            self.paired = False
            self.interleaved = False
        elif fastqi:
            if fastq1 or fastq2:
                raise TypeError("Got too many FASTQ files")
            self.inputs: dict[str, Path] = {self.KEYF_INTER: fastqi}
            self.paired = True
            self.interleaved = True
        elif fastq1:
            if not fastq2:
                raise TypeError("Got fastq1 but not fastq2")
            self.inputs: dict[str, Path] = {self.KEYF_MATE1: fastq1,
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
                     + ", ".join(f"{k} = {v} (type '{type(v).__name__}')"
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
    def parent(self) -> Path:
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
        refs: set[str] = set()
        # Get the reference name from each FASTQ that names a reference.
        for fq in self.paths:
            try:
                refs.add(fq.ref)
            except AttributeError:
                pass
        if len(refs) == 0:
            # No FASTQ files named a reference.
            return None
        if len(refs) > 1:
            raise ValueError(f"Reference names of input files {self} disagree: "
                             + " ≠ ".join(refs))
        return list(refs)[0]

    def is_compatible_fasta(self, fasta: Path, one_ref_fasta: bool):
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
        return tuple(chain(*[(self.BOWTIE2_FLAGS[key], fq)
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
                    fq_files: Iterable[Path],
                    one_ref: bool,
                    phred_enc: int):
        """ Yield a FastqUnit for each given path to a FASTQ file of
        single-end or interleaved paired-end reads. """
        for fq in fq_files:
            if fq.suffix not in path.FQ_EXTS:
                logger.warning(f"Skipping non-FASTQ file: '{fq}'")
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
                                       fq_files=Path(fq_dir).iterdir(),
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
        fq_dirs: tuple[Path]
            Directories containing the FASTQ files

        Return
        ------
        Iterable[FastqUnit]
            One for each FASTQ pair in the directory
        """
        for fq_dir in map(Path, fq_dirs):
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
                    logger.warning(f"Skipping non-FASTQ file: '{fq_file}'")
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
                                   fq_files=map(Path, fq_files),
                                   phred_enc=phred_enc,
                                   one_ref=False)

    @classmethod
    def _from_sample_pairs(cls, /, *,
                           fq1_files: Iterable[str],
                           fq2_files: Iterable[str],
                           phred_enc: int):
        for fq1_file, fq2_file in zip(fq1_files, fq2_files, strict=True):
            fq_paths = {cls.KEYF_MATE1: Path(fq1_file),
                        cls.KEYF_MATE2: Path(fq2_file)}
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
            # This should never happen; catching just in case.
            raise TypeError(f"Got extra keyword arguments: {extras}")

    def __str__(self):
        return f"{self.kind} {' and '.join(map(str, self.inputs.values()))}"


def run_fastqc(fq_unit: FastqUnit, out_dir: Path, extract: bool):
    """ Run FASTQC on the given FASTQ unit. """
    logger.info(f"Began FASTQC of {fq_unit}")
    # FASTQC command, including whether to extract automatically
    cmd = [FASTQC_CMD, "--extract" if extract else "--noextract"]
    # FASTQC output directory
    out_dir.mkdir(parents=True, exist_ok=True)
    logger.debug(f"Ensured directory: {out_dir}")
    cmd.extend(["-o", out_dir])
    # Input FASTQ files
    cmd.extend(fq_unit.inputs.values())
    # Run FASTQC
    run_cmd(cmd)
    logger.info(f"Ended FASTQC of {fq_unit}")


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
    logger.info(f"Began trimming {fq_inp}")
    # Cutadapt command
    cmd = [CUTADAPT_CMD, "--cores", n_procs]
    # Quality trimming
    if cut_nextseq:
        cut_qnext = max(cut_q1, cut_q2)
        if cut_q1 != cut_q2:
            logger.warning("NextSeq trimming takes one quality level, but got "
                           f"two ({cut_q1} and {cut_q2}); using {cut_qnext}")
        cmd.extend(["--nextseq-trim", cut_qnext])
    else:
        cmd.extend(["-q", cut_q1, "-Q", cut_q2])
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
    logger.debug(f"Ensured directory: {fq_out.parent}")
    # Run Cutadapt.
    run_cmd(cmd)
    logger.info(f"Ended trimming {fq_inp}; output {fq_out}")


def run_bowtie2(fq_inp: FastqUnit,
                index_pfx: Path,
                sam_out: Path, *,
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
    logger.info(f"Began aligning {fq_inp} to {index_pfx}")
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
    cmd.append(f"--{bt2_orient}")
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
    logger.debug(f"Ensured directory: {sam_out.parent}")
    # Run alignment.
    run_cmd(cmd)
    logger.info(f"Ended aligning {fq_inp}; output {sam_out}")
