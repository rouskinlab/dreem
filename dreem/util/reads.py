from abc import ABC, abstractmethod
from collections import namedtuple
import itertools
import logging
import os
import re
from functools import cached_property
from typing import BinaryIO

from dreem.util import path
from dreem.util.cli import DEFAULT_LOCAL, DEFAULT_UNALIGNED, DEFAULT_DISCORDANT, DEFAULT_MIXED, DEFAULT_DOVETAIL, \
    DEFAULT_CONTAIN, DEFAULT_FRAG_LEN_MIN, DEFAULT_FRAG_LEN_MAX, DEFAULT_N_CEILING, DEFAULT_SEED_INTERVAL, \
    DEFAULT_GAP_BAR, DEFAULT_SEED_SIZE, DEFAULT_EXTENSIONS, DEFAULT_RESEED, DEFAULT_PADDING, DEFAULT_ALIGN_THREADS, \
    MATCH_BONUS, MISMATCH_PENALTY, N_PENALTY, REF_GAP_PENALTY, READ_GAP_PENALTY, IGNORE_QUALS
from dreem.util.cli import DEFAULT_MIN_BASE_QUALITY, DEFAULT_ILLUMINA_ADAPTER, DEFAULT_MIN_OVERLAP, DEFAULT_MAX_ERROR, \
    DEFAULT_INDELS, DEFAULT_NEXTSEQ_TRIM, DEFAULT_DISCARD_TRIMMED, DEFAULT_DISCARD_UNTRIMMED, DEFAULT_MIN_LENGTH, \
    DEFAULT_SCORE_MIN
from dreem.util.dflt import BUFFER_LENGTH, NUM_PROCESSES
from dreem.util.excmd import FASTQC_CMD, CUTADAPT_CMD, BOWTIE2_CMD, \
    BOWTIE2_BUILD_CMD, SAMTOOLS_CMD, run_cmd
from dreem.util.seq import FastaParser

# General parameters
DEFAULT_INTERLEAVED = False
SAM_HEADER = b"@"
SAM_ALIGN_SCORE = b"AS:i:"
SAM_EXTRA_SCORE = b"XS:i:"
FASTQ_REC_LENGTH = 4
DEFAULT_MIN_MAPQ = 30

# FastQC parameters
DEFAULT_EXTRACT = False


class FastqUnit(object):
    def __init__(self,
                 inputs: tuple[path.SampleReadsInFilePath |
                               path.DemultReadsInFilePath] |
                         tuple[path.SampleReads1InFilePath,
                               path.SampleReads2InFilePath] |
                         tuple[path.DemultReads1InFilePath,
                               path.DemultReads2InFilePath],
                 interleaved: bool,
                 phred_enc: int):
        self._inputs = inputs
        self._interleaved = interleaved
        self._phred_enc = phred_enc

    @property
    def n_files(self):
        if 1 <= (n_files := len(self._inputs)) <= 2:
            return n_files
        raise TypeError(self._inputs)

    @property
    def interleaved(self):
        if self._interleaved and self.n_files != 1:
            raise TypeError(self._inputs)
        return self._interleaved

    @property
    def paired(self):
        return self.interleaved or self.n_files == 2

    @property
    def demult(self):
        input_types = tuple(map(type, self._inputs))
        if input_types in [(path.DemultReadsInFilePath,),
                           (path.DemultReads1InFilePath,
                            path.DemultReads2InFilePath)]:
            return True
        if input_types in [(path.SampleReadsInFilePath,),
                           (path.SampleReads1InFilePath,
                            path.SampleReads2InFilePath)]:
            return False
        raise TypeError(input_types)

    @property
    def phred_enc(self):
        if 0 <= self._phred_enc < 128:
            return self._phred_enc
        raise ValueError(self._phred_enc)

    @property
    def phred_arg(self):
        return f"--phred{self.phred_enc}"

    @property
    def paths(self):
        return self._inputs

    @cached_property
    def sample(self):
        if self.n_files == 1:
            return self.paths[0].sample
        if self.n_files == 2:
            if ((sample1 := self.paths[0].sample)
                    != (sample2 := self.paths[1].sample)):
                raise ValueError(f"Samples differ: '{sample1}' and '{sample2}'")
            return sample1
        raise ValueError(self.n_files)

    @cached_property
    def ref(self):
        if self.n_files == 1:
            if isinstance(self.paths[0], path.DemultReadsInFilePath):
                return self.paths[0].ref
        if self.n_files == 2:
            path1, path2 = self.paths
            if (isinstance(path1, path.DemultReads1InFilePath)
                    and isinstance(path2, path.DemultReads2InFilePath)):
                if (ref1 := path1.ref) != (ref2 := path2.ref):
                    raise ValueError(f"Refs differ: '{ref1}' and '{ref2}'")
                return path1.ref
        raise TypeError(self.paths)

    @property
    def cutadapt_input_args(self):
        return self.paths

    @property
    def _bowtie2_flags(self):
        if self.n_files == 2:
            return "-1", "-2"
        if self.interleaved:
            return "--interleaved",
        return "-U",

    @property
    def bowtie2_inputs(self):
        return tuple(itertools.chain(*map(list, zip(self._bowtie2_flags,
                                                    self.paths,
                                                    strict=True))))

    @classmethod
    def wrap(cls, *,
             phred_enc: int,
             **fastq_args):
        """
        Return a new FastqUnit instance by formatting the FASTQ file arguments
        into the first input argument for FastqUnit, and verifying that the
        set of FASTQ arguments given was valid. Exactly one of the following
        sets of FASTQ arguments must be given (extras raise an error):
        - fastqs (for one FASTQ file of single-end reads)
        - fastqi (for one FASTQ file of interleaved, paired-end reads)
        - fastq1, fastq2 (for one file each of paired-end reads, mates 1 and 2)

        Parameters
        ----------
        fastqs: SampleReadsInFilePath | OneRefReadsInFilePath | None
            Path to a FASTQ file of single-end reads, or None if no such file.
        fastqi: SampleReadsInFilePath | OneRefReadsInFilePath | None
            Path to a FASTQ file of paired-end reads interleaved in one file,
            or None if no such file.
        fastq1: SampleReads1InFilePath | OneRefReads1InFilePath | None
            Path to a FASTQ file of first mates from paired-end reads, or None
            if no such file. Note: If given, must also give a matched fastq2.
        fastq2: SampleReads2InFilePath | OneRefReads2InFilePath | None
            Path to a FASTQ file of second mates from paired-end reads, or None
            if no such file. Note: If given, must also give a matched fastq1.
        phred_enc: int
            The offset for encoding Phred scores as ASCII characters in the
            FASTQ file(s). Modern Illumina sequencers use +33, while older ones
            have used +64.

        Returns
        -------
        FastqUnit
            A new FastqUnit instance from the given arguments.
        """
        """
        Return the arguments 'inputs' and 'interleaved' for FastqUnit.__init__
        given the values of each FASTQ argument: 'fastqs', 'fastqi', 'fastq1',
        and 'fastq2'.

        Parameters
        ----------
        fq_args: tuple[SampleReadsInFilePath, None, None, None] |
             tuple[OneRefReadsInFilePath, None, None, None] |
             tuple[None, SampleReadsInFilePath, None, None] |
             tuple[None, OneRefReadsInFilePath, None, None] |
             tuple[None, None, SampleReads1InFilePath,
                               SampleReads2InFilePath] |
             tuple[None, None, OneRefReads1InFilePath,
                               OneRefReads2InFilePath]
            Tuple of the four arguments fastqs, fastqi, fastq1, and fastq2
            (explained in the docstring of FastqUnit.wrap), in that order.

        Returns
        -------
        tuple
            Tuple of one or two FASTQ input files, depending on fq_args:
            - If fastqs given, then (fastqs,)
            - If fastqi given, then (fastqi,)
            - If fastq1 and fastq2 given, then (fastq1, fastq2)
        bool
            Whether the FASTQ file contains interleaved paired-end reads.
            - If fastqi given, then True
            - If fastqs or fastq1 and fastq2 given, then False

        Raises
        ------
        TypeError
            If no FASTQ files are given, fastq1 is given without fastq2 (or vice
            versa), or more than one FASTQ file is given (except if both fastq1
            and fastq2 are given).
        """
        match fastq_args:
            case {"fastqs": path.SampleReadsInFilePath() |
                            path.DemultReadsInFilePath() as fastqs}:
                inputs = fastqs,
                interleaved = False
            case {"fastqi": path.SampleReadsInFilePath() |
                            path.DemultReadsInFilePath() as fastqi}:
                inputs = fastqi,
                interleaved = True
            case {"fastq1": path.SampleReads1InFilePath() |
                            path.DemultReads1InFilePath() as fastq1,
                  "fastq2": path.SampleReads2InFilePath() |
                            path.DemultReads2InFilePath() as fastq2}:
                inputs = fastq1, fastq2
                interleaved = False
            case _:
                raise TypeError(fastq_args)
        return FastqUnit(inputs, interleaved, phred_enc)


def _get_demultiplexed_mates(fq_dir: path.SampleInDirPath):
    mates1: dict[str, path.DemultReads1InFilePath] = dict()
    mates2: dict[str, path.DemultReads2InFilePath] = dict()
    for fq_file in os.listdir(fq_dir.path):
        fq_path = fq_dir.path.joinpath(fq_file)
        is1 = any(fq_file.endswith(ext) for ext in path.FQ1_EXTS)
        is2 = any(fq_file.endswith(ext) for ext in path.FQ2_EXTS)
        if is1 and is2:
            raise ValueError(f"FASTQ path matched both mates: '{fq_path}'")
        if is1:
            fq = path.DemultReads1InFilePath.parse_path(fq_path)
            if fq.ref in mates1:
                raise ValueError(f"Got >1 FASTQ file for ref '{fq.ref}'")
            mates1[fq.ref] = fq
        elif is2:
            fq = path.DemultReads2InFilePath.parse_path(fq_path)
            if fq.ref in mates2:
                raise ValueError(f"Got >1 FASTQ file for ref '{fq.ref}'")
            mates2[fq.ref] = fq
        else:
            raise ValueError(f"FASTQ matched no valid format: '{fq_file}'")
    return mates1, mates2


def _get_mate_pairs_by_ref(mates1: dict[str, path.DemultReads1InFilePath],
                           mates2: dict[str, path.DemultReads2InFilePath],
                           phred_enc: int) -> dict[str, FastqUnit]:
    refs1 = set(mates1)
    refs2 = set(mates2)
    refs = refs1 | refs2
    if missing := refs - (refs1 & refs2):
        raise ValueError(f"Missing mate 1 or 2 for refs {', '.join(missing)}")
    pairs: dict[str, FastqUnit] = dict()
    for ref in refs:
        pairs[ref] = FastqUnit.wrap(fastq1=mates1[ref],
                                    fastq2=mates2[ref],
                                    phred_enc=phred_enc)
    return pairs


def get_demultiplexed_fastq_pairs(fq_dir: path.SampleInDirPath, phred_enc: int):
    mates1, mates2 = _get_demultiplexed_mates(fq_dir)
    return _get_mate_pairs_by_ref(mates1, mates2, phred_enc)


class ReadsFileBase(ABC):
    partition = path.TEMP_DIR
    module = path.MOD_ALN
    step = ""
    ext = ""

    def __init__(self, top_dir: path.TopDirPath) -> None:
        self.top = top_dir

    @property
    @abstractmethod
    def sample(self):
        raise NotImplementedError

    @property
    def demult(self):
        raise NotImplementedError

    @cached_property
    def output_dir(self):
        fields = {"top": self.top.top,
                  "partition": self.partition,
                  "module": self.module,
                  "sample": self.sample}
        if self.partition == path.TEMP_DIR:
            return path.SampleTempDirPath(**fields, step=self.step)
        if self.partition == path.OUTPUT_DIR:
            return path.SampleOutDirPath(**fields)
        raise ValueError(self.partition)

    @cached_property
    @abstractmethod
    def output(self):
        raise NotImplementedError

    def setup(self):
        self.output_dir.path.mkdir(parents=True, exist_ok=True)

    @abstractmethod
    def run(self):
        raise NotImplementedError

    @abstractmethod
    def clean(self):
        raise NotImplementedError


class FastqBase(ReadsFileBase):
    def __init__(self, top_dir: path.TopDirPath, fastq: FastqUnit):
        super().__init__(top_dir)
        self.fastq = fastq

    @property
    def sample(self):
        return self.fastq.sample

    @property
    def demult(self):
        return self.fastq.demult

    def qc(self, extract: bool = DEFAULT_EXTRACT):
        cmd = [FASTQC_CMD]
        if extract:
            cmd.append("--extract")
        cmd.extend(map(str, self.fastq.paths))
        run_cmd(cmd)


class FastqTrimmer(FastqBase):
    step = path.ALN_TRIM

    @cached_property
    def _output_fastqs(self):
        fields = self.output_dir.dict()
        return tuple(path.ReadsInToReadsTemp.map_inst(fq, **fields,
                                                      preserve_type=True)
                     for fq in self.fastq.paths)

    @cached_property
    def output(self):
        return FastqUnit(self._output_fastqs,
                         self.fastq.interleaved,
                         self.fastq.phred_enc)

    @property
    def _cutadapt_output_flags(self):
        return "-o", "-p"

    @property
    def _cutadapt_output_args(self):
        return tuple(itertools.chain(*zip(self._cutadapt_output_flags,
                                          self._output_fastqs,
                                          strict=False)))

    def _cutadapt(self,
                  qual1: int = DEFAULT_MIN_BASE_QUALITY,
                  qual2: int = 0,
                  adapters15: tuple[str] = (),
                  adapters13: tuple[str] = (DEFAULT_ILLUMINA_ADAPTER,),
                  adapters25: tuple[str] = (),
                  adapters23: tuple[str] = (DEFAULT_ILLUMINA_ADAPTER,),
                  min_overlap: int = DEFAULT_MIN_OVERLAP,
                  max_error: float = DEFAULT_MAX_ERROR,
                  indels: bool = DEFAULT_INDELS,
                  nextseq_trim: bool = DEFAULT_NEXTSEQ_TRIM,
                  discard_trimmed: bool = DEFAULT_DISCARD_TRIMMED,
                  discard_untrimmed: bool = DEFAULT_DISCARD_UNTRIMMED,
                  min_length: bool = DEFAULT_MIN_LENGTH,
                  cores: int = NUM_PROCESSES):
        cmd = [CUTADAPT_CMD]
        if cores >= 0:
            cmd.extend(["--cores", cores])
        if nextseq_trim:
            if qual1 > 0:
                cmd.extend(["--nextseq-trim", qual1])
        else:
            if qual1 > 0:
                cmd.extend(["-q", qual1])
            if qual2 > 0:
                cmd.extend(["-Q", qual2])
        adapters = {"g": adapters15, "a": adapters13,
                    "G": adapters25, "A": adapters23}
        for arg, adapter in adapters.items():
            if adapter and (self.fastq.paired or arg.islower()):
                for adapt in adapter:
                    cmd.extend([f"-{arg}", adapt])
        if min_overlap >= 0:
            cmd.extend(["-O", min_overlap])
        if max_error >= 0:
            cmd.extend(["-e", max_error])
        if not indels:
            cmd.append("--no-indels")
        if discard_trimmed:
            cmd.append("--discard-trimmed")
        if discard_untrimmed:
            cmd.append("--discard-untrimmed")
        if min_length:
            cmd.extend(["-m", min_length])
        cmd.extend(["--report", "minimal"])
        if self.fastq.interleaved:
            cmd.append("--interleaved")
        cmd.extend(self._cutadapt_output_args)
        cmd.extend(self.fastq.cutadapt_input_args)
        run_cmd(cmd)
        return self.output

    def run(self, **kwargs):
        self.setup()
        return self._cutadapt(**kwargs)

    def clean(self):
        for out_path in self.output.paths:
            out_path.path.unlink()


class FastqAligner(FastqBase):
    step = path.ALN_ALIGN
    ext = path.SAM_EXT

    def __init__(self,
                 top_dir: path.TopDirPath,
                 fastq: FastqUnit,
                 fasta: path.RefsetSeqInFilePath | path.OneRefSeqTempFilePath):
        super().__init__(top_dir, fastq)
        if isinstance(fasta, path.RefsetSeqInFilePath):
            if self.demult:
                raise TypeError("Got a multi-FASTA but a demultiplexed FASTQ")
        elif isinstance(fasta, path.OneRefSeqTempFilePath):
            if not self.demult:
                raise TypeError("Got a single-FASTA but a multiplexed FASTQ")
        else:
            raise TypeError(fasta)
        self.fasta = fasta

    @property
    def fasta_prefix(self):
        return self.fasta.path.with_suffix("")

    @cached_property
    def output(self):
        # fasta provides either 'refset' or 'ref'
        # output_dir provides 'top', 'partition', 'module', 'step', and 'sample'
        # ext provides the file extension
        fields = {**self.fasta.dict(),
                  **self.output_dir.dict(),
                  path.EXT_KEY: self.ext}
        outputs = list()
        for fq in self.fastq.paths:
            temp_inst = path.ReadsInToAlignmentTemp.map_inst(fq, **fields)
            in_type = path.AlignmentInToAlignmentTemp.inverse().map_type(
                type(temp_inst))
            in_inst = in_type.parse_path(temp_inst.path)
            outputs.append(in_inst)
        if not outputs:
            raise ValueError("No output files")
        if any(output != outputs[0] for output in outputs[1:]):
            raise ValueError("Inconsistent output files")
        return outputs[0]

    def _bowtie2_build(self):
        """ Build an index of a reference genome using Bowtie 2. """
        cmd = [BOWTIE2_BUILD_CMD, "-q", self.fasta, self.fasta_prefix]
        run_cmd(cmd)

    def _bowtie2(self,
                 local=DEFAULT_LOCAL,
                 unaligned=DEFAULT_UNALIGNED,
                 discordant=DEFAULT_DISCORDANT,
                 mixed=DEFAULT_MIXED,
                 dovetail=DEFAULT_DOVETAIL,
                 contain=DEFAULT_CONTAIN,
                 score_min=DEFAULT_SCORE_MIN,
                 frag_len_min=DEFAULT_FRAG_LEN_MIN,
                 frag_len_max=DEFAULT_FRAG_LEN_MAX,
                 n_ceil=DEFAULT_N_CEILING,
                 gap_bar=DEFAULT_GAP_BAR,
                 seed_size=DEFAULT_SEED_SIZE,
                 seed_interval=DEFAULT_SEED_INTERVAL,
                 extensions=DEFAULT_EXTENSIONS,
                 reseed=DEFAULT_RESEED,
                 padding=DEFAULT_PADDING,
                 threads=DEFAULT_ALIGN_THREADS):
        cmd = [BOWTIE2_CMD]
        cmd.extend(self.fastq.bowtie2_inputs)
        cmd.extend(["-x", self.fasta_prefix])
        cmd.extend(["-S", self.output])
        cmd.append(self.fastq.phred_arg)
        cmd.append("--xeq")
        cmd.extend(["--ma", MATCH_BONUS])
        cmd.extend(["--mp", MISMATCH_PENALTY])
        cmd.extend(["--np", N_PENALTY])
        cmd.extend(["--rfg", REF_GAP_PENALTY])
        cmd.extend(["--rdg", READ_GAP_PENALTY])
        if score_min:
            cmd.extend(["--score-min", score_min])
        if local:
            cmd.append("--local")
        if not unaligned:
            cmd.append("--no-unal")
        if not discordant:
            cmd.append("--no-discordant")
        if not mixed:
            cmd.append("--no-mixed")
        if dovetail:
            cmd.append("--dovetail")
        if not contain:
            cmd.append("--no-contain")
        if frag_len_min:
            cmd.extend(["-I", frag_len_min])
        if frag_len_max:
            cmd.extend(["-X", frag_len_max])
        if n_ceil:
            cmd.extend(["--n-ceil", n_ceil])
        if gap_bar:
            cmd.extend(["--gbar", gap_bar])
        if seed_size:
            cmd.extend(["-L", seed_size])
        if seed_interval:
            cmd.extend(["-i", seed_interval])
        if extensions:
            cmd.extend(["-D", extensions])
        if reseed:
            cmd.extend(["-R", reseed])
        if padding:
            cmd.extend(["--dpad", padding])
        if threads:
            cmd.extend(["-p", threads])
        if IGNORE_QUALS:
            cmd.append("--ignore-quals")
        run_cmd(cmd)
        return self.output

    def run(self, **kwargs):
        self.setup()
        self._bowtie2_build()
        return self._bowtie2(**kwargs)

    def clean(self):
        self.output.path.unlink()


class XamBase(ReadsFileBase):
    def __init__(self,
                 top_dir: path.TopDirPath,
                 xam: path.RefsetAlignmentInFilePath |
                      path.OneRefAlignmentInFilePath):
        super().__init__(top_dir)
        self.xam = xam

    @property
    def sample(self):
        return self.xam.sample

    @property
    def demult(self):
        if isinstance(self.xam, path.RefsetAlignmentInFilePath):
            return False
        if isinstance(self.xam, path.OneRefAlignmentInFilePath):
            return True
        raise TypeError(self.xam)

    @cached_property
    def output(self):
        if self.partition == path.TEMP_DIR:
            mapper = path.AlignmentInToAlignmentTemp
        elif self.partition == path.OUTPUT_DIR:
            mapper = path.AlignmentInToAlignmentOut
        else:
            raise ValueError(self.partition)
        return mapper.map_inst(self.xam,
                               **self.output_dir.dict(),
                               ext=self.ext,
                               preserve_type=True)

    @cached_property
    def xam_index(self):
        return self.xam.replace(ext=path.BAI_EXT)

    def create_index(self):
        cmd = [SAMTOOLS_CMD, "index", self.xam]
        run_cmd(cmd)
        if not self.xam_index.path.is_file():
            raise FileNotFoundError(self.xam_index.path)
        return self.xam_index

    def clean(self):
        self.output.path.unlink(missing_ok=True)


class SamRemoveEqualMappers(XamBase):
    step = path.ALN_REM
    ext = path.SAM_EXT

    pattern_a = re.compile(SAM_ALIGN_SCORE + rb"(\d+)")
    pattern_x = re.compile(SAM_EXTRA_SCORE + rb"(\d+)")

    _MIN_SAM_FIELDS = 11
    _MAX_SAM_FLAG = 4095  # 2^12 - 1

    @staticmethod
    def _get_score(line: bytes, ptn: re.Pattern[bytes]):
        return (float(match.groups()[0])
                if (match := ptn.search(line)) else None)

    @classmethod
    def _is_best_alignment(cls, line: bytes):
        return ((score_x := cls._get_score(line, cls.pattern_x)) is None
                or score_x < cls._get_score(line, cls.pattern_a))

    @classmethod
    def _read_is_paired(cls, line: bytes):
        info = line.split()
        if len(info) < cls._MIN_SAM_FIELDS:
            raise ValueError(f"Invalid SAM line:\n{line.decode()}")
        flag = int(info[1])
        if 0 <= flag <= cls._MAX_SAM_FLAG:
            return bool(flag % 2)
        raise ValueError(f"Invalid SAM flag: {flag}")

    def _iter_paired(self, sam: BinaryIO, line: bytes):
        for line2 in sam:
            if self._is_best_alignment(line) or self._is_best_alignment(line2):
                yield b"".join((line, line2))
            line = sam.readline()

    def _iter_single(self, sam: BinaryIO, line: bytes):
        while line:
            if self._is_best_alignment(line):
                yield line
            line = sam.readline()

    def _remove_equal_mappers(self, buffer_length=BUFFER_LENGTH):
        with (open(self.xam.path, "rb") as sami,
              open(self.output.path, "wb") as samo):
            # Copy the header from the input to the output SAM file.
            while (line := sami.readline()).startswith(SAM_HEADER):
                samo.write(line)
            if line:
                if self._read_is_paired(line):
                    lines = self._iter_paired(sami, line)
                else:
                    lines = self._iter_single(sami, line)
                while text := b"".join(itertools.islice(lines, buffer_length)):
                    samo.write(text)
        return self.output

    def run(self):
        logging.info("\nRemoving Reads Mapping Equally to Multiple Locations"
                     f" in {self.xam}\n")
        self.setup()
        return self._remove_equal_mappers()


class XamSorter(XamBase):
    def _sort(self, name: bool):
        cmd = [SAMTOOLS_CMD, "sort"]
        if name:
            cmd.append("-n")
        cmd.extend(["-o", self.output, self.xam])
        run_cmd(cmd)
        return self.output

    def run(self, name: bool = False):
        logging.info(f"\nSorting {self.xam} by Reference and Coordinate\n")
        self.setup()
        return self._sort(name)


class BamAlignSorter(XamSorter):
    step = path.ALN_SORT
    ext = path.BAM_EXT


class SamVectorSorter(XamSorter):
    module = path.MOD_VEC
    step = path.VEC_SORT
    ext = path.SAM_EXT


class BamSplitter(XamBase):
    partition = path.OUTPUT_DIR
    step = path.ALN_SPLIT
    ext = path.BAM_EXT

    def __init__(self,
                 top_dir: path.TopDirPath,
                 xam: path.RefsetAlignmentInFilePath |
                      path.OneRefAlignmentInFilePath,
                 fasta: path.RefsetSeqInFilePath |
                        path.OneRefSeqTempFilePath):
        super().__init__(top_dir, xam)
        if isinstance(fasta, path.RefsetSeqInFilePath):
            if self.demult:
                raise TypeError("Got multi-FASTA but demultiplexed BAM")
        elif isinstance(fasta, path.OneRefSeqTempFilePath):
            if not self.demult:
                raise TypeError("Got single-FASTA but multiplexed BAM")
        else:
            raise TypeError(fasta)
        self.fasta = fasta

    @cached_property
    def refs(self):
        return tuple(ref for ref, _ in FastaParser(self.fasta.path).parse())

    def _get_cmd(self, output: path.OneRefAlignmentOutFilePath, ref: str = ""):
        cmd = [SAMTOOLS_CMD, "view"]
        if self.ext == path.BAM_EXT:
            cmd.append("-b")
        cmd.extend(("-o", output, self.xam))
        if ref:
            cmd.append(ref)
        return cmd

    def _extract_one_ref(self, ref: str):
        output = path.OneRefAlignmentOutFilePath(**self.output_dir.dict(),
                                                 ref=ref, ext=self.ext)
        cmd = self._get_cmd(output, ref)
        run_cmd(cmd)
        return output

    def _split_xam(self):
        self.create_index()
        return tuple(map(self._extract_one_ref, self.refs))

    def _move_xam(self):
        if self.ext == self.xam.ext:
            os.rename(self.xam.path, self.output.path)
        else:
            cmd = self._get_cmd(self.output)
            run_cmd(cmd)
        return self.output,

    def run(self) -> tuple[path.OneRefAlignmentOutFilePath, ...]:
        logging.info(f"\nSplitting {self.xam} into Individual References\n")
        self.setup()
        if self.demult:
            return self._move_xam()
        else:
            return self._split_xam()

    def clean(self):
        return


class BamVectorSelector(XamBase):
    module = path.MOD_VEC
    step = path.VEC_SELECT
    ext = path.BAM_EXT

    def __init__(self,
                 top_dir: path.TopDirPath,
                 xam: path.RefsetAlignmentInFilePath |
                      path.OneRefAlignmentInFilePath,
                 ref: str,
                 first: int,
                 last: int):
        super().__init__(top_dir, xam)
        self.ref = ref
        self.first = first
        self.last = last

    @staticmethod
    def ref_coords(ref: str, first: int, last: int):
        return f"{ref}:{first}-{last}"

    def _select(self):
        cmd = [SAMTOOLS_CMD, "view", "-h", "-o", self.output, self.xam,
               self.ref_coords(self.ref, self.first, self.last)]
        run_cmd(cmd)
        return self.output

    def run(self):
        self.setup()
        return self._select()


'''
primer1 = "CAGCACTCAGAGCTAATACGACTCACTATA"
primer1rc = "TATAGTGAGTCGTATTAGCTCTGAGTGCTG"
primer2 = "TGAAGAGCTGGAACGCTTCACTGA"
primer2rc = "TCAGTGAAGCGTTCCAGCTCTTCA"
adapters5 = (primer1, primer2rc)
adapters3 = (primer2, primer1rc)
'''
