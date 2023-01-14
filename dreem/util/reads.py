from functools import cached_property
import itertools
from io import BufferedReader
import logging
import os
import re

from dreem.util.cmd import SAMTOOLS_CMD, run_cmd
from dreem.util.dflt import BUFFER_LENGTH
from dreem.util.seq import FastaParser
from dreem.util.path import BasePath, FastaInPath, XamTempPath

from functools import cached_property
import itertools
import logging
import os

from typing import Tuple, Dict

from dreem.util.cmd import FASTQC_CMD, CUTADAPT_CMD, BOWTIE2_CMD, \
    BOWTIE2_BUILD_CMD, run_cmd
from dreem.util.dflt import NUM_PROCESSES, PHRED_ENCODING
from dreem.util.cli import DEFAULT_MIN_BASE_QUALITY, DEFAULT_ILLUMINA_ADAPTER, DEFAULT_MIN_OVERLAP, DEFAULT_MAX_ERROR, DEFAULT_INDELS, DEFAULT_NEXTSEQ_TRIM, DEFAULT_DISCARD_TRIMMED, DEFAULT_DISCARD_UNTRIMMED, DEFAULT_MIN_LENGTH, DEFAULT_SCORE_MIN
from dreem.util.cli import DEFAULT_LOCAL, DEFAULT_UNALIGNED, DEFAULT_DISCORDANT, DEFAULT_MIXED, DEFAULT_DOVETAIL, DEFAULT_CONTAIN, DEFAULT_FRAG_LEN_MIN, DEFAULT_FRAG_LEN_MAX, DEFAULT_N_CEILING, DEFAULT_SEED_INTERVAL, DEFAULT_GAP_BAR, DEFAULT_SEED_SIZE, DEFAULT_EXTENSIONS, DEFAULT_RESEED, DEFAULT_PADDING, DEFAULT_ALIGN_THREADS, MATCH_BONUS, MISMATCH_PENALTY, N_PENALTY, REF_GAP_PENALTY, READ_GAP_PENALTY, IGNORE_QUALS
from dreem.util.path import SampleTempPath, FastaInPath, XamTempPath, TEMP_DIR, MOD_ALN, MOD_VEC, ALN_TRIM, ALN_ALIGN, ALN_REM, ALN_SORT, ALN_SPLIT, VEC_SELECT, VEC_SORT
from dreem.util.path import FastqInPath, Fastq1InPath, Fastq2InPath, FastqDemultiPath, Fastq1DemultiPath, Fastq2DemultiPath, SampleOutPath, XamTempPath, XamInPath, XamOutPath, XamIndexInPath, XamIndexTempPath, XamIndexSegment
from dreem.util.path import FastqOutPath, Fastq1OutPath, Fastq2OutPath, SAM_EXT, BAM_EXT, OUTPUT_DIR, XamSegment, PathError


# General parameters
DEFAULT_INTERLEAVED = False
SAM_HEADER = b"@"
SAM_ALIGN_SCORE = b"AS:i:"
SAM_EXTRA_SCORE = b"XS:i:"
FASTQ_REC_LENGTH = 4
DEFAULT_MIN_MAPQ = 30

# FastQC parameters
DEFAULT_EXTRACT = False


def _get_demultiplexed_mates(sample_path: SampleOutPath):
    mates1 = dict()
    mates2 = dict()
    for fq_file in os.listdir(sample_path.path):
        fq_path = str(sample_path.path.joinpath(fq_file))
        try:
            fq1 = Fastq1DemultiPath.parse(fq_path)
        except PathError:
            fq1 = None
        try:
            fq2 = Fastq2DemultiPath.parse(fq_path)
        except PathError:
            fq2 = None
        if fq1 and fq2:
            raise ValueError("Failed to determine whether FASTQ is mate 1 or 2:"
                             f" '{fq_path}'")
        if fq1:
            ref = fq1.fastq1.name
            if ref in mates1:
                raise ValueError(f"More than one FASTQ file for ref '{ref}'")
            mates1[ref] = fq1
        elif fq2:
            ref = fq2.fastq2.name
            if ref in mates2:
                raise ValueError(f"More than one FASTQ file for ref '{ref}'")
            mates2[ref] = fq2
        else:
            raise ValueError("FASTQ file name did not match any valid formats:"
                             f" '{fq_path}'")
    return mates1, mates2


def _get_mate_pairs_by_ref(mates1: Dict[str, Fastq1DemultiPath],
                           mates2: Dict[str, Fastq2DemultiPath]):
    refs1 = set(mates1)
    refs2 = set(mates2)
    refs = refs1 | refs2
    if (missing := refs - (refs1 & refs2)):
        raise ValueError("Missing mate 1 or 2 for refs "
                         f"{', '.join(missing)}")
    pairs: Dict[bytes, FastqUnit] = dict()
    for ref in refs:
        mate1 = str(mates1[ref].path)
        mate2 = str(mates2[ref].path)
        pairs[ref] = FastqUnit(fastq1=mate1, fastq2=mate2,
                               output=True, demultiplexed=True)
    return pairs


def get_demultiplexed_fastq_pairs(fq_dir: str):
    sample_path = SampleOutPath.parse(fq_dir)
    mates1, mates2 = _get_demultiplexed_mates(sample_path)
    pairs = _get_mate_pairs_by_ref(mates1, mates2)
    return pairs


class FastqUnit(object):
    _dm_classes = (FastqDemultiPath, FastqDemultiPath,
                   Fastq1DemultiPath, Fastq2DemultiPath)
    _out_classes = (FastqOutPath, FastqOutPath, Fastq1OutPath, Fastq2OutPath)
    _in_classes = (FastqInPath, FastqInPath, Fastq1InPath, Fastq2InPath)

    def __init__(self,
                 fastqu: str = "", fastqi: str = "",
                 fastq1: str = "", fastq2: str = "",
                 output: bool = False, demultiplexed: bool = False,
                 phred_encoding: int = PHRED_ENCODING):
        self._fqu = fastqu
        self._fqi = fastqi
        self._fq1 = fastq1
        self._fq2 = fastq2
        self._out = output
        self._demulti = demultiplexed
        self._phred = phred_encoding
    
    @property
    def _fq_paths(self):
        return {
            "fastqu": self._fqu,
            "fastqi": self._fqi,
            "fastq1": self._fq1,
            "fastq2": self._fq2,
        }

    @property
    def kwargs(self):
        return {
            **self._fq_paths,
            "output": self._out,
            "demultiplexed": self._demulti,
            "phred_encoding": self._phred,
        }

    @property
    def phred(self):
        return f"--phred{self._phred}"
    
    def _get_in_classes(self):
        if self._demulti:
            return self._dm_classes
        if self._out:
            return self._out_classes
        return self._in_classes
    
    def _get_out_classes(self):
        return self._dm_classes if self._demulti else self._out_classes
    
    @cached_property
    def _flags(self):
        is_unpaired = bool(self._fqu)
        is_interleaved = bool(self._fqi)
        has_paired1 = bool(self._fq1)
        has_paired2 = bool(self._fq2)
        if has_paired1 != has_paired2:
            raise ValueError("Cannot give fastq1 or fastq2 without the other")
        if (count := is_unpaired + is_interleaved + has_paired1) != 1:
            if count:
                raise ValueError("Too many FASTQ files given")
            raise ValueError("No FASTQ files given")
        return is_unpaired, is_interleaved, has_paired1, has_paired2
        
    @property
    def paired(self):
        return not self._flags[0]
    
    @property
    def interleaved(self):
        return self._flags[1]
    
    @property
    def two_files(self):
        return self._flags[2]
    
    @property
    def _in_classes_filtered(self):
        return tuple(cls for cls, flag in zip(self._get_in_classes(),
                                              self._flags, strict=True)
                     if flag)

    @property
    def _out_classes_filtered(self):
        return tuple(cls for cls, flag in zip(self._get_out_classes(),
                                              self._flags, strict=True)
                     if flag)
    
    @property
    def _fq_items_filtered(self):
        return tuple(item for item, flag in zip(self._fq_paths.items(),
                                                self._flags, strict=True)
                     if flag)
    
    @property
    def _fq_args_filtered(self):
        return tuple(arg for arg, _ in self._fq_items_filtered)

    @property
    def _fq_paths_filtered(self):
        return tuple(path for _, path in self._fq_items_filtered)
    
    @property
    def fqs(self):
        return [cls.parse(path) for cls, path in
                zip(self._in_classes_filtered, self._fq_paths_filtered,
                    strict=True)]
    
    @property
    def paths(self):
        return [fq.path for fq in self.fqs]
    
    def _get_single_value(self, values: list, name: str):
        if not values:
            raise ValueError(f"No {name} for {self}")
        if len(values) > 2:
            raise ValueError(f"Too many {name}s for {self}")
        if len(values) == 2 and values[0] != values[1]:
            raise ValueError(f"Inconsistent {name}s for {self}")
        return values[0]

    @property
    def sample(self) -> str:
        return self._get_single_value(
            [fq.last.name for fq in self.fqs], "sample")
    
    @property
    def ref(self) -> bytes:
        if self._demulti:
            return self._get_single_value(
                [fq.last.name for fq in self.fqs], "ref")
        raise ValueError("Cannot get ref from non-demultiplexed FASTQs")
    
    @property
    def cutadapt_input_args(self):
        return self.paths
    
    @property
    def cutadapt_output_args(self):
        args = list()
        for flag, output in zip(("-o", "-p"), self.paths, strict=False):
            args.extend([flag, output])
        return args

    @property
    def bowtie2_input_args(self):
        unpaired, interleaved, mate1, mate2 = self._flags
        if unpaired:
            flags = ("-U",)
        elif interleaved:
            flags = ("--interleaved",)
        else:
            assert mate1 and mate2
            flags = ("-1", "-2")
        return list(itertools.chain(*zip(flags, self.paths, strict=True)))
    
    def rename(self, **kwargs: str):
        renamed = dict()
        for fq, fq_arg, cls in zip(self.fqs, self._fq_args_filtered,
                                   self._out_classes_filtered, strict=True):
            fq_kwargs = {**fq.kwargs, **kwargs}
            if self._demulti:
                fq_kwargs["sample"] = self.sample
            renamed[fq_arg] = cls(**fq_kwargs).path
        return self.__class__(**renamed)


class ReadsFileBase(object):
    partition = TEMP_DIR
    module = MOD_ALN
    _step = ""

    def __init__(self, base_path: BasePath) -> None:
        self.base_path = base_path
    
    @property
    def step(self):
        if self._step:
            return self._step
        raise NotImplementedError
    
    def setup(self):
        raise NotImplementedError

    def run(self):
        raise NotImplementedError
    
    def clean(self):
        raise NotImplementedError


class FastqBase(ReadsFileBase):
    def __init__(self, base_path: BasePath, fq_unit: FastqUnit) -> None:
        super().__init__(base_path)
        self.fq_unit = fq_unit
    
    @property
    def paired(self):
        return self.fq_unit.paired

    @property
    def _path_args(self):
        return {
            "base": self.base_path,
            "partition": self.partition,
            "module": self.module,
            "temp_step": self.step,
            "sample": self.fq_unit.sample
        }
    
    @property
    def output_dir(self):
        return SampleTempPath(**self._path_args)
    
    def qc(self, extract: bool = DEFAULT_EXTRACT):
        cmd = [FASTQC_CMD]
        if extract:
            cmd.append("--extract")
        cmd.extend(self.fq_unit.paths)
        run_cmd(cmd)

    def setup(self):
        self.output_dir.path.mkdir(parents=True, exist_ok=True)


class FastqTrimmer(FastqBase):
    _step = ALN_TRIM

    @property
    def output(self):
        return self.fq_unit.rename(**self._path_args)

    def _cutadapt(self,
                  qual1: int=DEFAULT_MIN_BASE_QUALITY,
                  qual2: int=0,
                  adapters15: Tuple[str]=(),
                  adapters13: Tuple[str]=(DEFAULT_ILLUMINA_ADAPTER,),
                  adapters25: Tuple[str]=(),
                  adapters23: Tuple[str]=(DEFAULT_ILLUMINA_ADAPTER,),
                  min_overlap: int=DEFAULT_MIN_OVERLAP,
                  max_error: float=DEFAULT_MAX_ERROR,
                  indels: bool=DEFAULT_INDELS,
                  nextseq_trim: bool=DEFAULT_NEXTSEQ_TRIM,
                  discard_trimmed: bool=DEFAULT_DISCARD_TRIMMED,
                  discard_untrimmed: bool=DEFAULT_DISCARD_UNTRIMMED,
                  min_length: bool=DEFAULT_MIN_LENGTH,
                  cores: int=NUM_PROCESSES):
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
            if adapter and (self.paired or arg.islower()):
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
        if self.fq_unit.interleaved:
            cmd.append("--interleaved")
        cmd.extend(self.output.cutadapt_output_args)
        cmd.extend(self.fq_unit.cutadapt_input_args)
        run_cmd(cmd)
        return self.output
    
    def run(self, **kwargs):
        self.setup()
        return self._cutadapt(**kwargs)
    
    def clean(self):
        for path in self.output.paths:
            path.unlink()


class FastqAligner(FastqBase):
    _step = ALN_ALIGN

    def __init__(self, base_dir: BasePath,
                 fastq: BasePath, refs_file: FastaInPath) -> None:
        super().__init__(base_dir, fastq)
        self.refs_file = refs_file
    
    @property
    def refs_name(self):
        return self.refs_file.path.stem
    
    @property
    def refs_prefix(self):
        return self.refs_file.path.with_suffix("")

    @property
    def output(self):
        kwargs = {**self._path_args,
                  "xam": XamSegment.format(self.refs_name)}
        return XamTempPath(**kwargs)

    def _bowtie2_build(self):
        """
        Build an index of a reference genome using Bowtie 2.
        :param ref: (str) path to the reference genome FASTA file
        :return: None
        """
        cmd = [BOWTIE2_BUILD_CMD, "-q", self.refs_file, self.refs_prefix]
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
        cmd.extend(self.fq_unit.bowtie2_input_args)
        cmd.extend(["-x", self.refs_prefix])
        cmd.extend(["-S", self.output])
        cmd.append(self.fq_unit.phred)
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
    ext = SAM_EXT

    def __init__(self, base_path: BasePath, xam_in: XamInPath | XamTempPath) -> None:
        super().__init__(base_path)
        self.xam_in = xam_in
    
    @property
    def xai_class(self):
        if isinstance(self.xam_in, XamInPath):
            return XamIndexInPath
        if isinstance(self.xam_in, XamTempPath):
            return XamIndexTempPath
        raise TypeError(self.xam_in)
    
    @classmethod
    def out_class(cls):
        if cls.partition == TEMP_DIR:
            return XamTempPath
        if cls.partition == OUTPUT_DIR:
            return XamOutPath
        raise NotImplementedError
    
    @property
    def sample(self):
        return self.xam_in.sample

    @property
    def name_in(self):
        return self.xam_in.xam.name
    
    @property
    def name_out(self):
        return self.xam_out.xam.name
    
    def _get_xam_out(self, name: str):
        kwargs = {
            "base": self.base_path,
            "partition": self.partition,
            "module": self.module,
            "sample": self.sample
        }
        if self.partition == TEMP_DIR:
            kwargs["temp_step"] = self.step
        kwargs["xam"] = XamSegment.format(name, ext=self.ext)
        return self.out_class()(**kwargs)
    
    @cached_property
    def xam_out(self):
        return self._get_xam_out(self.name_in)
    
    @property
    def output_dir(self):
        return self.xam_out.path.parent
    
    def _make_output_dir(self):
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    @property
    def xam_index(self):
        xai_seg = XamIndexSegment.format(self.xam_in.last.name)
        args = self.xam_in.args[:-1] + (xai_seg,)
        return self.xai_class(*args)

    def create_index(self):
        cmd = [SAMTOOLS_CMD, "index", self.xam_in]
        run_cmd(cmd)
        xam_index = self.xam_index
        assert xam_index.path.is_file()
        return xam_index
    
    def clean(self):
        self.xam_out.path.unlink()


class SamRemoveEqualMappers(XamBase):
    step = ALN_REM

    pattern_a = re.compile(SAM_ALIGN_SCORE + rb"(\d+)")
    pattern_x = re.compile(SAM_EXTRA_SCORE + rb"(\d+)")

    @staticmethod
    def _get_score(line: bytes, ptn: re.Pattern[bytes]):
        return (float(match.groups()[0])
                if (match := ptn.search(line)) else None)
    
    @classmethod
    def _is_best_alignment(cls, line: bytes):
        return ((score_x := cls._get_score(line, cls.pattern_x)) is None
                or score_x < cls._get_score(line, cls.pattern_a))
    
    @staticmethod
    def _read_is_paired(line: bytes):
        info = line.split()
        if len(info) < 11:  # FIXME: define magic number
            raise ValueError(f"Invalid SAM line:\n{line.decode()}")
        flag = int(info[1])
        if 0 <= flag < 2**12:  # FIXME: define magic number 
            return bool(flag % 2)
        raise ValueError(f"Invalid SAM flag: {flag}")
    
    def _iter_paired(self, sam: BufferedReader, line: bytes):
        for line2 in sam:
            if self._is_best_alignment(line) or self._is_best_alignment(line2):
                yield b"".join((line, line2))
            line = sam.readline()
    
    def _iter_single(self, sam: BufferedReader, line: bytes):
        while line:
            if self._is_best_alignment(line):
                yield(line)
            line = sam.readline()

    def _remove_equal_mappers(self, buffer_length=BUFFER_LENGTH):
        with (open(self.xam_in.path, "rb") as sami,
              open(self.xam_out.path, "wb") as samo):
            # Copy the header from the input to the output SAM file.
            while (line := sami.readline()).startswith(SAM_HEADER):
                samo.write(line)
            if line:
                paired = self._read_is_paired(line)
                iter_sam = self._iter_paired if paired else self._iter_single
                lines = iter_sam(sami, line)
                while text := b"".join(itertools.islice(lines, buffer_length)):
                    samo.write(text)
        return self.xam_out
    
    def run(self):
        logging.info("\nRemoving Reads Mapping Equally to Multiple Locations"
                     f" in {self.xam_in}\n")
        self._make_output_dir()
        return self._remove_equal_mappers()


class XamSorter(XamBase):
    def _sort(self, name: bool):
        cmd = [SAMTOOLS_CMD, "sort"]
        if name:
            cmd.append("-n")
        cmd.extend(["-o", self.xam_out, self.xam_in])
        run_cmd(cmd)
        return self.xam_out
    
    def run(self, name: bool = False):
        logging.info(f"\nSorting {self.xam_in} by Reference and Coordinate\n")
        self._make_output_dir()
        return self._sort(name)


class BamAlignSorter(XamSorter):
    step = ALN_SORT
    ext = BAM_EXT


class SamVectorSorter(XamSorter):
    module = MOD_VEC
    step = VEC_SORT
    ext = SAM_EXT


class BamSplitter(XamBase):
    partition = OUTPUT_DIR
    step = ALN_SPLIT
    ext = BAM_EXT

    def __init__(self, base_path: BasePath,
                 xam_in: XamInPath | XamTempPath,
                 refs_file: FastaInPath) -> None:
        super().__init__(base_path, xam_in)
        self._refs_file = refs_file
    
    @cached_property
    def refs(self):
        return [ref for ref, _ in FastaParser(self._refs_file.path).parse()]
    
    def _extract_one_ref(self, ref: bytes):
        xam_out = self._get_xam_out(ref)
        cmd = [SAMTOOLS_CMD, "view", "-b", "-o", xam_out, self.xam_in, ref]
        run_cmd(cmd)
        return xam_out
    
    def _split_bam(self):
        self.create_index()
        return list(map(self._extract_one_ref, self.refs)) 
    
    def run(self):
        logging.info(f"\nSplitting {self.xam_in} into Individual References\n")
        self._make_output_dir()
        return self._split_bam()


class BamVectorSelector(XamBase):
    module = MOD_VEC
    step = VEC_SELECT
    ext = BAM_EXT
    
    @staticmethod
    def ref_coords(ref: str, first: int, last: int):
        return f"{ref}:{first}-{last}"

    def _select(self, ref: str, first: int, last: int):
        cmd = [SAMTOOLS_CMD, "view", "-h", "-o", self.xam_out, self.xam_in,
               self.ref_coords(ref, first, last)]
        run_cmd(cmd)
        return self.xam_out
    
    def run(self, ref: str, first: int, last: int):
        logging.info(f"\nSelecting {self.xam_in} ref {ref}: {first}-{last}\n")
        self._make_output_dir()
        return self._select(ref, first, last)



'''
primer1 = "CAGCACTCAGAGCTAATACGACTCACTATA"
primer1rc = "TATAGTGAGTCGTATTAGCTCTGAGTGCTG"
primer2 = "TGAAGAGCTGGAACGCTTCACTGA"
primer2rc = "TCAGTGAAGCGTTCCAGCTCTTCA"
adapters5 = (primer1, primer2rc)
adapters3 = (primer2, primer1rc)
'''
