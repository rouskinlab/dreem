from __future__ import annotations
from dataclasses import dataclass, asdict, astuple
import itertools
import os
from pathlib import Path
import re
from typing import Callable, Dict, Tuple, Optional, Any


# Partitions
OUTPUT_DIR = "output"
TEMP_DIR = "temp"
PARTITIONS = (OUTPUT_DIR, TEMP_DIR)

# Module directories
MOD_DMX = "demultiplexing"
MOD_ALN = "alignment"
MOD_VEC = "vectoring"
MOD_CLS = "clustering"
MOD_AGG = "aggregation"
MODULES = (MOD_DMX, MOD_ALN, MOD_VEC, MOD_CLS, MOD_AGG)

# Alignment steps
ALN_TRIM = "align_1_trim"
ALN_ALIGN = "align_2_align"
ALN_REM = "align_3_rem"
ALN_SORT = "align_4_sort"
ALN_SPLIT = "align_5_split"
ALN_STEPS = (ALN_TRIM, ALN_ALIGN, ALN_REM, ALN_SORT, ALN_SPLIT)

# Vectoring steps
VEC_SELECT = "vector_1_select"
VEC_SORT = "vector_2_sort"
VEC_STEPS = (VEC_SELECT, VEC_SORT)

TEMP_STEPS = ALN_STEPS + VEC_STEPS

# File extensions
FASTA_EXTS = (".fasta", ".fa")
FQ_EXTS = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
FQ_PATTERNS = ("_r{}{}", "_mate{}{}", "_{}_sequence{}")
FQ_CASES = (str.lower, str.upper)
FQ1_EXTS = tuple(case(ptrn).format(1, ext) for case, ptrn, ext in
                 itertools.product(FQ_CASES, FQ_PATTERNS, FQ_EXTS))
FQ2_EXTS = tuple(case(ptrn).format(2, ext) for case, ptrn, ext in
                 itertools.product(FQ_CASES, FQ_PATTERNS, FQ_EXTS))
SAM_EXT = ".sam"
BAM_EXT = ".bam"
CRAM_EXT = ".cram"
XAM_EXTS = (SAM_EXT, BAM_EXT, CRAM_EXT)
XAI_EXTS = (f"{BAM_EXT}.bai",)


# Path functions

def sanitize(path: Any):
    return os.path.realpath(os.path.normpath(os.path.abspath(str(path))))


def switch_directory(old_path: str, new_dir: str):
    return os.path.join(new_dir, os.path.basename(old_path))


def try_remove(file: str):
    try:
        os.remove(file)
    except OSError:
        pass


# 

class PathError(BaseException):
    pass


# Path segment classes


class PathSegment(object):
    def __init__(self, segment: Any) -> None:
        self._assert_valid(segment)
        self._segment = str(segment)
    
    @classmethod
    def _assert_valid(cls, segment: Any):
        if not cls.segment_is_valid(segment):
            raise PathError(f"Invalid segment for {cls}: '{segment}'")
    
    def __str__(self) -> str:
        return self._segment
    
    def __repr__(self) -> str:
        return f"'{self}'"
    
    @staticmethod
    def segment_is_valid(segment: Any):
        return True


class BaseSegment(PathSegment):
    def __init__(self, segment: Any) -> None:
        super().__init__(sanitize(segment))
    
    segment_is_valid = staticmethod(lambda seg: os.path.isdir(str(seg)))
    
    @property
    def path(self):
        return Path(self._segment)


class SubSegment(PathSegment):
    segment_is_valid = staticmethod(lambda seg: os.sep not in str(seg))
    
    @property
    def name(self):
        return self._segment


class PatternedSegment(SubSegment):
    _pattern: str = "(.+)"
    _dtypes: Tuple[type] = (str,)
    _exts: Tuple[str] = ()
    _fmt = staticmethod(lambda ext: "{}" + ext)

    @classmethod
    def _get_patterns(cls):
        return (re.compile(f"^{cls._pattern}({ext.replace('.', '[.]')})$")
                for ext in cls._exts)

    @classmethod
    def _cast_parsed_groups(cls, parsed_groups: tuple) -> tuple:
        return tuple(dtype(group) for dtype, group
                     in zip(cls._dtypes, parsed_groups, strict=True))
    
    @classmethod
    def _clean_ext(cls, ext: str) -> str:
        if not cls._exts:
            raise PathError(f"No extensions for {cls}")
        if not ext:
            ext = cls._exts[0]
        elif ext not in cls._exts:
            raise PathError(f"Invalid extension for {cls}: '{ext}'")
        return ext
    
    @classmethod
    def _format_noparse(cls, *args, ext: str) -> str:
        segment = cls._fmt(cls._clean_ext(ext)).format(*args)
        return segment
    
    @classmethod
    def _parse_noformat(cls, segment: Any) -> tuple:
        for pattern in cls._get_patterns():
            if (match := pattern.match(str(segment))):
                fields = cls._cast_parsed_groups(match.groups()[:-1])
                ext = match.groups()[-1]
                return fields, ext
        raise PathError(f"Failed to parse segment '{segment}' for class {cls}")
    
    @classmethod
    def format(cls, *args, ext: str = "") -> str:
        segment = cls._format_noparse(*args, ext=ext)
        assert cls._parse_noformat(segment) == (args, cls._clean_ext(ext))
        return segment
    
    @classmethod
    def parse(cls, segment: Any) -> tuple:
        fields, ext = cls._parse_noformat(segment)
        assert str(segment) == cls._format_noparse(*fields,
                                                   ext=cls._clean_ext(ext))
        return fields

    @classmethod
    def segment_is_valid(cls, segment: Any):
        if not super().segment_is_valid(segment):
            return False
        try:
            cls.parse(segment)
        except PathError:
            return False
        else:
            return True
    
    @property
    def fields(self):
        return self.parse(self._segment)
    
    @property
    def name(self):
        return self._fmt("").format(*self.fields)


class FastaSegment(PatternedSegment):
    _exts = FASTA_EXTS


class FastqSegment(PatternedSegment):
    _exts = FQ_EXTS


class Fastq1Segment(FastqSegment):
    _exts = FQ1_EXTS


class Fastq2Segment(FastqSegment):
    _exts = FQ2_EXTS


class XamSegment(PatternedSegment):
    _exts = XAM_EXTS


class XamIndexSegment(PatternedSegment):
    _exts = XAI_EXTS


class PartitionSegment(SubSegment):
    segment_is_valid = staticmethod(lambda seg: seg in PARTITIONS)


class ModuleSegment(SubSegment):
    segment_is_valid = staticmethod(lambda seg: seg in MODULES)


class TempStepSegment(SubSegment):
    segment_is_valid = staticmethod(lambda seg: seg in TEMP_STEPS)


class SampleSegment(SubSegment):
    pass


class RefSegment(SubSegment):
    pass


class RegionSegment(PatternedSegment):
    _pattern = "([0-9]+)-([0-9]+)"
    _dtypes = (int, int)
    _exts = ("",)
    _fmt = staticmethod(lambda ext: "{}-{}")


class MutVectorReportSegment(PatternedSegment):
    _pattern = "([0-9]+)-([0-9]+)_report"
    _dtypes = (int, int)
    _exts = (".txt",)
    _fmt = staticmethod(lambda ext: "{}-{}_report" + ext)


class MutVectorBatchSegment(PatternedSegment):
    _pattern = "vectors_([0-9]+)"
    _dtypes = (int,)
    _exts = (".orc",)
    _fmt = staticmethod(lambda ext: "vectors_{}" + ext)


# Segment descriptor class

class SegmentDescriptor(object):
    _priv_chr = "_"

    segment_types: Dict[str, type] = {
        "base": BaseSegment,
        "fasta": FastaSegment,
        "fastq": FastqSegment,
        "fastq1": Fastq1Segment,
        "fastq2": Fastq2Segment,
        "partition": PartitionSegment,
        "module": ModuleSegment,
        "temp_step": TempStepSegment,
        "sample": SampleSegment,
        "ref": RefSegment,
        "region": RegionSegment,
        "xam": XamSegment,
        "xai": XamIndexSegment,
        "mv_report": MutVectorReportSegment,
        "mv_batch": MutVectorBatchSegment,
    }

    def __set_name__(self, owner, name: str):
        if name.startswith(self._priv_chr):
            raise ValueError(f"name '{name}' starts with '{self._priv_chr}'")
        self.public_name = name
        self.private_name = f"{self._priv_chr}{name}"
        self.segment_type = self.segment_types[self.public_name]
    
    def __get__(self, obj, objtype=None):
        return getattr(obj, self.private_name)
    
    def __set__(self, obj, value):
        setattr(obj, self.private_name, self.segment_type(value))


# Paths

@dataclass
class BasePath(object):
    base: SegmentDescriptor = SegmentDescriptor()

    @classmethod
    def parse(cls, path: Any):
        if cls is BasePath:
            return cls(path)
        head, tail = os.path.split(str(path))
        parent, = cls.__bases__
        segments = parent.parse(head).args + (tail,)
        return cls(*segments)
    
    @property
    def segments(self):
        return astuple(self)
    
    @property
    def named_segments(self):
        return asdict(self)
    
    @property
    def args(self):
        return tuple(map(str, self.segments))
    
    @property
    def kwargs(self):
        return {name: str(seg) for name, seg in self.named_segments.items()}
    
    @property
    def path(self):
        return Path(*self.args)
    
    @property
    def last(self):
        """ Get the last segment of the path. """
        return astuple(self)[-1]
    
    def __str__(self) -> str:
        return str(self.path)


@dataclass
class PartitionPath(BasePath):
    partition: SegmentDescriptor = SegmentDescriptor()


@dataclass
class ModulePath(PartitionPath):
    module: SegmentDescriptor = SegmentDescriptor()


@dataclass
class SampleOutPath(ModulePath):
    sample: SegmentDescriptor = SegmentDescriptor()


@dataclass
class SampleInPath(BasePath):
    sample: SegmentDescriptor = SegmentDescriptor()


@dataclass
class TempStepPath(ModulePath):
    temp_step: SegmentDescriptor = SegmentDescriptor()


@dataclass
class SampleTempPath(TempStepPath):
    sample: SegmentDescriptor = SegmentDescriptor()


# Alignment Paths


@dataclass
class FastaInPath(BasePath):
    fasta: SegmentDescriptor = SegmentDescriptor()


@dataclass
class FastaSampleOutPath(SampleOutPath):
    fasta: SegmentDescriptor = SegmentDescriptor()


@dataclass
class FastqInPath(BasePath):
    fastq: SegmentDescriptor = SegmentDescriptor()


@dataclass
class Fastq1InPath(BasePath):
    fastq1: SegmentDescriptor = SegmentDescriptor()


@dataclass
class Fastq2InPath(BasePath):
    fastq2: SegmentDescriptor = SegmentDescriptor()


@dataclass
class FastqDemultiPath(SampleOutPath):
    fastq: SegmentDescriptor = SegmentDescriptor()


@dataclass
class Fastq1DemultiPath(SampleOutPath):
    fastq1: SegmentDescriptor = SegmentDescriptor()


@dataclass
class Fastq2DemultiPath(SampleOutPath):
    fastq2: SegmentDescriptor = SegmentDescriptor()


@dataclass
class FastqOutPath(SampleTempPath):
    fastq: SegmentDescriptor = SegmentDescriptor()


@dataclass
class Fastq1OutPath(SampleTempPath):
    fastq1: SegmentDescriptor = SegmentDescriptor()


@dataclass
class Fastq2OutPath(SampleTempPath):
    fastq2: SegmentDescriptor = SegmentDescriptor()


@dataclass
class XamTempPath(SampleTempPath):
    xam: SegmentDescriptor = SegmentDescriptor()


@dataclass
class XamOutPath(SampleOutPath):
    xam: SegmentDescriptor = SegmentDescriptor()


@dataclass
class XamIndexTempPath(SampleTempPath):
    xai: SegmentDescriptor = SegmentDescriptor()


# Vectoring Paths

@dataclass
class XamInPath(SampleInPath):
    xam: SegmentDescriptor = SegmentDescriptor()


@dataclass
class XamIndexInPath(SampleInPath):
    xai: SegmentDescriptor = SegmentDescriptor()


@dataclass
class RefPath(SampleOutPath):
    ref: SegmentDescriptor = SegmentDescriptor()


@dataclass
class RegionPath(RefPath):
    region: SegmentDescriptor = SegmentDescriptor()


@dataclass
class MutVectorReportPath(RefPath):
    mv_report: SegmentDescriptor = SegmentDescriptor()


@dataclass
class MutVectorBatchPath(RegionPath):
    mv_batch: SegmentDescriptor = SegmentDescriptor()
