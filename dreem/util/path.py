"""
Path module of DREEM
================================================================================
Version     0.0.2
Modified    2023-01-30
Author      Matty


Purpose
================================================================================
Several modules in DREEM produce files used by other modules.
For example, the align module creates alignment map files,
from which the vector module generates mutation vector files,
from which the aggregate module computes statistics and generates plots.

A module that writes a file and another module that reads the file must agree on
- the path to the file (so that the second module can find the file)
- the meaning of each part of the path (so that the second module can determine,
  e.g. the names of the sample and the reference from the path)

Although these path conventions could be written separately in each module,
this strategy is not ideal for several reasons:
- It would risk inconsistencies among the modules, potentially causing bugs.
- Changing the conventions would require modifying every module, which would
  not only be tedious but also increase the risk of inconsistencies.
- Defining all of the conventions in one place would reduce the size of the
  code base, improving readability, maintainability, and distribution.  

The path module defines the conventions for all file paths in one central
location, to be referenced by all other modules that use the file system.


Concepts
================================================================================

Paths and path segments
-----------------------
Every path is represented as a collection of path segments.
Each segment is the space between two path separators, except for the first
("top") segment of a path, which can contain separators.

For example, this path to an alignment map file
/home/rfranklin/tmv/exp58/output/alignment/dms2/tmv-genome.bam
would be represented as the following segments:
- top (full path of the experiment's top directory): /home/rfranklin/tmv/exp58/
- partition (finished output or temporary files): output
- module (DREEM module producing the results): alignment
- sample (name of the sample from the experiment): dms2
- ref (name of the reference sequence the sample was aligned to): tmv-genome

Classes of paths
----------------
(Description of what class you would use to represent several different kinds of paths)

Distinction between fields and 


Usage
================================================================================

Creating a path of a specific type using its class
--------------------------------------------------
An instance of a specific type of path can be created in three ways:
1. Calling the class directly, giving the names and values of the path segments
as keyword arguments:
>>> bam_path = OneRefAlignmentOutFilePath(top='/home/rfranklin/tmv/exp58',
...                                    partition='output', module='alignment',
...                                    sample='dms2', ref='tmv-rna', ext='.bam')
>>> print(bam_path)
/home/rfranklin/tmv/exp58/output/alignment/dms2/tmv-genome.bam

2. Calling the class directly, giving the names and values of the path segments
as a dictionary of keyword arguments:
>>> bam_fields = {'top': '/home/rfranklin/tmv/exp58', 'partition': 'output',
...               'module': 'alignment', 'sample': 'dms2', 'ref': 'tmv-rna',
...               'ext': '.bam'}
>>> bam_path = OneRefAlignmentOutFilePath(**bam_fields)
>>> print(bam_path)
/home/rfranklin/tmv/exp58/output/alignment/dms2/tmv-genome.bam

3. Parsing the path from a string (or from any other object whose __str__ method
returns a valid path, such as a pathlib.Path instance):
>>> path = '/home/rfranklin/tmv/exp58/output/alignment/sample/dms2/tmv-rna.bam'
>>> bam_path = OneRefAlignmentOutFilePath.parse_path(path)
>>> print(bam_path.dict())
{'top': '/home/rfranklin/tmv/exp58', 'partition': 'output',
'module': 'alignment', 'sample': 'dms2', 'ref': 'tmv-genome', 'ext': '.bam'}

Creating a path by inferring the type from the fields
-----------------------------------------------------
This module also provides a convenience function called make_path, which
behaves identically to the first two examples (of calling the class directly),
except that instead of the code specifying the type of path to be created,
make_path infers the type from the keywords arguments. For example:

A few notes on the function make_path:
- There is no equivalent of parse_path() for make_path because it is impossible
  to uniquely assign fields to segments of an arbitrary path if neither the
  structure nor the values of the fields are known.
- A TypeError is raised if there is no class of path whose field names match the
  names of the keyword arguments given to make_path.
- In most cases, this function suffices and alleviates the need to import every
  class of path that will be used in a module.


Implementation Details
================================================================================

Representation of paths and path segments
-----------------------------------------
In this module, every path is represented as a sequence of path segments. Every
kind of path is represented by a subclass of BasePath, and every kind of segment
is represented by a subclass of BaseSeg. For example, one subclass of BasePath
represents the full path to a mutation vector report file, and is made of path
segments that represent the name of the report file, the name of the directory
in which it resides, and so on up.

Responsibilities of paths and path segments
-------------------------------------------
Each path segment class is responsible for storing information contained in the
segment (i.e. the value of each field), generating a string representation of
the segment from the values of the fields, and parsing string representations
of the segment to determine the values of the fields encoded within. Each path
class is responsible for storing a series of path segments in the proper order,
ensuring that the path has all expected and no unexpected segments, ensuring
that the values of all of its segments are consistent with each other, and
generating and parsing string representations of the entire path by calling the
corresponding method of each segment and assembling the results.

Path segment classes and inheritance structure
----------------------------------------------
Every path segment class represents a specific component of a path, such as a
directory named after a sample or the name of a mutation vector batch file.

The base class of every path segment is BaseSeg, which does not represent any
specific path segment and should not be instantiated. Every other type of path
segment is a subclass of BaseSeg. BaseSeg is itself decorated with @dataclass
but contains no fields. However, it provides all shared methods for path segment
subclasses (all of which are also decorated with @dataclass):
- Initialization: __init__ (from @dataclass), __post_init__, validate
- Field access: fields, field_names, astuple, asdict, replace, _get_dtypes
- Interface: parse_seg, format_seg

The subclasses of BaseSeg that represent directories are structured as follows:

Subclass of BaseSeg     Fields          Constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
BaseSeg                 -               -
    TopSeg              top             directory must exist
    SubSeg              -               no field may contain a path separator
        PartSeg         partition       must be a valid partition name
        ModSeg          module          must be a valid module name
        TempStepSeg     step            must be a valid step name
        SampleSeg       sample          -
        RefsetSeg       refset          -
        RefSeg          ref             -
        StructSeg       -               -
            RegionSeg   first, last     1 <= first (int) <= last (int)
            ExtenSeg    ext             must be a valid extension for the class
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Most segments for file types inherit from both ExtenSeg and another class that
provides one or more fields that are part of the file name:

Subclass of ExtenSeg    Also subclass of
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MutVectorBatchSeg       -
MutVectorReportSeg	    RegionSeg
XamSeg                  -
    XamMixedSeg         RefsetSeg
    XamSplitSeg         RefSeg
FastaSeg	            -
    FastaMixedSeg	    RefsetSeg
    FastaSplitSeg	    RefSeg
FastqSeg                -
    FastqSampleSeg      SampleSeg
    FastqDemultSeg      RefSeg
Fastq1Seg               -
    Fastq1SampleSeg     SampleSeg
    Fastq1DemultSeg     RefSeg
Fastq2Seg               -
    Fastq2SampleSeg     SampleSeg
    Fastq2DemultSeg     RefSeg
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Implementation of path and path segment classes as Pydantic models
------------------------------------------------------------------


"""


# Imports ######################################################################

from __future__ import annotations
from functools import cache
from inspect import getmembers, isclass, signature
import itertools
import os
from pathlib import Path
import re
from string import ascii_letters, digits
import sys
from typing import Any, ClassVar, Iterable

from pydantic import BaseModel, root_validator, validator


# Constants ####################################################################

# Valid/invalid characters in fields
VALID_CHARS = ascii_letters + digits + "_~=+-"
VALID_FIELD_PATTERN = f"([{VALID_CHARS}]+)"
VALID_CHARS_SET = set(VALID_CHARS)
VALID_CHARS_TOP = set(VALID_CHARS + os.path.sep)
VALID_CHARS_EXT = set(VALID_CHARS + ".")

TOP_KEY = "top"

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
EXT_KEY = "ext"
EXT_PATTERN = "([.].+)"
LOG_EXTS = (".log",)
FASTA_EXTS = (".fasta", ".fa")
FQ_EXTS = (".fastq", ".fq", ".fastq.gz", ".fq.gz")
FQ_PAIRED_EXTS_TEMPLATES = ("_R{}{}", "_mate{}{}", "_{}_sequence{}")
FQ1_EXTS = tuple(template.format(1, ext) for template, ext in
                 itertools.product(FQ_PAIRED_EXTS_TEMPLATES, FQ_EXTS))
FQ2_EXTS = tuple(template.format(2, ext) for template, ext in
                 itertools.product(FQ_PAIRED_EXTS_TEMPLATES, FQ_EXTS))
FQ1_PATTERNS = tuple(VALID_FIELD_PATTERN + f"({ext})" for ext in FQ1_EXTS)
FQ2_PATTERNS = tuple(VALID_FIELD_PATTERN + f"({ext})" for ext in FQ2_EXTS)
SAM_EXT = ".sam"
BAM_EXT = ".bam"
CRAM_EXT = ".cram"
XAM_EXTS = (SAM_EXT, BAM_EXT, CRAM_EXT)
BAI_EXT = f"{BAM_EXT}.bai"
XAI_EXTS = (BAI_EXT,)
XAMI_EXTS = XAM_EXTS + XAI_EXTS

VALID_CHARS_MAP = {TOP_KEY: VALID_CHARS_TOP,
                   EXT_KEY: VALID_CHARS_EXT}


# Path functions ###############################################################

def sanitize(path: Any):
    return os.path.realpath(os.path.normpath(os.path.abspath(str(path))))


# Exceptions ###################################################################

class PathError(Exception):
    """ Base exception class for all errors specific to the Path module. """
    def __init__(self, msg):
        super().__init__(msg)


class FieldCountError(PathError):
    def __init__(self, class_: type[BaseModel], n_fields: int, n_expect: int):
        msg = f"{class_} has {n_fields} fields, but expected {n_expect}"
        super().__init__(msg)


class PathIsNotDirectoryError(PathError):
    def __init__(self, path: Any):
        msg = f"Path does not exist or is not a directory: '{path}'"
        super().__init__(msg)


class FieldContainsPathSeparatorError(PathError):
    def __init__(self, field: str, value: Any):
        msg = f"Field {field}='{value}' contains path separator '{os.path.sep}'"
        super().__init__(msg)


class FieldIsEmptyError(PathError):
    def __init__(self, field: str, value: Any):
        msg = f"Field {field}='{value}' is an empty string"
        super().__init__(msg)


class FieldContainsIllegalCharacterError(PathError):
    def __init__(self, field: str, value: Any):
        msg = f"Field {field}='{value}' contains an illegal character."
        super().__init__(msg)


class InvalidCategoricalFieldError(PathError):
    def __init__(self, value: Any, valid: tuple):
        msg = (f"Value '{value}' was not one of the valid values "
               f"({', '.join(valid)})")
        super().__init__(msg)


class InvalidExtensionError(PathError):
    def __init__(self, ext: Any, exts: tuple):
        msg = f"Invalid extension '{ext}'; expected one of {', '.join(exts)}"
        super().__init__(msg)


class MissingExtensionError(PathError):
    def __init__(self, class_):
        msg = f"{class_} has no extension field ('{EXT_KEY}')"
        super().__init__(msg)


class InvalidNumericalFieldError(PathError):
    def __init__(self, value: Any, criterion: str):
        msg = f"Numeric value '{value}' was not {criterion}"
        super().__init__(msg)


class PathSegmentFormatError(PathError):
    def __init__(self, seg: StructSeg):
        msg = f"Failed to format segment of {seg.__class__}: {seg.dict()}"
        super().__init__(msg)


class PathSegmentParseError(PathError):
    def __init__(self, seg: Any, cls: type[StructSeg]):
        msg = f"Failed to parse segment '{seg}' for {cls} '{cls.pattern_str}'"
        super().__init__(msg)


class PathParseError(PathError):
    def __init__(self, path: Any, cls: type[BasePath]):
        msg = f"Failed to parse path '{path}' for {cls}"
        super().__init__(msg)


class PathOrderError(PathError):
    pass


# Path segment classes #########################################################

class BaseSeg(BaseModel):
    """ Abstract base class for all classes that represent segments of a path.
    Should not be instantiated. """

    @classmethod
    def fields(cls):
        return cls.__fields__

    @classmethod
    def n_fields(cls):
        return len(cls.fields())

    @classmethod
    def keys(cls):
        return list(cls.fields().keys())

    @classmethod
    def build_dict(cls, values: Iterable):
        try:
            return dict(zip(cls.keys(), values, strict=True))
        except ValueError:
            raise FieldCountError(cls, cls.n_fields(), 1)

    def values(self):
        return list(self.dict().values())

    @classmethod
    def parse_seg(cls, segment: Any):
        """ Return a new instance of the class given a segment as a str-like
        argument. Non-structured segments (including this base class) must have
        exactly one field. """
        return cls(**cls.build_dict([segment]))
    
    def format_seg(self):
        """ Return a string representation of the first field of the segment.
        Non-structured segments (including this base class) must have exactly
        one field, else a NotImplementedError will be raised. """
        try:
            value, = self.values()
        except ValueError:
            raise FieldCountError(self.__class__, self.n_fields(), 1)
        return str(value)

    def __str__(self):
        """ Alias for ```self.format()``` """
        return self.format_seg()


class TopSeg(BaseSeg):
    """ Class representing the top-level working directory of DREEM. All
    temporary and final output files will be located in this directory. This
    directory must exist at the time DREEM is launched. """
    top: str

    @validator(TOP_KEY)
    def top_dir_must_exist(cls, top):
        if not os.path.isdir(top):
            raise NotADirectoryError(top)
        return top


class SubSeg(BaseSeg):
    @root_validator(pre=True)
    def fields_must_not_be_empty(cls, values: dict[str, Any]):
        """ Validate that no field except 'ext' has an empty value. """
        for key, value in values.items():
            if key != EXT_KEY and not str(value):
                raise FieldIsEmptyError(key, value)
        return values

    @root_validator(pre=True)
    def fields_must_not_contain_path_separator(cls, values: dict[str, Any]):
        """ Validate that no field except 'top' contains a path separator. """
        for key, value in values.items():
            if key != TOP_KEY and os.path.sep in str(value):
                raise FieldContainsPathSeparatorError(key, value)
        return values

    @root_validator(pre=True)
    def fields_must_contain_only_valid_characters(cls, values: dict[str, Any]):
        """ Validate that all fields contain no invalid field characters. """
        for key, value in values.items():
            if set(value) - VALID_CHARS_MAP.get(key, VALID_CHARS_SET):
                raise FieldContainsIllegalCharacterError(key, value)
        return values


class SampleSeg(SubSeg):
    """ Segment for a directory named after a sample. """
    sample: str


class RefsetSeg(SubSeg):
    """ Segment for a directory named after a set of reference sequences. """
    refset: str


class OneRefSeg(SubSeg):
    """ Segment for a directory named after one reference sequence. """
    ref: str


class PartSeg(SubSeg):
    """ Segment for a partition directory ('output' or 'temp'). """
    partition: str

    @validator("partition")
    def partition_is_valid(cls, partition):
        if partition not in PARTITIONS:
            raise InvalidCategoricalFieldError(partition, PARTITIONS)
        return partition


class ModSeg(SubSeg):
    """ Segment for a directory containing the outputs of a module. """
    module: str

    @validator("module")
    def module_is_valid(cls, module):
        if module not in MODULES:
            raise InvalidCategoricalFieldError(module, MODULES)
        return module


class TempStepSeg(SubSeg):
    """ Segment for a directory containing the temporary files of a step. """
    step: str

    @validator("step")
    def step_is_valid(cls, step):
        if step not in TEMP_STEPS:
            raise InvalidCategoricalFieldError(step, TEMP_STEPS)
        return step


class BatchSeg(SubSeg):
    """ Segment for a batch. """
    batch: int

    @validator("batch")
    def batch_is_not_negative(cls, batch):
        if batch < 0:
            raise InvalidNumericalFieldError(batch, ">= 0")
        return batch


class StructSeg(SubSeg):
    """
    Segment representing any path segment whose name contains an internal
    structure with one or more fields. For example, a segment could encode the
    name of a sample and the year, month, and day on which it was generated.
    This base class provides all the methods of structured segments but does
    not itself have any fields.

    Class attributes
    ----------------
    format_str: str
        Specify how to represent the path segment as a string, using Python's
        string formatting mini-language. The string is created by calling
        ```self.format_str.format(*self.dict())```.

        Note: ```format_str``` must be compatible with ```pattern_str```
    pattern_str: str
        Specify how to parse a string representing the path to determine the
        values of the fields. The string is parsed using ```re.match```,
        and ```match.groups()``` is used to get the values of the fields.

        Note: ```pattern_str``` must be compatible with ```format_str```

    Examples
    --------
    >>> class SampleBatchSeg(StructSeg, BatchSeg, SampleSeg):
    ...     format_str = "s-{}_b-{}"
    ...     pattern_str = f"s-{VALID_FIELD_PATTERN}_b-([0-9]+)"
    >>> seg = SampleBatchSeg(sample="tmv", batch=37)
    >>> assert seg.format_seg() == 's-tmv_b-37'
    """

    format_str: ClassVar[str] = ""
    pattern_str: ClassVar[str] = ""

    def _format_seg(self):
        """ Return the result of calling the ```format``` method of
        ```format_str``` with the field values as arguments. """
        return self.format_str.format(*self.values())

    def format_seg(self) -> str:
        """
        Return a string representation of the path segment by inserting the
        values of the segment's fields into its format specification.

        Returns
        -------
        str
            A representation of the path segment as a string.

        Raises
        ------
        PathSegmentFormatError
            if the resulting string cannot be parsed into a new instance whose
            fields match those of the current instance.
        """
        # Create the formatted string.
        segstr = self._format_seg()
        try:
            # Try to parse that string into a new instance, then confirm that
            # the fields of the new instance match those of the self instance.
            if self._parse_seg(segstr) == self:
                # If the fields match, then the formatted string is valid.
                return segstr
        except PathSegmentParseError:
            pass
        # If the parsing raised an error or the fields did not match, then an
        # error occurred during formatting.
        raise PathSegmentFormatError(self)

    @classmethod
    def _parse_seg(cls, segstr: Any):
        """ Return a new instance of the class by parsing a ```str```-like
        ```segstr``` with ```re.match``` (using ```pattern_str``` as the
        regular expression pattern), then passing the parsed groups to the
        class constructor. """
        # Try to parse segstr with re.match, using pattern_str as the pattern.
        if match := re.match(cls.pattern_str, str(segstr)):
            # If segstr matched the pattern, then create keyword arguments from
            # the matched values, and use them to initialize a new instance.
            return cls(**cls.build_dict(match.groups()))
        # If segstr did not match the pattern, then the parsing failed.
        raise PathSegmentParseError(segstr, cls)

    @classmethod
    def parse_seg(cls, segstr: Any):
        """
        Return a new instance of the class by parsing the values of the fields
        from a string and passing them to the class constructor.

        Parameters
        ----------
        segstr: Any
            Literal representation of the segment to parse as a ```str```
            (or as any other object ```x``` for which ```str(x)``` is a valid
            path segment, e.g. an instance of ```pathlib.Path```).
        
        Returns
        -------
        seginst: StructSeg
            New instance of a structured path segment with fields parsed from
            ```segstr```. Satisfies ```str(seginst) == str(segstr)```.

        Raises
        ------
        PathSegmentParseError
            if the string representation of the newly created instance of the
            class did not match ```segstr```.
        """
        # Create a new instance of the class by parsing segstr.
        seginst = cls._parse_seg(segstr)
        # Confirm that formatting the new instance yields the original segstr.
        if seginst._format_seg() != str(segstr):
            # If not, an error occurred during parsing.
            raise PathSegmentParseError(segstr, cls)
        # If so, the new instance was parsed correctly.
        return seginst


class RegionSeg(StructSeg):
    """
    Segment for a directory of a region of a reference sequence.

    Fields
    ------
    first: int
        The first coordinate in the region; 1-indexed, inclusive.
    last: int
        The last coordinate in the region; 1-indexed, inclusive.
    """
    first: int
    last: int

    format_str = "{}-{}"
    pattern_str = "([0-9]+)-([0-9]+)"

    @root_validator()
    def first_last_are_ordered(cls, values):
        """ Validate that 1 <= first <= last """
        first, last = map(int, (values["first"], values["last"]))
        if first < 1:
            raise InvalidNumericalFieldError(first, ">= 1")
        if last < first:
            raise InvalidNumericalFieldError(last, f">= first ({first})")
        return values


class FileSeg(StructSeg):
    """
    Segment representing a file name with an extension.

    Fields
    ------
    ext: str
        The file extension. Must begin with a period.

    Class attributes
    ----------------
    exts: tuple[str]
        All valid file extensions for the class. Each must begin with a period.
    """
    ext: str
    exts: ClassVar[tuple[str, ...]] = tuple()

    format_str = "{}{}"
    pattern_str = VALID_FIELD_PATTERN + EXT_PATTERN

    @validator(EXT_KEY)
    def valid_file_extension(cls, ext):
        """ Validate the file extension (```ext```). It must be an element of
        the class attribute ```exts```. """
        if ext not in cls.exts:
            raise InvalidExtensionError(ext, cls.exts)
        return ext

    @classmethod
    def keys_no_ext(cls):
        """ Return a list of field keys without the 'ext' field. """
        # Create a new list of the field keys.
        keys = cls.keys()
        try:
            # Remove 'ext' from the list.
            keys.remove(EXT_KEY)
        except ValueError:
            # This should not happen, but in case 'ext' was not in keys,
            # raise an error.
            raise MissingExtensionError(cls)
        # Return the list of keys without 'ext'.
        return keys

    @classmethod
    def _parse_seg(cls, segstr: Any):
        if all(ext[0] not in VALID_CHARS_SET for ext in cls.exts):
            # If the extension starts with a character that is not valid if it
            # occurs in a field (e.g. the most common starting character '.'),
            # then the extension can be distinguished from the other fields by
            # the regular expression matching defined in the superclass.
            return super()._parse_seg(segstr)
        # If the extension starts with a character that is valid if it occurs
        # in a field (e.g. '_R1.fq' is such an extension because it starts
        # with '_', which may occur in a field), then we need to match the
        # segment minus the extension to the pattern minus its extension.
        segstr_ = str(segstr)
        # Remove the final group, denoted by (), from the pattern string. This
        # group captures the extension, so it is removed.
        pattern_trunc = cls.pattern_str[:cls.pattern_str.rindex("(")]
        for ext in cls.exts:
            # Find the extension in exts that segstr ends with.
            if segstr_.endswith(ext):
                # Remove that extension from segstr.
                segstr_trunc = segstr_[:-len(ext)]
                assert segstr_ == segstr_trunc + ext
                # See if the truncated segment matches the truncated pattern.
                if match := re.match(pattern_trunc, segstr_trunc):
                    # If the truncated strings matched, then add the extension
                    # back to the end of the matched groups; we already know it
                    # matched because segstr_.endswith(ext) returned True.
                    return cls(**cls.build_dict(match.groups() + (ext,)))
        # If segstr did not match any extension, or did match an extension but
        # failed to match any patterns after removing the extension, then the
        # parsing failed.
        raise PathSegmentParseError(segstr, cls)


class MutVectorReportFileSeg(FileSeg, RegionSeg):
    """ Segment for a mutation vector report file. """
    format_str = "{}-{}_report{}"
    pattern_str = f"([0-9]+)-([0-9]+)_report{EXT_PATTERN}"
    exts = (".txt",)


class MutVectorBatchFileSeg(FileSeg, BatchSeg):
    """ Segment for a mutation vector batch file. """
    format_str = "vectors_{}{}"
    pattern_str = "vectors_([0-9]+)" + EXT_PATTERN
    exts = (".orc",)


class AbstractRefFileSeg(FileSeg):
    """ Abstract segment for a FASTA file. Should not be instantiated. """
    exts = FASTA_EXTS


class RefsetSeqFileSeg(AbstractRefFileSeg, RefsetSeg):
    """ Segment for a FASTA file that can contain one or more references. The
    refset field is the name of the file, minus the file extension. """


class OneRefFileSeg(AbstractRefFileSeg, OneRefSeg):
    """ Segment for a FASTA file that contains exactly one reference sequence.
    The ref field is the reference name (on the first line, after '>'). """


class AbstractReadsFileSeg(FileSeg):
    """ Abstract segment for a standalone (single-end or paired-interleaved)
    FASTQ file. Should not be instantiated. """
    exts = FQ_EXTS


class SampleReadsFileSeg(AbstractReadsFileSeg, SampleSeg):
    """ Segment for a standalone FASTQ file containing all reads from a
    sample. """


class DemultReadsFileSeg(AbstractReadsFileSeg, OneRefSeg):
    """ Segment for a standalone FASTQ file containing all reads from one
    construct (i.e. reference sequence) in one sample. """


class AbstractReads1FileSeg(FileSeg):
    """ Abstract segment for a paired-end FASTQ file containing mate-1 reads.
    Should not be instantiated. """
    exts = FQ1_EXTS


class SampleReads1FileSeg(AbstractReads1FileSeg, SampleSeg):
    """ Segment for a paired-end FASTQ file containing all mate-1 reads from a
    sample. """


class DemultReads1FileSeg(AbstractReads1FileSeg, OneRefSeg):
    """ Segment for a paired-end FASTQ file containing all mate-1 reads from one
    construct (i.e. reference sequence) in one sample. """


class AbstractReads2FileSeg(FileSeg):
    """ Abstract segment for a paired-end FASTQ file containing mate-2 reads.
    Should not be instantiated. """
    exts = FQ2_EXTS


class SampleReads2FileSeg(AbstractReads2FileSeg, SampleSeg):
    """ Segment for a paired-end FASTQ file containing all mate-2 reads from a
    sample. """


class DemultReads2FileSeg(AbstractReads2FileSeg, OneRefSeg):
    """ Segment for a paired-end FASTQ file containing all mate-2 reads from one
    construct (i.e. reference sequence) in one sample. """


class AbstractAlignmentFileSeg(FileSeg):
    """ Segment representing both alignment map files (.sam, .bam, .cram)
    and index files (.bam.bai) """
    exts = XAMI_EXTS


class RefsetAlignmentFileSeg(AbstractAlignmentFileSeg, RefsetSeg):
    """ Segment representing an ailgnment map file resulting from aligning a
    set of reads (i.e. FASTQ) to a set of references (i.e. FASTA). """


class OneRefAlignmentFileSeg(AbstractAlignmentFileSeg, OneRefSeg):
    """ Segment representing an alignment map file produced by splitting another
    alignment map file into one file for each reference. """


# Path classes #################################################################


class BasePath(BaseModel):
    """ Abstract base class for all classes that represent full, absolute paths.
    Should not be instantiated. """

    @classmethod
    def _segment_types(cls):
        """ Return a list of the type of every segment in the path. """
        seg_types: list[type[BaseSeg]] = list()
        # Search the base classes for path segments in reverse order.
        for base in reversed(cls.__bases__):
            if issubclass(base, BasePath):
                # The base class represents a path and thus also has base
                # classes representing path segments: add them to the list.
                seg_types.extend(base._segment_types())
            elif issubclass(base, BaseSeg):
                # The base class is a path segment but not a path or path
                # scaffold class. Thus, it should be in the list of segments.
                seg_types.append(base)
            else:
                # Otherwise, the base class neither is a nor contains any
                # path segments. Skip it.
                pass
        return seg_types

    @classmethod
    @cache
    def segment_types(cls) -> tuple[type[BaseSeg], ...]:
        """ Return a tuple of the type of every segment in the path, in order
        from the beginning to the end of the path. """
        seg_types = tuple(cls._segment_types())
        if seg_types and seg_types[0] is not TopSeg:
            raise PathOrderError(f"First segment of {cls} is not {TopSeg}")
        return seg_types

    def segments(self) -> tuple[BaseSeg, ...]:
        """ Return a tuple of an instance of every segment in the path, in order
        from the beginning to the end of the path. """
        return tuple(seg_type(**{key: getattr(self, key)
                                 for key in seg_type.keys()})
                     for seg_type in self.segment_types())

    @classmethod
    def _parse_fields(cls, seg_types: list[type[BaseSeg]], path: str):
        """
        Return a dict of the names and values of the fields encoded within a
        string representation of a path.

        Parameters
        ----------
        seg_types: list[type[BaseSeg]]
            A list of the types of segments in the path, in order from beginning
            to end of the path. The first item must be ```TopSeg```.
        path: str
            The path to parse, represented as a string.

        Returns
        -------
        dict[str, Any]
            A dict of the names (keys) and values (values) of the fields encoded
            in the path string.

        Raises
        ------
        PathParseError
            if any part of the path remains to be parsed after all the types of
            segments have been used.
        """
        if seg_types:
            # If at least one more segment type needs to be parsed, then get the
            # type that needs to be parsed now (from the end of seg_types).
            seg_type = seg_types.pop()
            if seg_types:
                # If any segment types need to be parsed after the current one
                # has been parsed, then parse the last segment of the path now
                # (with the current segment type), and parse everything above
                # it during the next recursive call (with the segment types
                # that need to be parsed after the current one).
                parse_next, parse_now = os.path.split(path)
            else:
                # If the current segment type is the last one that needs to be
                # parsed, then parse the entire remaining path with the current
                # segment type, and there is no part of the path to parse next.
                parse_next, parse_now = "", path
            # Parse the current segment with the current segment type, and the
            # next segment(s) with the remaining segment type(s), then merge.
            return {**cls._parse_fields(seg_types, parse_next),
                    **seg_type.parse_seg(parse_now).dict()}
        else:
            # No segment types still need to be parsed.
            if path:
                # Any part of the path that has not yet been parsed cannot be,
                # since no segment types are left to parse it.
                raise PathParseError(path, cls)
            # Return a dict with no fields, signifying that nothing remains to
            # be parsed.
            return {}

    @classmethod
    def parse_path(cls, path: Any):
        """
        Return a new instance of the class by parsing a given path.

        Parameters
        ----------
        path: Any
            A ```str``` representing the path to be parsed (or any other
            object ```x``` for which ```str(x)``` returns a valid path,
            e.g. an instance of ```pathlib.Path```).

        Returns
        -------
        BasePath
            A new instance of this class that represents the given path.

        Raises
        ------
        PathParseError
            if the newly created instance of the class yields a string
            representation that does not match the ```path``` argument.
        """
        pathstr = str(path)
        pathinst = cls(**cls._parse_fields(list(cls.segment_types()), pathstr))
        if str(pathinst) != pathstr:
            raise PathParseError(path, cls)
        return pathinst

    @property
    def path(self):
        """ Return a ```pathlib.Path``` instance representing this path. """
        return Path(*map(str, self.segments()))

    def format_path(self):
        """ Return a ```str``` representing this path. """
        return str(self.path)

    def __str__(self):
        """ Alias for ```self.format_path()``` """
        return self.format_path()

    @root_validator(pre=False)
    def override_str(cls, values):
        """
        Copy the ```__str__``` method of ```BasePath``` to each subclass. This
        step is necessary because all path classes have ```BasePath``` as their
        lowest-level base class and then inherit from various subclasses of
        ```BaseSeg``` on top. Thus, any method of ```BasePath``` that has the
        same name as a method of ```BaseSeg``` will be overridden in the derived
        path class. By design, no methods share names if such sharing can be
        avoided (hence methods that perform a similar function in path and
        segment classes, such as ```parse_path``` and ```parse_seg```, are given
        different names). One exception is ```__str__```: since it cannot be
        renamed, it must be copied from ```BasePath``` to each subclass.
        """
        method_names = ("__str__",)
        for method_name in method_names:
            setattr(cls, method_name, getattr(BasePath, method_name))
        return values


# General directory paths

class TopDirPath(TopSeg, BasePath):
    def upto(self, field_name: str, include: bool = True):
        """
        Return a new instance of a path class representing the path up to the
        field ```field_name```.

        Parameters
        ----------
        field_name: str
            Name of the field up to which the new path should go.
        include: bool
            Whether to include the field ```field_name``` in the new path
            (default: ```True```).

        """
        keys_upto = (ks := self.keys())[:(n := ks.index(field_name) + include)]
        values_upto = self.values()[:n]
        return assemble(**dict(zip(keys_upto, values_upto)))
    
    def replace(self, **changes):
        """
        Return a new instance of a path class with some fields changed.

        Parameters
        ----------
        changes
            Keyword arguments giving the changes to make to the path. Each key
            is a field name and each value is the new value of the field.
        
        Returns
        -------
        TopDirPath
            A new instance of the same class with fields replaced by the values
            in ```changes```.
        """
        return self.__class__(**{**self.dict(), **changes})


class PartitionDirPath(PartSeg, TopDirPath):
    pass


class ModuleDirPath(ModSeg, PartitionDirPath):
    pass


class TempDirPath(TempStepSeg, ModuleDirPath):
    pass


# Sample, Reference, and Region directory paths

class SampleInDirPath(SampleSeg, TopDirPath):
    pass


class SampleTempDirPath(SampleSeg, TempDirPath):
    pass


class SampleOutDirPath(SampleSeg, ModuleDirPath):
    pass


class RefInDirPath(OneRefSeg, SampleInDirPath):
    pass


class RefTempDirPath(OneRefSeg, SampleTempDirPath):
    pass


class RefOutDirPath(OneRefSeg, SampleOutDirPath):
    pass


class RegionInDirPath(RegionSeg, RefInDirPath):
    pass


class RegionTempDirPath(RegionSeg, RefTempDirPath):
    pass


class RegionOutDirPath(RegionSeg, RefOutDirPath):
    pass


# Reference sequence (FASTA) file paths

class RefsetSeqInFilePath(RefsetSeqFileSeg, TopDirPath):
    pass


class OneRefSeqTempFilePath(OneRefFileSeg, TempDirPath):
    pass


# Read (FASTQ) file paths

class SampleReadsInFilePath(SampleReadsFileSeg, TopDirPath):
    pass


class SampleReads1InFilePath(SampleReads1FileSeg, TopDirPath):
    pass


class SampleReads2InFilePath(SampleReads2FileSeg, TopDirPath):
    pass


class SampleReadsTempFilePath(SampleReadsFileSeg, TempDirPath):
    pass


class SampleReads1TempFilePath(SampleReads1FileSeg, TempDirPath):
    pass


class SampleReads2TempFilePath(SampleReads2FileSeg, TempDirPath):
    pass


class DemultReadsInFilePath(DemultReadsFileSeg, SampleInDirPath):
    pass


class DemultReads1InFilePath(DemultReads1FileSeg, SampleInDirPath):
    pass


class DemultReads2InFilePath(DemultReads2FileSeg, SampleInDirPath):
    pass


class DemultReadsTempFilePath(DemultReadsFileSeg, SampleTempDirPath):
    pass


class DemultReads1TempFilePath(DemultReads1FileSeg, SampleTempDirPath):
    pass


class DemultReads2TempFilePath(DemultReads2FileSeg, SampleTempDirPath):
    pass


class DemultReadsOutFilePath(DemultReadsFileSeg, SampleOutDirPath):
    pass


class DemultReads1OutFilePath(DemultReads1FileSeg, SampleOutDirPath):
    pass


class DemultReads2OutFilePath(DemultReads2FileSeg, SampleOutDirPath):
    pass


# Alignment (SAM/BAM/CRAM) file paths

class RefsetAlignmentInFilePath(RefsetAlignmentFileSeg, SampleInDirPath):
    pass


class RefsetAlignmentTempFilePath(RefsetAlignmentFileSeg, SampleTempDirPath):
    pass


class RefsetAlignmentOutFilePath(RefsetAlignmentFileSeg, SampleTempDirPath):
    pass


class OneRefAlignmentInFilePath(OneRefAlignmentFileSeg, SampleInDirPath):
    pass


class OneRefAlignmentTempFilePath(OneRefAlignmentFileSeg, SampleTempDirPath):
    pass


class OneRefAlignmentOutFilePath(OneRefAlignmentFileSeg, SampleOutDirPath):
    pass


# Vectoring file paths

class MutVectorBatchFilePath(MutVectorBatchFileSeg, RegionOutDirPath):
    pass


class MutVectorReportFilePath(MutVectorReportFileSeg, RefOutDirPath):
    pass


# Path managing functions ######################################################

def is_path_class(item: type):
    """
    Return whether ```item``` is a class of path that can be instantiated. It
    must be a subclass of both ```BasePath``` (to provide the methods for
    handling paths) and ```BaseSeg``` (to provide the segment(s) of the path,
    since every path that can be instantiated contains at least one segment).
    """
    return (isclass(item)
            and issubclass(item, BasePath)
            and issubclass(item, BaseSeg))


@cache
def _get_path_classes() -> dict[str, type[TopDirPath]]:
    """
    Return a ```dict``` of every valid, instantiable class of path defined in
    this module. Its keys are the class names (type ```str```) and its values
    the class objects (type ```type```).
    """
    return dict(getmembers(sys.modules[__name__], is_path_class))


def _sorted_fields(fields: Iterable[str]):
    """ Given an iterable of field names, return a tuple of the names in
    alphabetical order. """
    return tuple(sorted(fields))


@cache
def _get_path_signatures():
    return {name: _sorted_fields(signature(cls).parameters)
            for name, cls in _get_path_classes().items()}


def _path_class_from_signature(**fields):
    names = _sorted_fields(fields)
    matches = [_get_path_classes()[name]
               for name, sig in _get_path_signatures().items()
               if sig == names]
    if matches:
        if len(matches) > 1:
            raise TypeError(f"Multiple matching signatures for fields {names}")
        return matches[0]
    raise TypeError(f"No matching signature for fields: {names}")


def assemble(**fields):
    return _path_class_from_signature(**fields)(**fields)


# Path converters ##############################################################

class PathTypeMapper(object):
    _mapping: dict[type[TopDirPath], type[TopDirPath]] = dict()
    _inverse: type[PathTypeMapper] | None = None

    @classmethod
    def items(cls):
        return cls._mapping.items()

    @classmethod
    def map_type(cls, input_type: type[TopDirPath]):
        """ Return the type to which a given type maps in the dictionary that
        this class defines. A thin wrapper around the dictionary's get-item
        method so that the dictionary does not need to be exposed. """
        try:
            return cls._mapping[input_type]
        except KeyError:
            raise ValueError(f"{input_type} is not an input type for {cls}")

    @classmethod
    def map_inst(cls, orig_inst: TopDirPath, preserve_type: bool = False,
                 **new_fields):
        """
        Cast an instance of a path of a given type to the output type to which
        it corresponds in this PathTypeMapper class. Optionally, add or modify
        the fields of the new type.

        Parameters
        ----------
        orig_inst: TopDirPath
            Instance of a path to be cast into a new type. The new type is
            determined by mapping the type of orig_inst via the dictionary of
            types that this class of PathTypeMapper defines.
        preserve_type: bool = False
            Whether to return an instance of the same type as ```orig_inst```.
            If False (the default), the returned type will be the one to which
            the input type maps. If True, the output is converted back to the
            input type without changing the path that is represented by the
            output. Note that the output path cannot always be represented by
            the input type; if it cannot, a parsing error will be raised.
        **new_fields
            Keyword arguments of fields to add to or modify in the new instance.
            Any fields of the new type that are not in the original type must be
            given, else initializing the new instance will fail. To catch any
            mistyped field names, an error will be raised if any fields are not
            in the new type.

        Returns
        -------
        TopDirPath
            Instance of a path in the new type, possibly with new fields.

        Raises
        ------
        TypeError
            if any key of **new_fields is not a valid field in the output type.
        """
        # Determine the type of the instance to return.
        orig_type = type(orig_inst)
        new_type = cls.map_type(orig_type)
        # Check if any given fields are not defined in the output type.
        new_keys = set(new_type.keys())
        if extras := set(new_fields) - new_keys:
            # If so, then in the best-case scenario, nothing happens; in the
            # worst-case scenario, a keyword was mistyped, potentially causing
            # a subtle bug. Raise an error to protect against the latter case.
            raise TypeError(f"Got keys not in {new_type}: {', '.join(extras)}")
        # Remove any fields that were defined in the input type but not in the
        # output type. This is not an error: any extra fields input type will
        # just get lumped into the catch-all 'top' field in the output type.
        orig_fields = {key: value for key, value in orig_inst.dict().items()
                       if key in new_keys}
        # Create the new instance by combining the original fields (minus any
        # that are not defined in the output type) with the new fields (if a
        # key appears in both original and new fields, the new value overrides
        # the original) and passing to the constructor of the new type.
        new_inst = new_type(**{**orig_fields, **new_fields})
        if preserve_type:
            # If the new path is to be cast back to the original data type,
            # then parse the new path using the original type. Note that it
            # could be impossible to parse the new path with the original
            # type, in which case an error will be raised during parsing.
            return orig_type.parse_path(new_inst.path)
        else:
            # Otherwise, return the new instance of the new type.
            return new_inst

    @classmethod
    def inverse(cls):
        """
        Return a new subclass of PathTypeMapper that maps the output types to
        the input types. Raise a TypeError if the mapping is not invertible.
        """
        if cls._inverse is None:
            mapping = dict()
            for type_in, type_out in cls.items():
                if type_out in mapping:
                    raise TypeError(f"Mapping of {cls} is not invertible")
                mapping[type_out] = type_in

            class InverseMapper(PathTypeMapper):
                _mapping = mapping
                _inverse = cls

            InverseMapper.__name__ = f"{cls.__name__}_InverseMapper"
            cls._inverse = InverseMapper

        return cls._inverse

    @staticmethod
    def chain(name: str,
              *mappers: type[PathTypeMapper]):
        """
        Return a new subclass of PathTypeMapper that chains together the given
        mappers, passing the output of one mapper into the next, and so on.
        """
        mapping = dict()
        if len(mappers) < 2:
            raise TypeError("At least two mappers must be given.")
        for type_in, type_out in mappers[0].items():
            for mapper in mappers[1:]:
                type_out = mapper.map_type(type_out)
            mapping[type_in] = type_out

        class MapInToOut(PathTypeMapper):
            _mapping = mapping

        MapInToOut.__name__ = name
        return MapInToOut

    @classmethod
    def compose(cls, other: type[PathTypeMapper], name: str):
        return cls.chain(name, cls, other)


class ReadsInToReadsTemp(PathTypeMapper):
    _mapping = {
        SampleReadsInFilePath: SampleReadsTempFilePath,
        SampleReads1InFilePath: SampleReads1TempFilePath,
        SampleReads2InFilePath: SampleReads2TempFilePath,
        DemultReadsInFilePath: DemultReadsTempFilePath,
        DemultReads1InFilePath: DemultReads1TempFilePath,
        DemultReads2InFilePath: DemultReads2TempFilePath,
    }


class ReadsInToAlignmentTemp(PathTypeMapper):
    _mapping = {
        SampleReadsInFilePath: RefsetAlignmentTempFilePath,
        SampleReads1InFilePath: RefsetAlignmentTempFilePath,
        SampleReads2InFilePath: RefsetAlignmentTempFilePath,
        DemultReadsInFilePath: OneRefAlignmentTempFilePath,
        DemultReads1InFilePath: OneRefAlignmentTempFilePath,
        DemultReads2InFilePath: OneRefAlignmentTempFilePath,
    }


class AlignmentInToAlignmentTemp(PathTypeMapper):
    _mapping = {
        OneRefAlignmentInFilePath: OneRefAlignmentTempFilePath,
        RefsetAlignmentInFilePath: RefsetAlignmentTempFilePath,
    }


class AlignmentTempToAlignmentOut(PathTypeMapper):
    _mapping = {
        OneRefAlignmentTempFilePath: OneRefAlignmentOutFilePath,
        RefsetAlignmentTempFilePath: RefsetAlignmentOutFilePath,
    }


AlignmentInToAlignmentOut = PathTypeMapper.chain("AlignmentInToAlignmentOut",
                                                 AlignmentInToAlignmentTemp,
                                                 AlignmentTempToAlignmentOut)
