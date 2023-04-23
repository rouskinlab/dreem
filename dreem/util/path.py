"""
Path module of DREEM
========================================================================
Purpose
========================================================================
Several modules in DREEM produce files used by other modules.
For example, the alignment module creates alignment map files,
from which the vectoring module generates mutation vector files,
from which the aggregate module computes statistics and generates plots.
Modules that pass files to each other must agree on
- the path to the file, so that the second module can find the file
- the meaning of each part of the path, so that the second module can
  parse information contained in the path
Although these path conventions could be written separately in each
module this strategy is not ideal for several reasons:
- It would risk inconsistencies among the modules, causing bugs.
- Changing the conventions would require modifying every module,
  which would be not only tedious but also risky for the first reason.
- Defining all the conventions in one place would reduce the size of the
  code base, improving readability, maintainability, and distribution.  
This module defines all file path conventions for all other modules.
Concepts
========================================================================
Paths and path segments
-----------------------
Every path is represented as a collection of path segments.
Each segment is the space between two path separators, except for the
first ("top") segment of a path, which can contain separators.
For example, this path to an alignment map file
```/home/rfranklin/tmv/exp58/output/alignment/dms2/tmv-genome.bam```
would be represented as the following segments:
- top (full path of the top-level directory): ```/home/rfranklin/tmv/exp58/```
- module (DREEM module producing the results): ```alignment```
- sample (name of the sample from the experiment): ```dms2```
- ref (name of the reference sequence the sample was aligned to): tmv-genome
Classes of paths
----------------
Usage
========================================================================
Creating a path of a specific type using its class
--------------------------------------------------
An instance of a specific type of path can be created in three ways:
1. Calling the class directly, giving the names and values of the path
segments as keyword arguments:
>>> xam_inp = OneRefAlignmentOutFilePath(top=os.getcwd(),
...                                       module=Module.ALIGN,
...                                       sample="dms2",
...                                       ref="tmv-rna",
...                                       ext=".bam")
>>> assert str(xam_inp) == os.getcwd() + "/alignment/dms2/tmv-rna.bam"
2. Calling the class directly, giving the names and values of the path
segments as a dictionary of keyword arguments:
>>> bam_fields = {"top": os.getcwd(), "module": Module.ALIGN,
...               "sample": "dms2", "ref": "tmv-rna", "ext": ".bam"}
>>> xam_inp = OneRefAlignmentOutFilePath(**bam_fields)
>>> assert str(xam_inp) == os.getcwd() + "/alignment/dms2/tmv-rna.bam"
3. Parsing the path from a string (or from any other object whose
   __str__ method returns a valid path, such as a pathlib.Path object):
>>> path = os.getcwd() + "/alignment/dms2/tmv-rna.bam"
>>> xam_inp = OneRefAlignmentOutFilePath.parse(path)
>>> assert xam_inp.dict() == {"top": os.getcwd(),
...                            "module": Module.ALIGN,
...                            "sample": "dms2",
...                            "ref": "tmv-rna",
...                            "ext": ".bam"}
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
========================================================================
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
Subclass of BaseSeg     Fields      Constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
BaseSeg                 -           -
    TopSeg              top         must have permission to read & write
    SubSeg              -           may not contain a path separator
        ModSeg          module      must be a valid module name
        StepStepSeg     step        must be a valid step name
        SampleSeg       sample      -
        RefsetSeg       refset      -
        RefSeg          ref         -
        StructSeg       -           -
            SectionSeg   end5, end3  1 ≤ end5 (int) ≤ end3 (int)
            FileSeg     ext         must be a valid file extension
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Most segments for file types inherit from both ExtenSeg and another class that
provides one or more fields that are part of the file name:
Subclass of ExtenSeg    Also subclass of
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
MutVectorBatchSeg       -
MutVectorReportSeg	    SectionSeg
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

# Imports ##############################################################

from __future__ import annotations
from enum import Enum
from functools import cache
from inspect import getmembers, isclass
import itertools
import os
import pathlib
import re
from string import ascii_letters, digits
import sys
from typing import Any, ClassVar, Iterable

from pydantic import BaseModel, Extra, NonNegativeInt, PositiveInt
from pydantic import root_validator, validator


# Constants ############################################################

# Valid/invalid characters in fields
VALID_CHARS = ascii_letters + digits + "_~=+-"
VALID_CHARS_SET = set(VALID_CHARS)
VALID_FIELD_PATTERN = f"([{VALID_CHARS}]+)"
VALID_FIELD_REGEX = re.compile(f"^{VALID_FIELD_PATTERN}$")

TOP_KEY = "top"
EXT_KEY = "ext"
EXTS_KEY = "exts"


class Module(Enum):
    DEMULT = "demultiplexing"
    ALIGN = "alignment"
    VECTOR = "vectoring"
    CLUSTER = "clustering"
    AGGREG = "aggregation"


class Step(Enum):
    ALIGN_REFS = "align-0_refs"
    ALIGN_TRIM = "align-1_trim"
    ALIGN_ALIGN = "align-2_align"
    ALIGN_DEDUP = "align-3_dedup"
    ALIGN_SORT = "align-4_sort"
    ALIGN_SPLIT = "align-5_split"
    VECTOR_BAMS = "vector-0_bams"
    VECTOR_SELECT = "vector-1_select"
    VECTOR_SORT = "vector-2_sort"


class Fastqc(Enum):
    QC_INPUT = "qc-inp"
    QC_TRIM = "qc-cut"


# File extensions

EXT_PATTERN = "([.].+)"
FASTA_EXTS = (".fasta", ".fna", ".fa")
BOWTIE2_INDEX_EXTS = (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                      ".rev.1.bt2", ".rev.2.bt2")
BOWTIE2_INDEX_PATTERN = "([.][1-4][.]bt2|[.]rev[.][1-2][.]bt2)"
FQ_EXTS = (".fastq", ".fq", ".fastq.gz", ".fq.gz",
           "_001.fastq", "_001.fq", "_001.fastq.gz", "_001.fq.gz")
FQ_PAIRED_EXTS_TEMPLATES = ("_R{}{}", "_mate{}{}", "_{}_sequence{}")
FQ1_EXTS = tuple(template.format(1, ext) for template, ext in
                 itertools.product(FQ_PAIRED_EXTS_TEMPLATES, FQ_EXTS))
FQ2_EXTS = tuple(template.format(2, ext) for template, ext in
                 itertools.product(FQ_PAIRED_EXTS_TEMPLATES, FQ_EXTS))
SAM_EXT = ".sam"
BAM_EXT = ".bam"
CRAM_EXT = ".cram"
XAM_EXTS = (SAM_EXT, BAM_EXT, CRAM_EXT)
BAI_EXT = f"{BAM_EXT}.bai"
XAI_EXTS = (BAI_EXT,)


# Exceptions ###########################################################

class PathError(Exception):
    """ Any error involving a path """


class PathTypeError(PathError, TypeError):
    """ Use of the wrong type of path or segment """


class PathValueError(PathError, ValueError):
    """ Invalid value of a path segment field """


# Path functions #######################################################

def sanitize(path: Any):
    return os.path.realpath(os.path.normpath(os.path.abspath(str(path))))


# Path segment classes #################################################

class BaseSeg(BaseModel):
    """ Abstract base class for all classes that represent segments of a path.
    Should not be instantiated. """

    class Config:
        extra = Extra.forbid

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
        return dict(zip(cls.keys(), values, strict=True))

    def values(self):
        return list(self.dict().values())

    @classmethod
    def parse(cls, segment: Any):
        """ Return a new instance of the class given a segment as a str-like
        argument. Non-structured segments (including this base class) must have
        exactly one field. """
        return cls(**cls.build_dict([segment]))

    def __str__(self):
        """ Return a string representation of the first field of the
        segment. Non-structured segments (including this base class)
        must have exactly one field. """
        try:
            value, = self.values()
        except ValueError:
            raise PathTypeError(f"str is undefined for {self.__class__} "
                                f"with {self.n_fields()} (≠ 1) fields.")
        return str(value)


class TopSeg(BaseSeg):
    """ Top-level working directory of DREEM. """
    top: str

    @validator(TOP_KEY)
    def sanitize_top_dir(cls, top: Any):
        return sanitize(top)


class SubSeg(BaseSeg):
    """ Any path level below the top level. """

    @root_validator(pre=True)
    def fields_must_match_valid_field_regex(cls, values: dict[str, Any]):
        for key, value in values.items():
            if isinstance(value, Enum):
                # For Enum types, the value must be obtained from the
                # value attribute, because simply calling str(value) on
                # the original type will yield the string {NAME}.{value}
                value = value.value
            if (key != TOP_KEY and key != EXT_KEY
                    and not VALID_FIELD_REGEX.match(str(value))):
                raise PathValueError(f"{cls} got invalid '{key}' value: '{value}'")
        return values


class SampleSeg(SubSeg):
    """ Segment for a sample. """
    sample: str


class RefsetSeg(SubSeg):
    """ Segment for a set of reference sequences. """
    refset: str


class OneRefSeg(SubSeg):
    """ Segment for one reference sequence. """
    ref: str


class ModSeg(SubSeg):
    """ Segment for the outputs of a module. """
    module: Module

    def __str__(self):
        return str(self.module.value)


class StepSeg(SubSeg):
    """ Segment for the temporary files of a step. """
    step: Step

    def __str__(self):
        return str(self.step.value)


class FastqcSeg(SubSeg):
    """ Segment for a directory of FASTQC reports. """
    fastqc: Fastqc

    def __str__(self):
        return str(self.fastqc.value)


class BatchSeg(SubSeg):
    """ Segment for a batch. """
    batch: NonNegativeInt


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
    ...     format_str = "samp-{}_batch-{}"
    ...     pattern_str = f"samp-{VALID_FIELD_PATTERN}_batch-([0-9]+)"
    >>> str(SampleBatchSeg(sample="tmv", batch=37))
    'samp-tmv_batch-37'
    """

    format_str: ClassVar[str] = ""
    pattern_str: ClassVar[str] = ""

    @classmethod
    def _parse(cls, seg_str: Any):
        """ Return a new instance of the class by parsing a ```str```-like
        ```segstr``` with ```re.match``` (using ```pattern_str``` as the
        regular expression pattern), then passing the parsed groups to the
        class constructor. """
        # Try to parse segstr with re.match, using pattern_str as the pattern.
        if match := re.match(cls.pattern_str, str(seg_str)):
            # If segstr matched the pattern, then create keyword arguments from
            # the matched values, and use them to initialize a new instance.
            return cls(**cls.build_dict(match.groups()))
        # If segstr did not match the pattern, then the parsing failed.
        raise PathValueError(
            f"Segment '{seg_str}' failed to match pattern {cls.pattern_str}")

    @classmethod
    def parse(cls, seg_str: Any):
        """
        Return a new instance of the class by parsing the values of the
        fields from a string and passing them to the class constructor.
        Parameters
        ----------
        seg_str: Any
            Literal representation of the segment as a ```str``` (or as
            any other object ```x``` for which ```str(x)``` is a valid
            path segment, e.g. an instance of ```pathlib.Path```).
        
        Returns
        -------
        StructSeg
            New instance of StructSeg with fields from ```seg_str```;
            satisfies ```str(seg_inst) == str(seg_str)```.
        Raises
        ------
        ValueError
            if the string representation of the newly created instance
            of the class did not match ```seg_str```.
        """
        # Create a new instance of the class by parsing seg_str.
        seg_inst = cls._parse(seg_str)
        # Confirm that formatting seg_inst yields the original seg_str.
        if (new_str := seg_inst._format()) != str(seg_str):
            # If not, an error occurred during parsing.
            raise PathValueError(
                f"The new instance was formatted as '{new_str}'"
                f"(failed to match input '{seg_str}')")
        # If so, the new instance was parsed correctly.
        return seg_inst

    def _format(self):
        """ Return the result of calling the ```format``` method of
        ```format_str``` with the field values as arguments. """
        return self.format_str.format(*self.values())

    def __str__(self) -> str:
        """
        Return a string representation of the path segment by inserting the
        values of the segment's fields into its format specification.
        """
        # Create the formatted string.
        seg_str = self._format()
        # Try to parse that string into a new instance, then confirm that the
        # fields of the new instance match those of the self instance.
        if (parse := self._parse(seg_str)) != self:
            # Raise an error if the fields did not match.
            raise PathValueError(
                f"String representation '{seg_str}' was parsed as "
                f"{repr(parse)} (failed to match {repr(self)})")
        # If the fields match, then the formatted string is valid.
        return seg_str


class SectionSeg(StructSeg):
    """
    Segment for a directory of a section of a reference sequence.
    Fields
    ------
    end5: int
        The 5'-most coordinate in the section; 1-indexed, inclusive.
    end3: int
        The 3'-most coordinate in the section; 1-indexed, inclusive.
    """

    end5: PositiveInt
    end3: PositiveInt

    format_str = "{}-{}"
    pattern_str = "([0-9]+)-([0-9]+)"

    @root_validator()
    def end5_le_end3(cls, values):
        """ Validate that end5 ≤ end3 """
        if values["end5"] > values["end3"]:
            raise PathValueError(f"Got end5 ({values['end5']}) "
                                 f"> end3 ({values['end3']})")
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
    pattern_str = f"^{VALID_FIELD_PATTERN}{EXT_PATTERN}$"

    @validator(EXT_KEY)
    def valid_file_extension(cls, ext):
        """ Validate the file extension (```ext```). It must be an
        element of the class attribute ```exts```. """
        if ext not in cls.exts:
            raise PathValueError(f"Invalid extension for {cls}: '{ext}'")
        return ext

    @classmethod
    def keys_no_ext(cls):
        """ Return a list of field keys without the 'ext' field. """
        # Create a new list of the field keys.
        keys = cls.keys()
        # Remove 'ext' from the list.
        keys.remove(EXT_KEY)
        # Return the list of keys without 'ext'.
        return keys

    @classmethod
    def _parse(cls, segstr: Any):
        if all(ext[0] not in VALID_CHARS_SET for ext in cls.exts):
            # If the extension starts with a character that is not valid if it
            # occurs in a field (e.g. the most common starting character '.'),
            # then the extension can be distinguished from the other fields by
            # the regular expression defined in the superclass.
            return super()._parse(segstr)
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
                # See if the truncated segment matches the truncated pattern.
                if match := re.match(pattern_trunc, segstr_trunc):
                    # If the truncated strings matched, then add the extension
                    # back to the end of the matched groups; we already know it
                    # matched because segstr_.endswith(ext) returned True.
                    return cls(**cls.build_dict(match.groups() + (ext,)))
        # If segstr did not match any extension, or did match an extension but
        # failed to match any patterns after removing the extension, then the
        # parsing failed.
        raise PathValueError(
            f"Segment '{segstr}' failed to match pattern "
            f"'{cls.pattern_str}' with any extension ({cls.exts})")


class MutVectorReportFileSeg(FileSeg):
    """ Segment for a mutation vector report file. """
    format_str = "report{}"
    pattern_str = f"report{EXT_PATTERN}"
    exts = (".json",)


class MutVectorBatchFileSeg(FileSeg, BatchSeg):
    """ Segment for a mutation vector batch file. """
    format_str = "{}{}"
    pattern_str = "([0-9]+)" + EXT_PATTERN
    exts = (".orc",)


class AbstractRefFileSeg(FileSeg):
    """ Abstract segment for a FASTA file. Do not instantiate. """
    exts = FASTA_EXTS


class RefsetSeqFileSeg(AbstractRefFileSeg, RefsetSeg):
    """ Segment for a FASTA file containing one or more references. The
    refset field is the name of the file, minus the file extension. """


class OneRefSeqFileSeg(AbstractRefFileSeg, OneRefSeg):
    """ Segment for a FASTA file that contains exactly one reference
    sequence. The ref field is the reference name (on the first line,
    after '>'). """


class AbstractBowtie2IndexFileSeg(FileSeg):
    """ Abstract segment for a Bowtie2 index. Do not instantiate. """
    exts = BOWTIE2_INDEX_EXTS
    pattern_str = f"^{VALID_FIELD_PATTERN}{BOWTIE2_INDEX_PATTERN}$"


class RefsetBowtie2IndexFileSeg(AbstractBowtie2IndexFileSeg, RefsetSeg):
    """ Segment for a Bowtie2 index of a FASTA file containing one or
    more references. """


class OneRefBowtie2IndexFileSeg(AbstractBowtie2IndexFileSeg, OneRefSeg):
    """ Segment for a Bowtie2 index of a FASTA file that contains
    exactly one reference sequence. """


class AbstractReadsFileSeg(FileSeg):
    """ Abstract segment for a standalone (single-end or paired-
    interleaved) FASTQ file. Do not instantiate. """
    exts = FQ_EXTS


class SampleReadsFileSeg(AbstractReadsFileSeg, SampleSeg):
    """ Segment for a standalone FASTQ file containing reads from one
    sample. """


class OneRefReadsFileSeg(AbstractReadsFileSeg, OneRefSeg):
    """ Segment for a standalone FASTQ file containing reads from one
    construct (i.e. a synthetic, barcoded reference) in one sample. """


class AbstractReads1FileSeg(FileSeg):
    """ Abstract segment for a paired-end FASTQ file containing mate-1 reads.
    Should not be instantiated. """
    exts = FQ1_EXTS


class SampleReads1FileSeg(AbstractReads1FileSeg, SampleSeg):
    """ Segment for a paired-end FASTQ file containing all mate-1 reads from a
    sample. """


class OneRefReads1FileSeg(AbstractReads1FileSeg, OneRefSeg):
    """ Segment for a paired-end FASTQ file containing all mate-1 reads from one
    reference (i.e. reference sequence) in one sample. """


class AbstractReads2FileSeg(FileSeg):
    """ Abstract segment for a paired-end FASTQ file containing mate-2 reads.
    Should not be instantiated. """
    exts = FQ2_EXTS


class SampleReads2FileSeg(AbstractReads2FileSeg, SampleSeg):
    """ Segment for a paired-end FASTQ file containing all mate-2 reads from a
    sample. """


class OneRefReads2FileSeg(AbstractReads2FileSeg, OneRefSeg):
    """ Segment for a paired-end FASTQ file containing all mate-2 reads
    from one construct (i.e. reference sequence) in one sample. """


class AbstractAlignmentFileSeg(FileSeg):
    """ Abstract segment for alignment map files (.sam, .bam, .cram) """
    exts = XAM_EXTS


class RefsetAlignmentFileSeg(AbstractAlignmentFileSeg, RefsetSeg):
    """ Segment for an ailgnment map file resulting from aligning a set
    of reads to a set of references """


class OneRefAlignmentFileSeg(AbstractAlignmentFileSeg, OneRefSeg):
    """ Segment for an alignment map file produced by splitting another
    alignment map file into one file for each reference """


class SectionAlignmentFileSeg(AbstractAlignmentFileSeg, SectionSeg):
    """ Segment for an alignment map file produced by taking a subset
    of a particular section from another alignment map file """
    format_str = "{}-{}{}"
    pattern_str = f"([0-9]+)-([0-9]+){EXT_PATTERN}"


class AbstractAlignmentIndexFileSeg(FileSeg):
    """ Abstract segment for index files of alignment maps """
    exts = XAI_EXTS


class RefsetAlignmentIndexFileSeg(AbstractAlignmentIndexFileSeg, RefsetSeg):
    """ Segment for an index file of an ailgnment map file resulting
    from aligning a set of reads to a set of references """


class OneRefAlignmentIndexFileSeg(AbstractAlignmentIndexFileSeg, OneRefSeg):
    """ Segment for an index file of an ailgnment map file resulting
    from aligning a set of reads to a set of references """


# Path classes #########################################################

class BasePath(BaseModel):
    """ Abstract base class for all classes that represent full,
    absolute paths. Do not instantiate. """

    @classmethod
    def fields(cls):
        return cls.__fields__

    @classmethod
    def keys(cls):
        return list(cls.fields().keys())

    @classmethod
    def _get_segment_types(cls):
        """ Return a list of the type of every segment in the path. """
        seg_types: list[type[BaseSeg]] = list()
        # Search the base classes for path segments.
        for base in cls.__bases__:
            if issubclass(base, BasePath):
                # The base class represents a path and thus also has
                # base classes of path segments: add them to the list.
                seg_types.extend(base._get_segment_types())
            elif issubclass(base, BaseSeg):
                # The base class is a path segment but not a path class.
                # Thus, it should be in the list of segments.
                seg_types.append(base)
            else:
                # Otherwise, the base class contains no path segments.
                pass
        return seg_types

    @classmethod
    def segment_types(cls) -> tuple[type[BaseSeg], ...]:
        """ Return a tuple of the type of every segment in the path, in
        order from the beginning to the end of the path. """
        seg_types = tuple(cls._get_segment_types())
        if seg_types and seg_types[0] is not TopSeg:
            raise PathTypeError(f"{cls} begins with {seg_types[0].__name__}, "
                                f"not {TopSeg.__name__}")
        return seg_types

    @classmethod
    def _parse_segments(cls,
                        remaining_seg_types: list[type[BaseSeg]],
                        path: str):
        """
        Return a dict of the names and values of the fields encoded
        within a string representation of a path.
        Parameters
        ----------
        remaining_seg_types: list[type[BaseSeg]]
            Types of segments in the path that have yet to be used to
            parse part of the path string. The segment types are ordered
            from beginning (left) to end (right) of the path. The first
            item (left-most segment type) must be ```TopSeg```.
        path: str
            Path to be parsed, represented as a string.
        Return
        ------
        dict[str, Any]
            A dict of the names (keys) and values (values) of the fields
            encoded in the path string.
        Raises
        ------
        PathTypeError
            if any part of the path remains to be parsed after all the
            types of segments have been consumed.
        """
        try:
            # Get the last (right-most) type of segment in the path that
            # has not been used to parse a segment of the path string.
            rightmost_seg_type = remaining_seg_types.pop()
        except IndexError:
            # Base case of this recursive function. The pop() operation
            # raised an IndexError because seg_types is empty. Thus, no
            # segment types remain.
            if path:
                # Any part of the path that has not yet been parsed
                # cannot be parsed, since no segment types are left.
                raise PathTypeError(f"No segments remain to parse '{path}'")
            # Return a dict with no fields, signifying that nothing
            # remains to be parsed.
            return dict()
        # At least one more segment type remains.
        if remaining_seg_types:
            # If any segment types will remain after the one that is
            # currently the right-most segment type, then parse only the
            # last (right-most) segment of the path string (i.e. after
            # the final path separator) using the right-most segment
            # type. The remaining segments (towards the left side of the
            # path) will be parsed with the remaining segment types.
            remaining_segs, rightmost_seg = os.path.split(path)
        else:
            # If the segment type that is currently the right-most type
            # is the last one remaining, then use it to parse the entire
            # remaining part of the path. To do so, consider all that
            # remains of the path to be the "right-most segment" (which
            # will be parsed by the right-most segment type), and thus
            # assign an empty string to the part of the path that will
            # remain after parsing this "right-most segment".
            remaining_segs, rightmost_seg = "", path
        # Ensure that part of the path still remains to be parsed.
        if not path:
            raise PathValueError("No part of the path remains to be parsed by "
                                 f"segment '{rightmost_seg_type.__class__}'")
        # Parse the right-most segment with the right-most segment type,
        # and parse the remaining segments with the remaining segment
        # types. This step will also catch any duplicate fields (which
        # are not supposed to occur) by raising a TypeError.
        return dict(**cls._parse_segments(remaining_seg_types, remaining_segs),
                    **rightmost_seg_type.parse(rightmost_seg).dict())

    @classmethod
    def parse(cls, path: Any):
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
        PathValueError
            if the newly created instance of the class yields a string
            representation that does not match the ```path``` argument.
        """
        path_inst = cls(**cls._parse_segments(list(cls.segment_types()),
                                              str(path)))
        if str(path_inst) != os.path.abspath(path):
            raise PathValueError(
                f"String representation of new path '{path_inst}' "
                f"failed to match input '{os.path.abspath(path)}'")
        return path_inst

    @property
    def segments(self) -> tuple[BaseSeg, ...]:
        """ Return a tuple of an instance of every segment in the path,
        in order from the beginning to the end of the path. """
        return tuple(seg_type(**{key: getattr(self, key)
                                 for key in seg_type.keys()})
                     for seg_type in self.segment_types())

    @property
    def path(self):
        """ Return a ```pathlib.Path``` instance for this path. """
        # Using @functools.cached_property instead of @functools.cache
        # would be more efficient but causes errors when applied to
        # methods of Pydantic models.
        return pathlib.Path(*map(str, self.segments))

    @property
    def parent(self):
        """ Return the parent directory of the path. """
        parent_fields = dict()
        for seg in self.segments[:-1]:
            parent_fields.update(seg.dict())
        return create(**parent_fields)

    def _get_new_fields(self, **fields):
        """ Computes the new fields obtained by merging given fields
        with the existing fields of this path. If a field occurs in both
        the given and the existing fields, the given field replaces the
        existing field. """
        return {**self.dict(), **fields}

    def edit(self, **fields):
        """ Return a new path instance based on this instance with some
        fields modified according to the keyword arguments. Set a field
        to ```None``` to delete it from the new path instance.
        Return
        ------
        BasePath
            Path: modified; Type: modified
        """
        return create(**self._get_new_fields(**fields))

    def cast(self, **fields):
        """ Return a new path instance that evaluates to the same path
        string and has been cast to a new type. The new type is found by
        merging the given fields with the existing fields of this path
        and determining which type matches those fields. Set a field to
        ```None``` to delete it from the new path instance. Note that,
        except for 'ext' (which determines file types), the values of
        the given fields do not affect the return value, since the type
        of path does not depend on the values of the fields, but rather
        on which fields it has (again, except for the field 'ext').
        Return
        ------
        BasePath
            Path: preserved; Type: modified
        """
        new_type = get_path_class_by_fields(**self._get_new_fields(**fields))
        return new_type.parse(self)

    def move(self, **fields):
        """ Return a new path instance that has the same type as the
        original but has a different path. Determine the new path based
        on the fields given as keyword arguments.
        Return
        ------
        BasePath
            Path: modified; Type: preserved
        """
        # First, create the new path (modified path and modified type)
        # using the edit method.
        new_path = self.edit(**fields)
        # Compare the types directly here rather than using isinstance
        # so that subclasses are not considered equivalent to their
        # parent classes. See the comments in __eq__ for more details.
        if new_path.__class__ is self.__class__:
            return new_path
        # If, as usually happens, the type of the new path does not
        # match the type of the previous path, then parse the new path
        # with the type of the previous path.
        return self.__class__.parse(new_path)

    def __str__(self):
        return str(self.path)

    def __eq__(self, other: Any):
        # Note that this comparison avoids using isinstance because a
        # subclass should not evaluate as equal to its parent class even
        # if the path strings are identical. This is because, unlike in
        # most cases where a subclass is a "more specific" version of
        # its parent class, in this module, object inheritance is used
        # to create the hierarchical structure of directories and files,
        # wherein a subclass in the code equates to a child directory in
        # the file system. Since a child directory or file is not just a
        # "more specific" kind of its parent directory, subclasses and
        # parent classes are considered different for equality test_input.
        return self.__class__ is other.__class__ and str(self) == str(other)


# General directory paths

class TopDirPath(BasePath, TopSeg):
    pass


class ModuleDirPath(TopDirPath, ModSeg):
    pass


class StepDirPath(ModuleDirPath, StepSeg):
    pass


# Sample, Reference, and Section directory paths

class AbstractSampleDirPath(BasePath, SampleSeg):
    """ Abstract directory named after a sample """


class SampleInDirPath(TopDirPath, AbstractSampleDirPath):
    """ Directory of input files, named after a sample """


class SampleStepDirPath(StepDirPath, AbstractSampleDirPath):
    """ Directory of temporary files, named after a sample """


class SampleOutDirPath(ModuleDirPath, AbstractSampleDirPath):
    """ Directory of output files, named after a sample """


class AbstractRefDirPath(BasePath, OneRefSeg):
    """ Abstract directory named after a reference """


class RefInDirPath(SampleInDirPath, AbstractRefDirPath):
    """ Directory of input files, named after a reference """


class RefStepDirPath(SampleStepDirPath, AbstractRefDirPath):
    """ Directory of temporary files, named after a reference """


class RefOutDirPath(SampleOutDirPath, AbstractRefDirPath):
    """ Directory of output files, named after a reference """


class AbstractSectionDirPath(BasePath, SectionSeg):
    """ Abstract directory named after a section """


class SectionInDirPath(RefInDirPath, AbstractSectionDirPath):
    """ Directory of input files, named after a section """


class SectionStepDirPath(RefStepDirPath, AbstractSectionDirPath):
    """ Directory of temporary files, named after a section """


class SectionOutDirPath(RefOutDirPath, AbstractSectionDirPath):
    """ Directory of output files, named after a section """


# Reference sequence (FASTA) and Bowtie2 index file paths

class AbstractRefsetSeqFilePath(BasePath, RefsetSeqFileSeg):
    """ Abstract FASTA file named after a set of references """


class RefsetSeqInFilePath(TopDirPath, AbstractRefsetSeqFilePath):
    """ Input FASTA file named after a set of references """


class RefsetSeqTempFilePath(StepDirPath, AbstractRefsetSeqFilePath):
    """ Temporary FASTA file named after a set of references """


class RefsetSeqOutFilePath(ModuleDirPath, AbstractRefsetSeqFilePath):
    """ Output FASTA file named after a set of references """


class AbstractRefsetBowtie2IndexFilePath(BasePath, RefsetBowtie2IndexFileSeg):
    """ Abstract Bowtie2 index file for a FASTA file named after a set
    of references """


class RefsetBowtie2IndexInFilePath(TopDirPath,
                                   AbstractRefsetBowtie2IndexFilePath):
    """ Bowtie2 index file for an input FASTA file named after a set of
    references """


class RefsetBowtie2IndexTempFilePath(StepDirPath,
                                     AbstractRefsetBowtie2IndexFilePath):
    """ Bowtie2 index file for a temporary FASTA file named after a set
    of references """


class RefsetBowtie2IndexOutFilePath(ModuleDirPath,
                                    AbstractRefsetBowtie2IndexFilePath):
    """ Bowtie2 index file for an output FASTA file named after a set of
    references """


class AbstractRefsetBowtie2IndexPrefix(BasePath, RefsetSeg):
    """ Abstract Bowtie2 index prefix for a FASTA file named after a set
    of references """


class RefsetBowtie2IndexInPrefix(TopDirPath,
                                 AbstractRefsetBowtie2IndexPrefix):
    """ Bowtie2 index prefix for an input FASTA file named after a set
    of references """


class RefsetBowtie2IndexTempPrefix(StepDirPath,
                                   AbstractRefsetBowtie2IndexPrefix):
    """ Bowtie2 index prefix for a temporary FASTA file named after a
    set of references """


class RefsetBowtie2IndexOutPrefix(ModuleDirPath,
                                  AbstractRefsetBowtie2IndexPrefix):
    """ Bowtie2 index prefix for an output FASTA file named after a set
    of references """


class AbstractOneRefSeqFilePath(BasePath, OneRefSeqFileSeg):
    """ Abstract FASTA file named after one reference """


class OneRefSeqInFilePath(TopDirPath, AbstractOneRefSeqFilePath):
    """ Input FASTA file named after one reference """


class OneRefSeqTempFilePath(StepDirPath, AbstractOneRefSeqFilePath):
    """ Temporary FASTA file named after one reference """


class OneRefSeqOutFilePath(ModuleDirPath, AbstractOneRefSeqFilePath):
    """ Output FASTA file named after one reference """


class AbstractOneRefBowtie2IndexFilePath(BasePath, OneRefBowtie2IndexFileSeg):
    """ Abstract Bowtie2 index file for a FASTA file named after one
    reference """


class OneRefBowtie2IndexInFilePath(TopDirPath,
                                   AbstractOneRefBowtie2IndexFilePath):
    """ Bowtie2 index file for an input FASTA file named after one
    reference """


class OneRefBowtie2IndexTempFilePath(StepDirPath,
                                     AbstractOneRefBowtie2IndexFilePath):
    """ Bowtie2 index file for a temporary FASTA file named after one
    reference """


class OneRefBowtie2IndexOutFilePath(ModuleDirPath,
                                    AbstractOneRefBowtie2IndexFilePath):
    """ Bowtie2 index file for an output FASTA file named after one
    reference """


class AbstractOneRefBowtie2IndexPrefix(BasePath, OneRefSeg):
    """ Abstract Bowtie2 index prefix for a FASTA file named after one
    reference """


class OneRefBowtie2IndexInPrefix(TopDirPath,
                                 AbstractOneRefBowtie2IndexPrefix):
    """ Bowtie2 index file for an input FASTA file named after one
    reference """


class OneRefBowtie2IndexTempPrefix(StepDirPath,
                                   AbstractOneRefBowtie2IndexPrefix):
    """ Bowtie2 index file for a temporary FASTA file named after one
    reference """


class OneRefBowtie2IndexOutPrefix(ModuleDirPath,
                                  AbstractOneRefBowtie2IndexPrefix):
    """ Bowtie2 index file for an output FASTA file named after one
    reference """


# Sequencing reads (FASTQ) file paths

class AbstractSampleReadsFilePath(BasePath, SampleReadsFileSeg):
    """ Abstract FASTQ file named after a sample """


class SampleReadsInFilePath(TopDirPath, AbstractSampleReadsFilePath):
    """ Input FASTQ file named after a sample """


class SampleReadsTempFilePath(StepDirPath, AbstractSampleReadsFilePath):
    """ Temporary FASTQ file named after a sample """


class SampleReadsOutFilePath(ModuleDirPath, AbstractSampleReadsFilePath):
    """ Output FASTQ file named after a sample """


class AbstractSampleReads1FilePath(BasePath, SampleReads1FileSeg):
    """ Abstract FASTQ mate 1 file named after a sample """


class SampleReads1InFilePath(TopDirPath, AbstractSampleReads1FilePath):
    """ Input FASTQ mate 1 file named after a sample """


class SampleReads1TempFilePath(StepDirPath, AbstractSampleReads1FilePath):
    """ Temporary FASTQ mate 1 file named after a sample """


class SampleReads1OutFilePath(ModuleDirPath, AbstractSampleReads1FilePath):
    """ Output FASTQ mate 1 file named after a sample """


class AbstractSampleReads2FilePath(BasePath, SampleReads2FileSeg):
    """ Abstract FASTQ mate 2 file named after a sample """


class SampleReads2InFilePath(TopDirPath, AbstractSampleReads2FilePath):
    """ Input FASTQ mate 2 file named after a sample """


class SampleReads2TempFilePath(StepDirPath, AbstractSampleReads2FilePath):
    """ Temporary FASTQ mate 2 file named after a sample """


class SampleReads2OutFilePath(ModuleDirPath, AbstractSampleReads2FilePath):
    """ Output FASTQ mate 2 file named after a sample """


class AbstractOneRefReadsFilePath(BasePath, OneRefReadsFileSeg):
    """ Abstract FASTQ file named after one reference """


class OneRefReadsInFilePath(SampleInDirPath, AbstractOneRefReadsFilePath):
    """ Input FASTQ file named after one reference """


class OneRefReadsTempFilePath(SampleStepDirPath, AbstractOneRefReadsFilePath):
    """ Temporary FASTQ file named after one reference """


class OneRefReadsOutFilePath(SampleOutDirPath, AbstractOneRefReadsFilePath):
    """ Output FASTQ file named after one reference """


class AbstractOneRefReads1FilePath(BasePath, OneRefReads1FileSeg):
    """ Abstract FASTQ mate 1 file named after one reference """


class OneRefReads1InFilePath(SampleInDirPath, AbstractOneRefReads1FilePath):
    """ Input FASTQ mate 1 file named after one reference """


class OneRefReads1TempFilePath(SampleStepDirPath, AbstractOneRefReads1FilePath):
    """ Temporary FASTQ mate 1 file named after one reference """


class OneRefReads1OutFilePath(SampleOutDirPath, AbstractOneRefReads1FilePath):
    """ Output FASTQ mate 1 file named after one reference """


class AbstractOneRefReads2FilePath(BasePath, OneRefReads2FileSeg):
    """ Abstract FASTQ mate 2 file named after one reference """


class OneRefReads2InFilePath(SampleInDirPath, AbstractOneRefReads2FilePath):
    """ Input FASTQ mate 2 file named after one reference """


class OneRefReads2TempFilePath(SampleStepDirPath, AbstractOneRefReads2FilePath):
    """ Temporary FASTQ mate 2 file named after one reference """


class OneRefReads2OutFilePath(SampleOutDirPath, AbstractOneRefReads2FilePath):
    """ Output FASTQ mate 2 file named after one reference """


# FASTQ quality control (FASTQC) paths

class FastqcOutDirPath(SampleOutDirPath, FastqcSeg):
    """ Output directory of FASTQC reports """


# Alignment map (SAM/BAM/CRAM) and index (BAM.BAI) file paths

class AbstractRefsetAlignmentFilePath(BasePath, RefsetAlignmentFileSeg):
    """ Abstract alignment map file named after a set of references """


class RefsetAlignmentInFilePath(SampleInDirPath,
                                AbstractRefsetAlignmentFilePath):
    """ Input alignment map file named after a set of references """


class RefsetAlignmentTempFilePath(SampleStepDirPath,
                                  AbstractRefsetAlignmentFilePath):
    """ Temporary alignment map file named after a set of references """


class RefsetAlignmentOutFilePath(SampleOutDirPath,
                                 AbstractRefsetAlignmentFilePath):
    """ Output alignment map file named after a set of references """


class AbstractRefsetAlignmentIndexFilePath(BasePath,
                                           RefsetAlignmentIndexFileSeg):
    """ Abstract alignment map index file named after a set of
    references """


class RefsetAlignmentIndexInFilePath(SampleInDirPath,
                                     AbstractRefsetAlignmentIndexFilePath):
    """ Input alignment map index file named after a set of
    references """


class RefsetAlignmentIndexTempFilePath(SampleStepDirPath,
                                       AbstractRefsetAlignmentIndexFilePath):
    """ Temporary alignment map index file named after a set of
    references """


class RefsetAlignmentIndexOutFilePath(SampleOutDirPath,
                                      AbstractRefsetAlignmentIndexFilePath):
    """ Output alignment map index file named after a set of
    references """


class AbstractOneRefAlignmentFilePath(BasePath, OneRefAlignmentFileSeg):
    """ Abstract alignment map file named after one reference """


class OneRefAlignmentInFilePath(SampleInDirPath,
                                AbstractOneRefAlignmentFilePath):
    """ Input alignment map file named after one reference """


class OneRefAlignmentTempFilePath(SampleStepDirPath,
                                  AbstractOneRefAlignmentFilePath):
    """ Temporary alignment map file named after one reference """


class OneRefAlignmentOutFilePath(SampleOutDirPath,
                                 AbstractOneRefAlignmentFilePath):
    """ Output alignment map file named after one reference """


class AbstractOneRefAlignmentIndexFilePath(BasePath,
                                           OneRefAlignmentIndexFileSeg):
    """ Abstract alignment map index file named after one reference """


class OneRefAlignmentIndexInFilePath(SampleInDirPath,
                                     AbstractOneRefAlignmentIndexFilePath):
    """ Input alignment map index file named after one reference """


class OneRefAlignmentIndexTempFilePath(SampleStepDirPath,
                                       AbstractOneRefAlignmentIndexFilePath):
    """ Temporary alignment map index file named after one reference """


class OneRefAlignmentIndexOutFilePath(SampleOutDirPath,
                                      AbstractOneRefAlignmentIndexFilePath):
    """ Output alignment map index file named after one reference """


class AbstractSectionAlignmentFilePath(BasePath, SectionAlignmentFileSeg):
    """ Abstract alignment map file named after a section """
    ref: str


class SectionAlignmentInFilePath(RefInDirPath,
                                 AbstractSectionAlignmentFilePath):
    """ Input alignment map file named after a section """


class SectionAlignmentTempFilePath(RefStepDirPath,
                                   AbstractSectionAlignmentFilePath):
    """ Temporary alignment map file named after a section """


class SectionAlignmentOutFilePath(RefOutDirPath,
                                  AbstractSectionAlignmentFilePath):
    """ Output alignment map file named after a section """


# Vectoring file paths

class MutVectorBatchFilePath(SectionOutDirPath, MutVectorBatchFileSeg):
    """ Output file of a batch of mutation vectors """


class MutVectorReportFilePath(SectionOutDirPath, MutVectorReportFileSeg):
    """ Output vectorization report file """


# Path collection and conversion functions #############################

def is_concrete_path_class(member: Any):
    """ Return whether ```member``` is a class of concrete path. """
    # To be a class of concrete path (i.e. instantiable, with at least
    # the top-level field and all essential path methods implemented),
    # a member of the module must 1) be a class and 2) be TopDirPath or
    # a subclass of TopDirPath (issubclass returns True for both), which
    # provides the top-level segment and all required methods of paths.
    return isclass(member) and issubclass(member, TopDirPath)


@cache
def _generate_name_to_path_class() -> dict[str, type[BasePath]]:
    """ Return a ```dict``` of every valid, instantiable class of path
    in this module, with class names as keys. """
    return dict(getmembers(sys.modules[__name__], is_concrete_path_class))


@cache
def _generate_fields_exts_to_path_class():
    """ Return a ```dict``` of every valid, instantiable class of path.
    Each key is a ```tuple``` of field names in alphabetical order.
    For directories, the field names uniquely determine the path class.
    For files, the extension also determines the path class; thus, each
    key of field names maps to another ```dict``` keyed by extensions.
    One entry each for a directory and set of files are illustrated:
    >>> _ = {
    ...     # Fields 'top' and 'module' map uniquely to ModuleDirPath.
    ...     ('module', 'top'): ModuleDirPath,
    ...     # Fields 'top', 'module', 'sample', 'ref', 'end5', 'end3',
    ...     # and 'ext' map to four file extensions that map to a total
    ...     # of two classes of files (```SectionAlignmentOutFilePath```
    ...     # and ```MutVectorReportFilePath```).
    ...     ('end3', 'end5', 'ext', 'module', 'ref', 'sample', 'top'): {
    ...         '.json': MutVectorReportFilePath,
    ...         '.bam': SectionAlignmentOutFilePath,
    ...         '.cram': SectionAlignmentOutFilePath,
    ...         '.sam': SectionAlignmentOutFilePath,
    ...     }
    ... }
    """
    # Dictionary that maps fields and extensions to path classes
    fe_to_cls = dict()
    # Iterate over every valid class of path.
    for cls in _generate_name_to_path_class().values():
        # Get the field names of the class in alphabetical order.
        fields = tuple(sorted(cls.keys()))
        if EXT_KEY in fields:
            # The class has a field for file extensions.
            if fields not in fe_to_cls:
                # If the dictionary of extensions for this set of fields
                # does not yet exist, initialize an empty one.
                fe_to_cls[fields]: dict[str, type[BasePath]] = dict()
            if not isinstance(fe_to_cls[fields], dict):
                raise TypeError(f"Expected type 'dict' at {fields}, but got "
                                f"'{type(fe_to_cls[fields]).__name__}'")
            # For each file extension of the class, add the class to the
            # dictionary using the extension as the key.
            for ext in getattr(cls, EXTS_KEY):
                if ext in fe_to_cls[fields]:
                    raise ValueError(f"{fe_to_cls[fields]} and {cls} have "
                                     f"identical fields {fields} and file "
                                     f"extension '{ext}'")
                fe_to_cls[fields][ext] = cls
        else:
            # The class does not have a field for file extensions.
            if fields in fe_to_cls:
                raise TypeError(f"Cannot insert {cls.__name__} at {fields}: "
                                f"{fe_to_cls[fields]} is already there")
            # Add the class to the dictionary, keyed only by the fields.
            fe_to_cls[fields] = cls
    return fe_to_cls


def get_path_class_by_name(name: str):
    """ Return the path class with the given name. """
    try:
        return _generate_name_to_path_class()[name]
    except KeyError:
        raise ValueError(f"No path class named '{name}'")


def get_path_class_by_fields_ext(fields: Iterable[str],
                                 ext: Iterable[str] | str | None = None):
    """ Return the path class with the given fields and, if applicable,
    file extension. If no file extension is given, the class is assumed
    to be a directory. """
    # Get the fields in standard form (tuple, alphabetical order).
    std_fields = tuple(sorted(fields))
    if ext is None:
        # Assume the class is a directory path.
        try:
            return _generate_fields_exts_to_path_class()[std_fields]
        except KeyError:
            raise ValueError(f"No class of directory with fields {std_fields}")
    # The class is not a directory: assume it is a file path.
    if isinstance(ext, str):
        # Search for the file path given a single file extension.
        try:
            return _generate_fields_exts_to_path_class()[std_fields][ext]
        except KeyError:
            raise ValueError(f"No class of file with fields {std_fields} and "
                             f"extension '{ext}'")
    # Search for the file path given multiple file extensions.
    # First, find the type of path corresponding to each extension.
    # Keep at most one of each type of path by storing them in a set
    # before converting the set to a list.
    classes = list({get_path_class_by_fields_ext(std_fields, e) for e in ext})
    # Verify that exactly one unique type of path corresponds to all
    # the given file extensions.
    if not classes:
        raise ValueError(f"No path class with fields {std_fields} and "
                         f"extensions {ext}: {classes}")
    if len(classes) > 1:
        raise ValueError(f"Multiple path classes with fields {std_fields} and "
                         f"extensions {ext}: {classes}")
    # Return the type of path corresponding to all the file extensions.
    return classes[0]


def clean_fields(fields: dict[str, Any]):
    """ Remove each field whose value is ```None```. """
    return {name: value for name, value in fields.items() if value is not None}


def get_path_class_by_fields(**fields):
    """ Return the path class that corresponds to the given fields. Like
    ```get_path_class_by_fields_ext```, but here all fields are given as
    keyword arguments rather than an iterable of fields and a separate
    argument for the optional file extension. """
    return get_path_class_by_fields_ext(clean_fields(fields),
                                        fields.get(EXT_KEY))


def create(**fields):
    """ Create a new path instance with the given fields. """
    return get_path_class_by_fields(**fields)(**clean_fields(fields))


# Ensure that there are no two paths with identical fields.
_generate_fields_exts_to_path_class()