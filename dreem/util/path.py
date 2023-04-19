"""
Path module of DREEM
========================================================================
Auth: Matty
Date: 2023-04-07

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
"""

# Imports ##############################################################

from __future__ import annotations
from collections import Counter
from functools import cache
from itertools import chain, product
from logging import getLogger
import os
import pathlib as pl
import re
from string import ascii_letters, digits, printable
from typing import Any, Iterable, Sequence


# Constants ############################################################

logger = getLogger(__name__)

# Valid/invalid characters in fields

PATH_CHARS = printable
STR_CHARS = ascii_letters + digits + "_.=+-"
STR_CHARS_SET = frozenset(STR_CHARS)
INT_CHARS = digits
PATH_PATTERN = f"([{PATH_CHARS}]+)"
STR_PATTERN = f"([{STR_CHARS}]+)"
INT_PATTERN = f"([{INT_CHARS}]+)"
RE_PATTERNS = {str: STR_PATTERN, int: INT_PATTERN, pl.Path: PATH_PATTERN}

MOD_DMPL = "demultiplex"
MOD_ALGN = "alignment"
MOD_VECT = "vectoring"
MOD_CLST = "clustering"
MOD_AGGR = "aggregate"
MODULES = MOD_DMPL, MOD_ALGN, MOD_VECT, MOD_CLST, MOD_AGGR

STEPS_FSQC = "qc-inp", "qc-trim"
STEPS_ALGN = ("align-0_refs", "align-1_trim", "align-2_align",
              "align-3_dedup", "align-4_sort", "align-5_split")
STEPS_VECT = "vector-0_bams",
STEPS = STEPS_FSQC + STEPS_ALGN + STEPS_VECT

# File extensions

CSV_EXT = ".csv"
ORC_EXT = ".orc"
JSON_EXT = ".json"
FASTA_EXTS = ".fasta", ".fna", ".fa"
BOWTIE2_INDEX_EXTS = (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2",
                      ".rev.1.bt2", ".rev.2.bt2")
FQ_EXTS = (".fastq", ".fq", ".fastq.gz", ".fq.gz",
           "_001.fastq", "_001.fq", "_001.fastq.gz", "_001.fq.gz")
FQ_PAIRED_EXTS_TEMPLATES = ("_R{}{}", "_mate{}{}", "_{}_sequence{}")
FQ1_EXTS = tuple(template.format(1, ext) for template, ext in
                 product(FQ_PAIRED_EXTS_TEMPLATES, FQ_EXTS))
FQ2_EXTS = tuple(template.format(2, ext) for template, ext in
                 product(FQ_PAIRED_EXTS_TEMPLATES, FQ_EXTS))
SAM_EXT = ".sam"
BAM_EXT = ".bam"
XAM_EXTS = SAM_EXT, BAM_EXT
BAI_EXT = f"{BAM_EXT}.bai"


# Path Exceptions ######################################################

class PathError(Exception):
    """ Any error involving a path """


class PathTypeError(PathError, TypeError):
    """ Use of the wrong type of path or segment """


class PathValueError(PathError, ValueError):
    """ Invalid value of a path segment field """


# Path Functions #######################################################

def sanitize(path: str | pl.Path):
    return pl.Path(os.path.realpath(os.path.normpath(os.path.abspath(path))))


# Path Fields ##########################################################

# Field validation functions

def validate_str(txt: str):
    if not isinstance(txt, str):
        raise PathTypeError(f"Expected 'str', got '{type(txt).__name__}'")
    if illegal := "".join(sorted(set(txt) - STR_CHARS_SET)):
        raise PathValueError(f"Text '{txt}' has illegal characters: {illegal}")


def validate_top(top: pl.Path):
    if not isinstance(top, pl.Path):
        raise PathTypeError(f"Expected 'Path', got '{type(top).__name__}'")
    if not top.parent.is_dir():
        raise PathValueError(f"Not a directory: {top.parent}")
    if top.is_file():
        raise PathValueError(f"File exists: {top}")


def validate_int(num: int):
    if num <= 0:
        raise PathValueError(f"Expected positive integer, got {num}")


VALIDATE = {int: validate_int,
            str: validate_str,
            pl.Path: validate_top}


# Field class

class Field(object):
    def __init__(self, /,
                 dtype: type[str | int | pl.Path],
                 options: Iterable | None = None,
                 is_ext: bool = False):
        self.dtype = dtype
        self.options = list() if options is None else list(options)
        if not all(isinstance(option, self.dtype) for option in self.options):
            raise PathTypeError("All options of a field must be of its type")
        self.is_ext = is_ext
        if self.is_ext:
            if self.dtype is not str:
                raise PathTypeError("Extension field must be type 'str', "
                                    f"but got type '{self.dtype.__name__}'")
            if not self.options:
                raise PathValueError("Extension field must have options")

    def validate(self, val: Any):
        if not isinstance(val, self.dtype):
            raise PathTypeError(f"Expected '{self.dtype.__name__}', "
                                f"but got '{type(val).__name__}'")
        if self.options and val not in self.options:
            raise PathValueError(
                f"Invalid option '{val}'; expected one of {self.options}")
        VALIDATE[self.dtype](val)

    def build(self, val: Any):
        """ Validate a value and return it as a string. """
        self.validate(val)
        return str(val)

    def parse(self, text: str) -> Any:
        """ Parse a value from a string, validate it, and return it. """
        try:
            val = self.dtype(text)
        except Exception as error:
            raise PathValueError(f"Failed to interpret '{text}' as type "
                                 f"'{self.dtype.__name__}': {error}")
        self.validate(val)
        return val

    def __str__(self):
        return f"Field({self.dtype.__name__})"


# Fields
TopField = Field(pl.Path)
NameField = Field(str)
ModField = Field(str, MODULES)
StepField = Field(str, STEPS)
IntField = Field(int)

# File extensions
# CsvExt = Field(str, [CSV_EXT], is_ext=True)
VecBatExt = Field(str, [ORC_EXT], is_ext=True)
VecRepExt = Field(str, [JSON_EXT], is_ext=True)
FastaExt = Field(str, FASTA_EXTS, is_ext=True)
FastaIndexExt = Field(str, BOWTIE2_INDEX_EXTS, is_ext=True)
FastqExt = Field(str, FQ_EXTS, is_ext=True)
Fastq1Ext = Field(str, FQ1_EXTS, is_ext=True)
Fastq2Ext = Field(str, FQ2_EXTS, is_ext=True)
XamExt = Field(str, XAM_EXTS, is_ext=True)
BamIndexExt = Field(str, [BAI_EXT], is_ext=True)


# Path Segments ########################################################

# Segment class

class Segment(object):
    def __init__(self, /,
                 field_types: dict[str, Field], *,
                 frmt: str | None = None):
        self.field_types = field_types
        if not self.field_types:
            raise PathValueError(f"Segment got no fields")
        # Verify that a field has the key EXT if and only if it is an
        # extension and is the last field in the segment.
        for i, (name, field) in enumerate(self.field_types.items(), start=1):
            if name == EXT:
                if not field.is_ext:
                    raise PathValueError(f"Field '{EXT}' is not an extension")
                if i != len(self.field_types):
                    raise PathValueError(f"Extension is not the last field")
            elif field.is_ext:
                raise PathValueError(f"Extension field has name '{name}'")
        # Determine the format string.
        if frmt is None:
            # Default format is to concatenate all the fields.
            frmt = "".join("{" + name + "}" for name, field
                           in self.field_types.items())
        self.frmt = frmt
        # Generate the pattern string using the format string, excluding
        # the extension (if any) because before parsing, it is removed
        # from the end of the string.
        patterns = {name: "" if field.is_ext else RE_PATTERNS[field.dtype]
                    for name, field in self.field_types.items()}
        self.ptrn = self.frmt.format(**patterns)

    def build(self, **vals: Any):
        # Verify that a value is given for every field, with no extras.
        if (v := sorted(vals.keys())) != (f := sorted(self.field_types.keys())):
            raise PathValueError(f"Expected fields {f}, but got {v}")
        # Validate the value passed to every field.
        fields = {name: field.build(vals[name])
                  for name, field in self.field_types.items()}
        # Return the formatted segment.
        return self.frmt.format(**fields)

    def parse(self, text: str):
        ext_val = None
        if (ext := self.field_types.get(EXT)) is not None:
            # If the segment has a file extension, then determine the
            # extension and remove it from the end of the text.
            for option in ext.options:
                if text.endswith(option):
                    # If the text ends with one of the extensions, use
                    # that extension. If it ends with more than one,
                    # use the longest matching extension.
                    if ext_val is None or len(option) > len(ext_val):
                        ext_val = option
            # Check whether the text ends in a valid file extension.
            if ext_val is None:
                raise PathValueError(f"Segment '{text}' lacks valid extension; "
                                     f"options are as follows: {ext.options}")
            # Remove the file extension from the end of the text.
            text = text[:-len(ext_val)]
        # Try to parse the text (with the extension, if any, removed).
        if not (match := re.match(self.ptrn, text)):
            raise PathValueError(f"Could not parse fields in text '{text}'"
                                 f"using pattern '{self.ptrn}'")
        vals = list(match.groups())
        # If there is an extension field, add its value back to the end
        # of the parsed values.
        if ext_val is not None:
            vals.append(ext_val)
        # Return a dict of the names of the fields in the segment and
        # their parsed values.
        fields = {name: field.parse(group) for (name, field), group
                  in zip(self.field_types.items(), vals, strict=True)}
        return fields

    def __str__(self):
        return f"Segment({', '.join(self.field_types)})"


# Field names
TOP = "top"

STEP = "step"
MOD = "module"
SAMP = "sample"
REF = "ref"
END5 = "end5"
END3 = "end3"
BATCH = "batch"
EXT = "ext"

# Directory segments
TopSeg = Segment({TOP: TopField})
ModSeg = Segment({MOD: ModField})
StepSeg = Segment({STEP: StepField})
SampSeg = Segment({SAMP: NameField})
RefSeg = Segment({REF: NameField})
SectSeg = Segment({END5: IntField, END3: IntField}, frmt="{end5}-{end3}")

# File segments
FastaSeg = Segment({REF: NameField, EXT: FastaExt})
FastaIndexSeg = Segment({REF: NameField, EXT: FastaIndexExt})
FastqSeg = Segment({SAMP: NameField, EXT: FastqExt})
Fastq1Seg = Segment({SAMP: NameField, EXT: Fastq1Ext})
Fastq2Seg = Segment({SAMP: NameField, EXT: Fastq2Ext})
DmFastqSeg = Segment({REF: NameField, EXT: FastqExt})
DmFastq1Seg = Segment({REF: NameField, EXT: Fastq1Ext})
DmFastq2Seg = Segment({REF: NameField, EXT: Fastq2Ext})
XamSeg = Segment({REF: NameField, EXT: XamExt})
BamIndexSeg = Segment({REF: NameField, EXT: BamIndexExt})
VecBatSeg = Segment({BATCH: IntField, EXT: VecBatExt})
VecRepSeg = Segment({EXT: VecRepExt}, frmt="report{ext}")


class Path(object):
    def __init__(self, /, *seg_types: Segment):
        self.seg_types = list(seg_types)
        if TopSeg in self.seg_types:
            raise PathValueError(f"{TopSeg} may not be given in seg_types")
        self.seg_types.insert(0, TopSeg)
        for i, seg_type in enumerate(self.seg_types, start=1):
            if EXT in seg_type.field_types and i != len(self.seg_types):
                raise PathValueError("Only the last segment can have a field "
                                     f"with an extension, but segment {i} does")

    def build(self, **fields: Any):
        """ Return a ```pathlib.Path``` instance by assembling the given
        ```fields``` into a full path. """
        # Copy the fields to avoid modifying the original argument.
        fields_left = fields.copy()
        # Build the new path one segment at a time.
        segments = list()
        for seg_type in self.seg_types:
            # For each type of segment in the path, try to get the names
            # and values of all fields of the segment.
            try:
                seg_fields = {name: fields_left.pop(name)
                              for name in seg_type.field_types}
            except KeyError as error:
                raise PathValueError(f"Missing field: {error}")
            # Generate a string representation of the segment using the
            # values of its fields, and add it to the growing path.
            segments.append(seg_type.build(**seg_fields))
        # Check whether any fields were given but not used by the path.
        if fields_left:
            raise PathValueError(f"Unexpected fields: {fields_left}")
        # Assemble the segment strings into a path, and return it.
        return pl.Path(*segments)

    def parse(self, path: str | pl.Path):
        """ Return the field names and values from a given path. """
        # Convert the given path into a canonical, absolute path.
        path = str(sanitize(path))
        # Get the field names and values one segment at a time.
        fields = dict()
        # Iterate from the deepest (last) to shallowest (first) segment.
        for seg_type in reversed(self.seg_types):
            if seg_type is TopSeg:
                # The top-most segment has been reached and must be used
                # to parse the entire remaining path.
                tail = path
            else:
                # The top-most segment of the path has not been reached.
                # Split off the deepest part of the path (tail), and
                # parse it using the current segment type.
                path, tail = os.path.split(path)
            # Verify that the entire path has not been consumed.
            if not tail:
                raise PathValueError(f"No path remains to parse {seg_type}")
            # Parse the deepest part of the path to obtain the fields,
            # and use them to update the field names and values.
            fields.update(seg_type.parse(tail))
        return fields


@cache
def create_path_type(*segment_types: Segment):
    """ Create and cache a Path instance from the segment types. """
    return Path(*segment_types)


def build(*segment_types: Segment, **field_values: Any):
    """ Return a ```pathlib.Path``` from the given segment types and
    field values. """
    return create_path_type(*segment_types).build(**field_values)


def parse(path: str | pl.Path, /, *segment_types: Segment):
    """ Return the fields of a path given as a ```str``` based on the
    segment types. """
    return create_path_type(*segment_types).parse(path)


def find_files(path: pl.Path, segments: Sequence[Segment]) -> list[Path]:
    """ List of all files that match a given sequence of path segments.
    The behavior depends on what ```path``` is:

    - If it is a file, then return a 1-item list containing ```path```
      if it matches the segments, otherwise an empty list.
    - If it is a directory, then search it recursively and return a
      (possibly empty) list of all files in the directory and its
      subdirectories that match the segments.
    - If it does not exist, then raise ```FileNotFoundError```.
    """
    if path.is_file():
        # Check if the given path is a file.
        try:
            # Determine if the file matches the segments.
            parse(path, *segments)
        except PathError:
            # If not, skip it.
            logger.debug(f"File {path} does not match pattern {segments}")
            return []
        else:
            # If so, return it.
            logger.debug(f"File {path} matches pattern {segments}")
            return [path]
    # Otherwise, assume it is a directory and search it for reports.
    logger.debug(f"Searching {path} for files matching {segments}")
    return list(chain(*map(find_files, path.iterdir())))


def find_files_chain(paths: Iterable[pl.Path], segments: Sequence[Segment]):
    """ Call ```find``` on every path in ```paths``` and return a flat
    list of all files matching the segments. """
    found: list[Path] = list()
    for path, count in Counter(paths).items():
        if count > 1:
            logger.warning(f"Path {path} given {count} times")
        try:
            len_init = len(found)
            found.extend(find_files(path, *segments))
            if len(found) == len_init:
                logger.warning(f"Found no files matching {segments} in {path}")
        except FileNotFoundError:
            logger.critical(f"Path does not exist: {path}")
        except Exception as error:
            logger.critical(f"Failed to search for files in {path}: {error}")
    return found
