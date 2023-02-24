from __future__ import annotations

import logging
import sys
from collections import defaultdict, namedtuple
from functools import cache, cached_property
from itertools import repeat, starmap
import os
import re
from datetime import datetime
from hashlib import md5
from multiprocessing import Pool
import time
from typing import Any, ClassVar, Sequence

import numpy as np
import pandas as pd
from pydantic import (BaseModel, Extra, Field, NonNegativeInt, NonNegativeFloat,
                      PositiveInt, StrictStr, validator, root_validator)

from ..util import path
from ..util.seq import (BLANK_INT, MATCH_INT, DELET_INT, INS_5_INT, INS_3_INT,
                        SUB_A_INT, SUB_C_INT, SUB_G_INT, SUB_T_INT, AMBIG_INT,
                        DNA, FastaParser)
from ..util.util import get_num_parallel
from ..vector.samread import SamReader
from ..vector.vector import SamRecord


RegionTuple = namedtuple("PrimerTuple", ["pos5", "pos3"])


class Region(object):
    """
    Represent a region of a reference sequence between two coordinates.

    Attributes
    ----------
    ref_seq: DNA
        Entire reference sequence (not of just the region).
    ref: str
        Name of the reference sequence.
    end5: int (1 ≤ end5 ≤ end3)
        Coordinate of the reference sequence at which the region's 5' end is
        located (1-indexed).
    end3: int (end5 ≤ end3 ≤ len(ref_seq))
        Coordinate of the reference sequence at which the region's 3' end is
        located (1-indexed; end3 itself is included in the region).

    Examples
    --------
    >>> seq = DNA(b"CATCTGGA")
    >>> name = "example"
    >>> region = Region(ref_seq=seq, ref=name, end5=1, end3=8)
    >>> assert region.region_seq == seq
    >>> region = Region(ref_seq=seq, ref=name, end5=1, end3=-1)
    >>> assert region.region_seq == seq
    >>> region = Region(ref_seq=seq, ref=name, end5=-8, end3=8)
    >>> assert region.region_seq == seq
    >>> region = Region(ref_seq=seq, ref=name, end5=3, end3=7)
    >>> assert region.region_seq == DNA(b"TCTGG")
    >>> region = Region(ref_seq=seq, ref=name, end5=-5, end3=-3)
    >>> assert region.region_seq == DNA(b"CTG")
    >>> try:
    ...     region = Region(ref_seq=seq, ref=name, end5=-9, end3=5)
    ...     assert False, "Failed to catch end5 < -len(ref_seq)"
    ... except ValueError:
    ...     pass
    >>> try:
    ...     region = Region(ref_seq=seq, ref=name, end5=6, end3=5)
    ...     assert False, "Failed to catch end3 < end5"
    ... except ValueError:
    ...     pass
    >>> try:
    ...     region = Region(ref_seq=seq, ref=name, end5=1, end3=9)
    ...     assert False, "Failed to catch end3 > len(ref_seq)"
    ... except ValueError:
    ...     pass
    """

    def __init__(self, *, ref_seq: DNA, ref: str, end5: int, end3: int):
        """
        Parameters
        ----------
        ref_seq: DNA
            Entire reference sequence (not of just the region).
        ref: str
            Name of the reference sequence.
        end5: int (-len(ref_seq) ≤ end5 ≤ len(ref_seq); end5 ≠ 0)
            Coordinate of the reference sequence at which the 5' end of the
            region is located. If positive, number the coordinates of the
            reference sequence 1, 2, ... starting at the 5' end (i.e. 1-based
            indexing). If negative, number the coordinates of the reference
            sequence -1, -2, ... starting at the 3' end (i.e. 1-based indexing
            from the other side), then convert to the corresponding (positive)
            1-based index from the 5' end.
        end3: int (-len(ref_seq) ≤ end3 ≤ len(ref_seq); end3 ≠ 0)
            Coordinate of the reference sequence at which the region's 3' end is
            located. Follows the same coordinate numbering convention as end5.
        """
        length = len(ref_seq)
        if end5 < 0:
            # If end5 is negative, find the corresponding positive coordinate.
            end5 += length + 1
        if end3 < 0:
            # If end3 is negative, find the corresponding positive coordinate.
            end3 += length + 1
        if not 1 <= end5 <= end3 <= length:
            raise ValueError("Must have 1 ≤ end5 ≤ end3 ≤ len(ref_seq), "
                             f"but got end5 = {end5}, end3 = {end3}, and "
                             f"len(ref_seq) = {length}")
        self.ref_seq = ref_seq
        self.ref = ref
        self.end5 = end5
        self.end3 = end3

    @property
    def path_fields(self):
        return {"ref": self.ref, "end5": self.end5, "end3": self.end3}

    @property
    def prof_fields(self):
        return {"ref": self.ref, "end5": self.end5, "end3": self.end3}

    @property
    def spanning(self) -> bool:
        """ Return whether the region spans the entire reference sequence. """
        return self.end5 == 1 and self.end3 == len(self.ref_seq)

    @property
    def region_seq(self) -> DNA:
        """ Return the sequence of the region of interest. """
        return self.ref_seq[self.end5 - 1: self.end3]

    @property
    def length(self):
        """ Return the length of the region of interest. """
        return self.end3 - self.end5 + 1

    @cached_property
    def positions(self):
        """ Return all positions in the region of interest as an NDArray. """
        return np.arange(self.end5, self.end3 + 1)

    @property
    def ref_coords(self):
        """ Return the name of the reference and the first and last positions
        of the region of interest; for equality testing and hashing. """
        return self.ref, self.end5, self.end3

    @staticmethod
    def pos_to_cols(positions: Sequence[int]):
        """ Convert positions to column names. """
        return list(map(str, positions))

    @staticmethod
    def cols_to_pos(columns: list[str]):
        """ Convert column names to positions. """
        return np.array(list(map(int, columns)))

    @cached_property
    def columns(self):
        """ Return the column names of the region. """
        return self.pos_to_cols(self.positions)

    @property
    def tag(self):
        """ Return a hashable identifier for the region. """
        return self.ref, self.end5, self.end3

    def __str__(self):
        return f"{self.ref}:{self.end5}-{self.end3}"


class RegionFinder(Region):
    """
    The 5' and 3' ends of a region can be given explicitly as integers, but if
    the sample is of an amplicon (i.e. generated by RT-PCR using site-specific
    primers), then it is often more convenient to enter the sequences of the
    PCR primers and have the software determine the coordinates. RegionFinder
    accepts 5' and 3' coordinates given as integers or primers, validates them,
    and stores the coordinates as integers. This process works as follows:

    end5
    - If given as a parameter, this value is used.
    - Else if fwd is given, first is based on the location of its 3' end.
      - Else first is set to 1.
        - last
          - If last is given as an argument, last is set to this value.
          - Else if rev is given, last is based on the location of the 5' end
            of its reverse complement.
          - Else last is set to the length of ref_seq.

    """

    def __init__(self, *, ref_seq: DNA | None, ref: str, primer_gap: int,
                 end5: int | None = None, end3: int | None = None,
                 fwd: DNA | None = None, rev: DNA | None = None):
        """
        Parameters
        ----------
        ref_seq: see superclass
        ref: see superclass
        primer_gap: int
            (For coordinates specified by fwd/rev only) Number of
            positions 3' of the forward primer and 5' of the reverse
            primer to exclude from the region. Coordinates within 1 - 2
            nucleotides of each primer may contain DMS reactivity
            artifacts. If primer_gap = 0, then end5 and end3 are set,
            respectively, to the coordinates immediately adjacent to
            (i.e. 1 nucleotide 3' and 5' of) the 3' end of the forward
            and reverse primers.
        end5: int | None (default: None)
            If given, behaves as in the superclass; otherwise, ignored.
        end3: int | None (default: None)
            If given, behaves as in the superclass; otherwise, ignored.
        fwd: DNA | None = None (default: None)
            (For amplicons only) Sequence of the forward PCR primer
            that was used to generate the amplicon
        rev: DNA | None = None (default: None)
            (For amplicons only) Sequence of the reverse PCR primer
            that was used to generate the amplicon (the actual sequence,
            not its reverse complement)
        """
        if primer_gap < 0:
            logging.warning("Primer gap must be ≥ 0: setting to 0")
            primer_gap = 0
        if ref_seq is None:
            raise ValueError(f"No sequence for reference named '{ref}'")
        offset = primer_gap + 1
        if end5 is None:
            # If first is to be determined from the fwd primer sequence,
            # the primer is aligned to the reference, and first is set to
            # the position (primer_gap + 1) downstream of its 3' coordinate.
            end5 = (1 if fwd is None
                    else self.locate(ref_seq, fwd).pos3 + offset)
        if end3 is None:
            # If last is to be determined from the rev primer sequence,
            # the reverse complement of the primer is aligned to the reference,
            # and last is set to the position (primer_gap + 1) upstream of its
            # 5' coordinate.
            end3 = (len(ref_seq) if rev is None
                    else self.locate(ref_seq, rev.rc).pos5 - offset)
        super().__init__(ref_seq=ref_seq, ref=ref, end5=end5, end3=end3)

    @staticmethod
    def locate(ref_seq: DNA, primer: DNA) -> RegionTuple:
        """
        Return the 5' and 3' positions (1-indexed) of a primer within a
        reference sequence. The primer must occur exactly once in the
        reference, otherwise an error is raised.

        Parameters
        ----------
        ref_seq: DNA
            Sequence of the entire reference (not just the region of interest)
        primer: DNA
            Sequence of the forward PCR primer or of the reverse complement of
            the reverse PCR primer
        
        Returns
        -------
        RegionTuple
            Named tuple of the first and last positions that the primer occupies
            in the reference sequence. Positions are 1-indexed and include the
            first and last coordinates.
        """
        matches = list(re.finditer(primer, ref_seq))
        if not matches:
            raise ValueError(f"Primer '{primer}' is not in ref '{ref_seq}'")
        if len(matches) > 1:
            raise ValueError(f"Primer '{primer}' occurs {len(matches)} times "
                             f"in ref '{ref_seq}'")
        # Add 1 to convert from 0-indexed (re.finditer) to 1-indexed (DREEM).
        pos5 = matches[0].start() + 1
        # No change is needed to convert from exclusive 0-indexed (re.finditer)
        # to inclusive 1-indexed (DREEM).
        pos3 = matches[0].end()
        return RegionTuple(pos5, pos3)


class MutationalProfile(Region):
    """
    Represent all reads from one sample that overlap a region of interest in a
    reference sequence.

    Fields
    ------
    sample: str
        Name of the sample
    """

    def __init__(self, *, sample: str, **kwargs):
        super().__init__(**kwargs)
        self.sample = sample

    @property
    def path_fields(self):
        return {"sample": self.sample, **super().path_fields}

    @property
    def prof_fields(self):
        return {"sample": self.sample, **super().prof_fields}

    @property
    def tag(self):
        return tuple([self.sample, *super().tag])

    def __str__(self):
        return f"{self.sample}:{super().__str__()}"


class VectorIO(MutationalProfile):
    """
    Read and write mutation vectors.

    Fields
    ------
    top_dir: str
        Path to the top-level directory in which all temporary and final
        output files will be written
    """

    def __init__(self, *, out_dir: path.TopDirPath, **kwargs):
        """
        Initialize a VectorIO object to read and write mutation vectors.
        
        Parameters
        ----------
        out_dir: str
        """
        super().__init__(**kwargs)
        self.out_dir = out_dir
        self.num_batches: int = 0
        self.num_vectors: int = 0
        self.checksums: list[str] = list()

    @property
    def path_fields(self):
        return {"top": self.out_dir.top,
                "module": path.Module.VECTOR,
                **super().path_fields}

    @property
    def prof_fields(self):
        return {"out_dir": self.out_dir,
                "num_batches": self.num_batches,
                "num_vectors": self.num_vectors,
                "checksums": self.checksums,
                **super().prof_fields}

    @property
    def report_path(self):
        return path.MutVectorReportFilePath(**self.path_fields, ext=".json")

    @property
    def batch_dir(self):
        return path.RegionOutDirPath(**self.path_fields)

    def get_mv_batch_path(self, batch: int):
        return path.MutVectorBatchFilePath(**self.path_fields,
                                           batch=batch,
                                           ext=".orc")

    @property
    def batch_nums(self):
        """ Return a range of all batch numbers. """
        return range(self.num_batches)

    @property
    def mv_batch_paths(self):
        """ Return the path of every mutation vector batch file. """
        return list(map(self.get_mv_batch_path, self.batch_nums))

    @classmethod
    def digest_file(cls, file_path: Any) -> str:
        """
        Compute the checksum of a file.
        
        Parameters
        ----------
        file_path: Any (path-like)
            Path of the file on which to compute the checksum. Can be
            any type that the open() function recognizes as a path.
        
        Returns
        -------
        str
            Checksum of the file (in hexadecimal)
        """
        with open(file_path, "rb") as f:
            digest = md5(f.read()).hexdigest()
        return digest


class VectorReport(BaseModel, VectorIO):
    """
    Read and write a report about a mutational profile, with information such as
    - the sample, reference, and region
    - number of mutation vectors
    - number of mutation vector batch files and their checksums
    - beginning and ending time, duration, and speed of vectoring

    Fields
    ------
    out_str: str
        The top-level directory of the output files, as a string; writing to and
        reading from JSON format does not work with TopDirPath objects directly.
    sample: str
        Name of the sample
    

    Examples
    --------
    >>> report = VectorReport(out_str=os.getcwd(), sample="dms2",
    ...                       ref="tmv-rna", end5=1, end3=20,
    ...                       ref_str=DNA(b"GTATTTTTACAACAATTACC"),
    ...                       num_vectors=10346, num_batches=2,
    ...                       checksums=["b47260fcad8be60470bee67862f187b4",
    ...                                  "098f40cfadc266ea5bc48ab2e18cdc95"],
    ...                       began=datetime.now(),
    ...                       ended=(time.sleep(1E-5), datetime.now())[-1])
    >>> report.ref_str
    'GTATTTTTACAACAATTACC'
    """

    class Config:
        allow_population_by_field_name = True
        extra = Extra.ignore

    # Fields
    out_str: StrictStr = Field(alias="Top-level output directory")
    sample: StrictStr = Field(alias="Sample name")
    ref: StrictStr = Field(alias="Reference name")
    end5: PositiveInt = Field(alias="5' end of region")
    end3: PositiveInt = Field(alias="3' end of region")
    ref_str: StrictStr = Field(alias="Reference sequence")
    num_vectors: NonNegativeInt = Field(alias="Number of vectors")
    num_batches: NonNegativeInt = Field(alias="Number of batches")
    checksums: list[StrictStr] = Field(alias="MD5 checksums of vector batches")
    began: datetime = Field(alias="Began vectoring")
    ended: datetime = Field(alias="Ended vectoring")
    duration: NonNegativeFloat = Field(default=float("nan"),
                                       alias="Duration of vectoring (s)")
    speed: NonNegativeFloat = Field(default=float("nan"),
                                    alias="Speed of vectoring (vectors/s)")

    # Format of dates and times in the report file
    dt_fmt: ClassVar[str] = "on %Y-%m-%d at %H:%M:%S.%f"

    @validator("out_str", pre=True)
    def convert_out_dir_to_str(cls, out_str: str | path.TopDirPath):
        """ Return top-level directory (TopDirPath) as a string.
        Must be str in order to write and load from JSON correctly. """
        if isinstance(out_str, str):
            return out_str
        if isinstance(out_str, path.TopDirPath):
            return out_str.top
        raise TypeError(out_str)

    @property
    def out_dir(self):
        """ Return top-level directory string (from JSON) as a TopDirPath.
        The methods from VectorIO expect top_dir to be of type TopDirPath. """
        return path.TopDirPath.parse(self.out_str)

    @validator("ref_str", pre=True)
    def convert_ref_seq_to_str(cls, ref_str: DNA):
        """ Return reference sequence (DNA) as a string.
        Must be str in order to write to and load from JSON correctly. """
        return str(ref_str)

    @property
    def ref_seq(self):
        """ Return reference sequence string (from JSON) as a DNA object.
        The methods from VectorIO expect ref_seq to be of type DNA. """
        return DNA(self.ref_str.encode())

    @root_validator(pre=False)
    def calculate_duration_and_speed(cls, values):
        """
        Calculate and return the duration and speed of vectoring:
        - duration: difference between ending and beginning time (sec)
        - speed: number of vectors processed per unit time (vectors/sec)
        """
        began = values["began"]
        ended = values["ended"]
        dt = ended - began
        # Convert seconds and microseconds (both int) to seconds (float)
        duration = dt.seconds + dt.microseconds / 1E6
        values["duration"] = duration
        if duration < 0.0:
            # Duration may not be negative.
            logging.error(f"Began at {began.strftime(cls.dt_fmt)}, but ended "
                          f"earlier, at {ended.strftime(cls.dt_fmt)}: "
                          "setting duration to 0 sec.")
            duration = 0.0
        num_vectors = values["num_vectors"]
        # Handle the unlikely case that duration == 0.0 by returning
        # - inf (i.e. 1 / 0) if at least 1 vector was processed
        # - nan (i.e. 0 / 0) if no vectors were processed
        # when calculating the speed of processing vectors.
        if duration == 0.0:
            logging.warning("Cannot compute speed because duration is 0 sec.")
            speed = float("inf" if num_vectors else "nan")
        else:
            speed = round(num_vectors / duration, 1)
        values["speed"] = speed
        return values

    @root_validator(pre=False)
    def num_batches_len_checksums(cls, values):
        num_batches = values["num_batches"]
        num_checksums = len(values["checksums"])
        if num_batches != num_checksums:
            raise ValueError(f"Numbers of batches ({num_batches}) and "
                             f"checksums ({num_checksums}) did not match.")
        return values

    @root_validator(pre=False)
    def region_is_valid(cls, values):
        # The initialization of this Region instance will raise an error if the
        # region is not valid.
        Region(ref_seq=values["ref_str"], ref=values["ref"],
               end5=values["end5"], end3=values["end3"])
        return values

    def find_invalid_batches(self, validate_checksums: bool):
        """ Return all the batches of mutation vectors that either do not exist
        or do not match their expected checksums. """
        missing = list()
        badsum = list()
        for file, checksum in zip(self.mv_batch_paths, self.checksums,
                                  strict=True):
            fpath = file.path
            if fpath.is_file():
                if validate_checksums and self.digest_file(fpath) != checksum:
                    # The batch file exists but does not match the checksum.
                    badsum.append(file)
            else:
                # The batch file does not exist.
                missing.append(file)
        return missing, badsum

    def assert_valid_batches(self, validate_checksums: bool):
        missing, badsum = self.find_invalid_batches(validate_checksums)
        if missing:
            raise FileNotFoundError(f"Missing vector batch files: {missing}")
        if badsum:
            raise ValueError(f"Batch files have bad checksums: {badsum}")

    def get_reader(self):
        return VectorReader.from_report(self)

    def save(self):
        text = self.json(by_alias=True)
        self.report_path.path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.report_path.path, "w") as f:
            f.write(text)

    @classmethod
    def load(cls, file: Any, validate_checksums: bool = True):
        """ Load a mutation vector report from a file. """
        file_path = path.MutVectorReportFilePath.parse(file)
        report = cls.parse_file(file_path.path)
        if (file_path.top != report.out_dir.top
                or file_path.sample != report.sample
                or file_path.ref != report.ref
                or file_path.end5 != report.end5
                or file_path.end3 != report.end3):
            raise ValueError(f"Report fields do not match path '{file_path}'")
        report.assert_valid_batches(validate_checksums)
        return report


class VectorWriter(VectorIO):
    """
    Compute mutation vectors for all reads from one sample mapping to one
    region of one reference sequence.
    """

    def __init__(self, *,
                 bam_path: path.OneRefAlignmentInFilePath,
                 batch_size: int,
                 min_qual: int,
                 rerun: bool,
                 **kwargs):
        super().__init__(sample=bam_path.sample,
                         ref=bam_path.ref,
                         **kwargs)
        self.bam_path = bam_path
        self.batch_size = batch_size
        self.min_qual = min_qual
        self.seq_bytes = bytes(self.region_seq)
        self.rerun = rerun
        self._vectorized = False

    def _write_report(self, began: datetime, ended: datetime):
        report = VectorReport(**{**self.prof_fields,
                                 "out_str": self.out_dir,
                                 "ref_str": self.ref_seq,
                                 "began": began,
                                 "ended": ended})
        report.save()
        return report

    def _write_batch(self,
                     read_names: list[str],
                     muts: tuple[bytearray, ...],
                     batch_num: int) -> tuple[path.MutVectorBatchFilePath, int]:
        """
        Write a batch of mutation vectors to an ORC file.

        ** Arguments **
        muts (NDArray) --> batch of mutation vectors in which each row is a
                           mutation vector and each column is a position in the
                           region of interest
        batch_num (int) -> non-negative integer label for the batch; every
                           mutational profile containing n batches includes all
                           batch numbers i in the range 0 <= i < n

        ** Returns **
        mv_file (str) <--- file path where the mutation vectors were written
        """
        # Convert the concatenated mutation vectors to a 2D NumPy array.
        muts_array = np.frombuffer(b"".join(muts), dtype=np.byte)
        n_records = len(muts_array) // self.length
        muts_array.resize((n_records, self.length))
        # Data must be converted to pd.DataFrame for PyArrow to write.
        # Explicitly set copy=False to copying the mutation vectors.
        df = pd.DataFrame(data=muts_array, index=read_names,
                          columns=self.columns, copy=False)
        mv_file = self.get_mv_batch_path(batch_num)
        df.to_orc(mv_file.path, index=True, engine="pyarrow")
        return mv_file, n_records

    def _vectorize_record(self, rec: SamRecord):
        """
        Compute the mutation vector of one record from a SAM file.

        ** Arguments **
        rec (SamRecord) --> SAM record for which to compute a mutation vector

        ** Returns **
        muts (bytearray) <- mutation vector
        """
        try:
            if rec.ref != self.ref:
                raise ValueError(f"Read '{rec.read_name}' had reference"
                                 f"'{rec.ref}' differing from profile "
                                 f"reference '{self.ref}'")
            muts = rec.vectorize(self.seq_bytes, self.end5, self.end3)
            if not any(muts):
                raise ValueError(f"Vector for read '{rec.read_name}' was blank")
        except ValueError as error:
            logging.error(f"Read '{rec.read1.qname.decode()}' failed to "
                          f"vectorize due to the following error: {error}")
            return "", bytearray()
        return rec.read_name, muts

    def _vectorize_batch(self, reader: SamReader, batch_num: int,
                         start: int, stop: int):
        """
        Generate a batch of mutation vectors and write them to a file.

        Parameters
        ----------
        reader: SamReader
            Reader for the SAM file for which to generate a batch of
            mutation vectors
        batch_num: int (≥ 0)
            Non-negative integer label for the batch
        start: int (≥ 0)
            Seek to this position in the SAM file, then start generating
            vectors; position must be the beginning of a line, else an
            error will be raised.
        stop: int (≥ start)
            Stop generating vectors upon reaching this position in the
            SAM file; position must be the beginning of a line, else an
            error will be raised.
        
        Returns
        -------
        int
            Number of records read from the SAM file between positions
            start and stop
        str
            MD5 checksum of the ORC file of the batch of vectors
        """
        if stop > start:
            with reader as readopen:
                # Use the SAM reader to generate the mutation vectors.
                # Collect them as a single, 1-dimensional bytes object.
                read_names, muts = zip(*map(self._vectorize_record,
                                            readopen.get_records(start, stop)))
                # For every read for which creating a mutation vector
                # failed, an empty string was returned as the read name
                # and an empty bytearray as the mutation vector. The
                # empty read names must be filtered out, while the empty
                # mutation vectors will not cause problems because,
                # being of length zero, they will effectively disappear
                # when all the vectors are concatenated into a 1D array.
                read_names = list(filter(None, read_names))
        else:
            read_names, muts = (), ()
        # Write the mutation vectors to a file and compute its checksum.
        mv_file, n_records = self._write_batch(read_names, muts, batch_num)
        checksum = self.digest_file(mv_file.path)
        return n_records, checksum

    def _vectorize_sam(self, temp_dir: path.TopDirPath, n_procs: int,
                       save_temp: bool, resume: bool):
        if self._vectorized:
            raise RuntimeError(f"Vectoring was already run for {self}.")
        self._vectorized = True
        # Open the primary SAM file reader to write the subset of SAM
        # records to a temporary SAM file and determine the number and
        # start/stop indexes of each batch of records in the file.
        with SamReader(temp_dir=temp_dir,
                       save_temp=save_temp,
                       resume=resume,
                       n_procs=n_procs,
                       xam_path=self.bam_path,
                       ref_name=self.ref,
                       end5=self.end5,
                       end3=self.end3,
                       spanning=self.spanning,
                       min_qual=self.min_qual) as reader:
            vectors_per_batch = max(1, self.batch_size // self.length)
            indexes = list(reader.get_batch_indexes(vectors_per_batch))
            starts = indexes[:-1]
            stops = indexes[1:]
            self.num_batches = len(starts)
            # Once the number of batches has been determined, a list of
            # new SAM file readers is created. Each is responsible for
            # converting one batch of SAM reads into one batch of
            # mutation vectors. Setting owner=False prevents the SAM
            # reader from creating a new SAM file (as well as from
            # deleting it upon exit); instead, it uses the file written
            # by the primary SAM reader, which is reader.sam_path.
            readers = [SamReader(temp_dir=temp_dir,
                                 save_temp=save_temp,
                                 resume=resume,
                                 n_procs=n_procs,
                                 xam_path=reader.sam_path,
                                 ref_name=self.ref,
                                 end5=self.end5,
                                 end3=self.end3,
                                 spanning=self.spanning,
                                 min_qual=self.min_qual,
                                 owner=False)
                       for _ in self.batch_nums]
            # The _vectorize_batch method requires a SamReader, batch
            # number, start index, and stop index; these are zipped into
            # this list of arguments.
            args = list(zip(readers, self.batch_nums, starts, stops,
                            strict=True))
            # Then pass the arguments to _vectorize_batch.
            if (pool_size := min(n_procs, self.num_batches)) > 1:
                with Pool(pool_size) as pool:
                    results = list(pool.starmap(self._vectorize_batch, args,
                                                chunksize=1))
            else:
                results = list(starmap(self._vectorize_batch, args))
            # results is a list containing, for each batch, a tuple of
            # the number of vectors in the batch and the checksum of
            # the batch. Gather them into num_vectors and checksums.
            for num_vectors, checksum in results:
                self.num_vectors += num_vectors
                self.checksums.append(checksum)

    @property
    def outputs_valid(self):
        try:
            VectorReport.load(self.report_path.path)
        except (FileNotFoundError, ValueError):
            return False
        else:
            return True

    def vectorize(self, temp_dir: path.TopDirPath, n_procs: int,
                  save_temp: bool, resume: bool):
        if self.rerun or not self.outputs_valid:
            self.batch_dir.path.mkdir(parents=True, exist_ok=True)
            logging.info(f"{self}: computing vectors")
            began = datetime.now()
            self._vectorize_sam(temp_dir, n_procs, save_temp, resume)
            ended = datetime.now()
            logging.info(f"{self}: writing report")
            self._write_report(began, ended)
            logging.info(f"{self}: finished")
        else:
            logging.warning(f"{self}: already finished. To rerun, add --rerun")
        return self.report_path


class VectorWriterSpawner(object):
    def __init__(self, *,
                 out_dir: path.TopDirPath,
                 temp_dir: path.TopDirPath,
                 refset_path: path.RefsetSeqInFilePath,
                 bam_paths: list[path.OneRefAlignmentInFilePath],
                 coords: list[tuple[str, int, int]],
                 primers: list[tuple[str, DNA, DNA]],
                 primer_gap: int,
                 cfill: bool,
                 batch_size: int,
                 parallel: bool,
                 max_procs: int,
                 min_phred: int,
                 phred_enc: int,
                 rerun: bool):
        self.out_dir = out_dir
        self.temp_dir = path.TopDirPath.parse(temp_dir)
        self.bam_paths = bam_paths
        self.refset_path = path.RefsetSeqInFilePath.parse(refset_path)
        self.coords = coords
        self.primers = primers
        if primer_gap < 0:
            logging.warning("Primer gap must be ≥ 0: setting to 0")
            primer_gap = 0
        self.primer_gap = primer_gap
        self.fill = cfill
        self.batch_size = batch_size
        self.parallel = parallel
        if max_procs < 1:
            logging.warning("Max CPUs must be ≥ 1: setting to 1")
            max_procs = 1
        self.max_procs = max_procs
        self.min_phred = min_phred
        self.phred_enc = phred_enc
        self.rerun = rerun

    @property
    def min_qual(self):
        """
        Return the minimum quality for a base in a read to be considered
        informative, as the ASCII integer corresponding to the character in the
        FASTQ file that is the minimum valid quality.

        Returns
        -------
        int
            The ASCII value corresponding to the character in the FASTQ file read
            quality string that represents the minimum quality for a base to be
            considered informative.

        Examples
        --------
        For example, if the minimum Phred score (```self.min_phred```) that is
        accepted as informative is 20, and the Phred encoding of the FASTQ file
        (```phred_enc```) is 33 (i.e. ASCII+33), then the minimum quality as an
        ASCII integer (```min_qual```) is 20 + 33 = 53, which is character '5'.
        If ```min_phred``` were 37, then ```min_qual``` would be 37 + 33 = 70,
        which is character 'F'.
        """
        return self.min_phred + self.phred_enc

    @cached_property
    def bams_per_sample(self):
        samples: dict[str, list[
            path.OneRefAlignmentInFilePath]] = defaultdict(list)
        for bam in self.bam_paths:
            samples[bam.sample].append(bam)
        return samples

    @property
    def samples(self):
        return set(self.bams_per_sample.keys())

    @cached_property
    def ref_seqs(self):
        seqs = dict(FastaParser(self.refset_path.path).parse())
        if not seqs:
            logging.critical(f"'{self.refset_path}' contained no sequences")
        return seqs

    @cached_property
    def regions(self):
        regions: dict[str, list[RegionFinder]] = defaultdict(list)

        def add_region(region: RegionFinder):
            if any(region == other for other in regions[region.ref]):
                raise ValueError(f"Duplicate region: {region.ref_coords}")
            regions[region.ref].append(region)

        for ref, first, last in self.coords:
            add_region(RegionFinder(ref_seq=self.ref_seqs.get(ref),
                                    ref=ref, end5=first, end3=last,
                                    primer_gap=self.primer_gap))
        for ref, fwd, rev in self.primers:
            add_region(RegionFinder(ref_seq=self.ref_seqs.get(ref),
                                    ref=ref, fwd=fwd, rev=rev,
                                    primer_gap=self.primer_gap))
        if self.fill:
            for ref, seq in self.ref_seqs.items():
                if ref not in regions:
                    add_region(RegionFinder(ref_seq=seq, ref=ref,
                                            primer_gap=self.primer_gap))
        return regions

    @cached_property
    def writers(self):
        writers: dict[tuple, VectorWriter] = dict()
        for bam in self.bam_paths:
            for region in self.regions[bam.ref]:
                if region.ref != bam.ref:
                    logging.error(f"Skipping region {region} of {bam.path} "
                                  "because its reference does not match that "
                                  f"of the BAM file ('{bam.ref}').")
                    continue
                writer = VectorWriter(out_dir=self.out_dir,
                                      bam_path=bam,
                                      ref_seq=self.ref_seqs[bam.ref],
                                      end5=region.end5,
                                      end3=region.end3,
                                      batch_size=self.batch_size,
                                      min_qual=self.min_qual,
                                      rerun=self.rerun)
                if writer.tag in writers:
                    logging.warning("Skipping duplicate mutational profile: "
                                    f"{writer}.")
                    continue
                writers[writer.tag] = writer
        return list(writers.values())

    def generate_profiles(self, save_temp: bool, resume: bool):
        n_profiles = len(self.writers)
        if n_profiles == 0:
            logging.critical("No BAM files and/or regions specified")
            return ()
        # Determine method of parallelization. Do not use the hybrid
        # feature, which would try to process multiple SAM files in
        # parallel and use multiple processes for each file. Python
        # multiprocessing.Pool forbids a daemon process (one for each
        # SAM file in parallel) from spawning additional processes,
        # which would be needed to use multiple processes on each file.
        n_tasks_parallel, n_procs_per_task = get_num_parallel(n_profiles,
                                                              self.max_procs,
                                                              self.parallel)
        # VectorWriter.vectorize takes two arguments, temp_dir and
        # n_procs, which are the same for all vector writers and are
        # thus repeated using itertools.repeat
        args = tuple(zip(self.writers,
                         repeat(self.temp_dir),
                         repeat(n_procs_per_task),
                         repeat(save_temp),
                         repeat(resume)))
        # Call the vectorize method of each writer, passing args.
        if n_tasks_parallel > 1:
            with Pool(n_tasks_parallel) as pool:
                report_files = tuple(pool.starmap(VectorWriter.vectorize, args,
                                                  chunksize=1))
        else:
            report_files = tuple(starmap(VectorWriter.vectorize, args))
        return report_files


class VectorReader(VectorIO):
    INDEX_COL = "__index_level_0__"

    def __init__(self, num_vectors: int, num_batches: int, **kwargs):
        super().__init__(**kwargs)
        self.num_vectors = num_vectors
        self.num_batches = num_batches

    @property
    def shape(self):
        return self.num_vectors, self.length

    def get_batch(self,
                  batch: int,
                  positions: Sequence[int] | None = None):
        """
        Return the mutation vectors from one batch. Optionally, select
        a subset of the columns of the mutation vectors.

        Parameters
        ----------
        batch: int (≥ 0)
            Number of the batch of mutation vectors
        positions: sequence[int] | None (default: None)
            If given, use only these positions from the mutation vectors;
            otherwise, use all positions.

        Return
        ------
        DataFrame
            Mutation vectors; each row is a vector indexed by its name,
            each column a position indexed by its positional number
        """
        if batch not in self.batch_nums:
            raise ValueError(f"Invalid batch number: {batch} "
                             f"(expected one of {list(self.batch_nums)})")
        columns = ([self.INDEX_COL] + self.pos_to_cols(positions) if positions
                   else None)
        vectors = pd.read_orc(self.get_mv_batch_path(batch).path,
                              columns=columns)
        vectors.set_index(self.INDEX_COL, drop=True, inplace=True)
        vectors.columns = positions if positions else self.positions
        return vectors

    def get_all_batches(self, positions: Sequence[int] | None = None):
        """
        Yield every batch of mutation vectors.

        Parameters
        ----------
        positions: sequence[int] | None (default: None)
            If given, use only these position from the mutation vectors;
            otherwise, use all positions.

        Yield
        -----
        DataFrame
            Mutation vectors; each row is a vector indexed by its name,
            each column a position indexed by its base and number
        """
        for batch in self.batch_nums:
            yield self.get_batch(batch, positions)

    def get_all_vectors(self, positions: Sequence[int] | None = None):
        """
        Return all mutation vectors for this vector reader. Note that
        reading all vectors could take more than the available memory
        and cause the program to crash. Thus, use this method only if
        all vectors will fit into memory. Otherwise, use the method
        ```get_all_batches``` to process the vectors in small batches.

        Parameters
        ----------
        positions: sequence[int] | None (default: None)
            If given, use only these position from the mutation vectors;
            otherwise, use all positions.

        Return
        ------
        DataFrame
            Mutation vectors; each row is a vector indexed by its name,
            each column a position indexed by its base and number
        """
        return pd.concat(self.get_all_batches(positions), axis=0)

    @staticmethod
    def _query_vectors(vectors: pd.DataFrame, query: int) -> pd.DataFrame:
        """
        Return a boolean array of the same shape as vectors where
        element i,j is True if and only if the byte at element i,j of
        vectors is both non-zero and a bitwise subset of query (a byte
        equal to query also counts as a subset).

        Parameters
        ----------
        vectors: DataFrame
            Mutation vectors
        query: int
            Byte to query

        Return
        ------
        DataFrame
            Boolean type DataFrame of the same shape as vectors where
            each element is True if the element at the same position in
            vectors matched the query and False otherwise
        """
        # Flag as False all bytes in vectors that have no bits set to 1,
        # since these bytes represent positions in vectors that were not
        # covered by reads and thus should not count as query matches.
        covered = vectors.astype(bool, copy=False)
        if query == AMBIG_INT:
            # If the query byte is all 1s (i.e. decimal 255), then the
            # next step (bitwise OR) will return True for every byte,
            # so the return value will equal that of covered. For the
            # sake of speed, return covered now.
            # Note that for the sake of speed and memory, covered is
            # a view to the SAME DataFrame as vectors. Thus, functions
            # that use covered should merely read it, NEVER modify it.
            return covered
        # Flag as True all bytes in vectors that are subsets of the
        # query byte (including, for now, bytes in vector that are
        # 00000000). In order for a byte in vectors to be a subset of
        # the query byte, every one of its bits that is set to 1 must
        # also be set to 1 in the query byte. Equivalently, the union
        # (bitwise OR) of the vector byte and query byte must equal the
        # query byte; if not, the vector byte would have had at least
        # one bit set to 1 that was not set to 1 in the query byte.
        return covered & np.equal(np.bitwise_or(vectors, query), query)

    @cache
    def count_query(self,
                    query: int,
                    positions: Sequence[int] | None = None) -> pd.Series:
        """
        Return, for each column in the mutational profile, the number of
        vectors that match the query.

        Parameters
        ----------
        query: int (0 ≤ query < 256)
            Query byte: to match, a byte in the vector must be equal to
            or a bitwise subset of the query byte, and must be non-zero.
        positions: sequence[int] | None (default: None)
            If given, use only these positions from the mutation vectors;
            otherwise, use all positions.

        Return
        ------
        Series
            Number of vectors matching the query at each position; index
            is the position in the vector.
        """
        # For each batch of vectors, get a DataFrame of boolean values
        # indicating query matches (True) and mismatches (False), sum
        # over axis 0 to count the number of matches at each position,
        # and sum these counts over all batches to get the total counts.
        counts = pd.Series(np.zeros_like(positions, dtype=int),
                           index=(self.positions if positions is None
                                  else positions))
        for vectors in self.get_all_batches(positions):
            counts += self._query_vectors(vectors, query).sum(axis=0)
        return counts

    @classmethod
    def from_report(cls, report: VectorReport):
        return cls(out_dir=report.out_dir,
                   sample=report.sample,
                   ref_seq=report.ref_seq,
                   ref=report.ref,
                   end5=report.end5,
                   end3=report.end3,
                   num_vectors=report.num_vectors,
                   num_batches=report.num_batches)

    @classmethod
    def from_report_file(cls, report_file):
        return cls.from_report(VectorReport.load(report_file))


class VectorTextTranslator(object):
    table = bytes.maketrans(*map(b"".join, zip(*[
        (i.to_bytes(length=1, byteorder=sys.byteorder),
         (b"." if i == BLANK_INT
          else b"~" if i == MATCH_INT
          else b"/" if i == DELET_INT
          else b"{" if i == (INS_5_INT | MATCH_INT)
          else b"}" if i == (INS_3_INT | MATCH_INT)
          else b"A" if i == SUB_A_INT
          else b"C" if i == SUB_C_INT
          else b"G" if i == SUB_G_INT
          else b"T" if i == SUB_T_INT
          else b"?"))
        for i in range(256)])))

    @classmethod
    def itertrans(cls, vectors: pd.DataFrame):
        for index, row in zip(vectors.index, vectors.values, strict=True):
            translated = row.tobytes(order='C').translate(cls.table)
            yield b"%b\t%b\n" % (index.encode(), translated)

    @classmethod
    def blocktrans(cls, vectors: pd.DataFrame):
        return b"".join(cls.itertrans(vectors))
