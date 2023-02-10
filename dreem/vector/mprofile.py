from __future__ import annotations

import logging
from collections import defaultdict, namedtuple
from functools import cached_property
from itertools import repeat, starmap
import os
import pathlib
import re
from datetime import datetime
from hashlib import file_digest
from multiprocessing import Pool
import time
from typing import Any, ClassVar

import numpy as np
import pandas as pd
from pydantic import (BaseModel, Extra, Field, NonNegativeInt, NonNegativeFloat, PositiveInt,
                      StrictStr, validator, root_validator)

from dreem.util.cli import ParallelChoice
from dreem.util import path
from dreem.util.seq import DNA, FastaParser
from dreem.vector.samview import SamViewer
from dreem.vector.vector import SamRecord


DEFAULT_BATCH_SIZE = 33_554_432  # 2^25 bytes ≈ 33.6 Mb
DEFAULT_MIN_PHRED = 25  # minimum Phred score to consider a base in a FASTQ file


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

    @property
    def positions(self):
        """ Return all positions in the region of interest as an NDArray. """
        return np.arange(self.end5, self.end3 + 1)

    @property
    def ref_coords(self):
        """ Return the name of the reference and the first and last positions
        of the region of interest; for equality testing and hashing. """
        return self.ref, self.end5, self.end3
    
    @property
    def columns(self):
        """ Return a tuple of the bases and coordinates in the region, each of
        the form '{base}{position}' (e.g. ['G13', 'C14', 'A15']). """
        return list(f"{chr(base)}{coord}" for base, coord
                    in zip(self.region_seq, self.positions, strict=True))


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

    def __init__(self, *, ref_seq: DNA | None, ref: str,
                 end5: int | None = None, end3: int | None = None,
                 fwd: DNA | None = None, rev: DNA | None = None,
                 primer_gap: int = 2):
        """
        Parameters
        ----------
        ref_seq: see superclass
        ref: see superclass
        end5: see superclass
        end3: see superclass
        fwd: DNA | None = None (optional)
            (For amplicons only) Sequence of the forward PCR primer that was
            used to generate the amplicon
        rev: DNA | None = None (optional)
            (For amplicons only) Sequence of the reverse PCR primer that was
            used to generate the amplicon
        primer_gap: int = 2 (optional)
            (For coordinates specified by fwd/rev only) Number of positions 3'
            of the forward primer and 5' of the reverse primer to exclude from
            the region. Coordinates within 1 - 2 nucleotides of each primer may
            contain DMS reactivity artifacts. If primer_gap = 0, then end5 and
            end3 are set, respectively, to the coordinates immediately adjacent
            to (i.e. 1 nucleotide 3' and 5' of) the 3' end of the forward and
            reverse primers. The default (primer_gap = 2) generally works well.
        """
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

    def __str__(self):
        return (f"Mutational Profile for sample '{self.sample}' "
                f"over reference '{self.ref}' region {self.end5}-{self.end3}")


class VectorIO(MutationalProfile):
    """
    Read and write mutation vectors.

    Fields
    ------
    top_dir: str
        Path to the top-level directory in which all temporary and final
        output files will be written
    """

    # Hashing algorithm to compute checksums of mutation vector files.
    digest_algo: ClassVar = "md5"

    def __init__(self, *, top_dir: path.TopDirPath, **kwargs):
        """
        Initialize a VectorIO object to read and write mutation vectors.
        
        Parameters
        ----------
        top_dir: str
        """
        super().__init__(**kwargs)
        self.top_dir = top_dir
        self.num_batches: int = 0
        self.num_vectors: int = 0
        self.checksums: list[str] = list()

    @property
    def path_fields(self):
        return {"top": self.top_dir.top, "partition": path.Partition.OUTPUT,
                "module": path.Module.VECTOR, **super().path_fields}

    @property
    def prof_fields(self):
        return {"top_dir": self.top_dir, "num_batches": self.num_batches,
                "num_vectors": self.num_vectors, "checksums": self.checksums,
                **super().prof_fields}

    @property
    def report_path(self):
        return path.MutVectorReportFilePath(**self.path_fields, ext=".json")

    @property
    def batch_dir(self):
        return path.RegionOutDirPath(**self.path_fields)

    def get_mv_batch_path(self, batch: int):
        return path.MutVectorBatchFilePath(**self.path_fields,
                                           batch=batch, ext=".orc")
    
    @property
    def batch_nums(self):
        """ List all the batch numbers. """
        return range(self.num_batches)

    @property
    def mv_batch_paths(self):
        return list(map(self.get_mv_batch_path, self.batch_nums))

    @classmethod
    def digest_file(cls, file_path: pathlib.Path) -> str:
        """
        Compute the checksum of a file.
        
        Parameters
        ----------
        file_path: Any
            Path of the file on which to compute the checksum. Can be any type
            as long as casting it to a string yields a valid path.
        
        Returns
        -------
        str
            Checksum of the file (in hexadecimal)
        """
        with open(file_path, "rb") as f:
            digest = file_digest(f, cls.digest_algo).hexdigest()
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
    top_str: str
        The top-level directory of the output files, as a string; writing to and
        reading from JSON format does not work with TopDirPath objects directly.
    sample: str
        Name of the sample
    

    Examples
    --------
    >>> report = VectorReport(top_str=os.getcwd(), sample="dms2",
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
    top_str: StrictStr = Field(alias="Top-level output directory")
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

    @validator("top_str", pre=True)
    def convert_top_dir_to_str(cls, top_str: str | path.TopDirPath):
        """ Return top-level directory (TopDirPath) as a string.
        Must be str in order to write and load from JSON correctly. """
        if isinstance(top_str, str):
            return top_str
        if isinstance(top_str, path.TopDirPath):
            return top_str.top
        raise TypeError(top_str)

    @property
    def top_dir(self):
        """ Return top-level directory string (from JSON) as a TopDirPath.
        The methods from VectorIO expect top_dir to be of type TopDirPath. """
        return path.TopDirPath.parse_path(self.top_str)

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
        - duration: difference between ending and beginning time (seconds)
        - speed: number of vectors processed per unit time (vectors/second)
        """
        began = values["began"]
        ended = values["ended"]
        dt = ended - began
        # Time delta objects store time as seconds and microseconds (both int).
        # These values must be converted to a float number of seconds like so:
        duration = dt.seconds + dt.microseconds / 1E6
        values["duration"] = duration
        if duration < 0:
            # Duration may not be negative.
            raise ValueError(f"Began {ended.strftime(cls.dt_fmt)}, but ended "
                             f"earlier, {began.strftime(cls.dt_fmt)}")
        num_vectors = values["num_vectors"]
        # Handle the unlikely case that duration == 0.0 by returning either
        # - inf (i.e. 1 / 0) if at least 1 vector was processed
        # - nan (i.e. 0 / 0) if no vectors were processed
        # when calculating the speed of processing vectors.
        speed = (float("inf" if num_vectors else "nan") if duration == 0.0
                 else round(num_vectors / duration, 1))
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

    def find_invalid_batches(self):
        """ Return all the batches of mutation vectors that either do not exist
        or do not match their expected checksums. """
        missing = list()
        badsum = list()
        for file, checksum in zip(self.mv_batch_paths, self.checksums,
                                  strict=True):
            try:
                if self.digest_file(file.pathstr) != checksum:
                    # The batch file exists but does not match the checksum.
                    badsum.append(file.pathstr)
            except FileNotFoundError:
                # The batch file does not exist.
                missing.append(file.pathstr)
        return missing, badsum

    def assert_valid_batches(self):
        missing, badsum = self.find_invalid_batches()
        if missing:
            raise FileNotFoundError(f"Missing vector batch files: {missing}")
        if badsum:
            raise ValueError(f"Batch files have bad checksums: {badsum}")

    def save(self):
        text = self.json(by_alias=True)
        self.report_path.path.parent.mkdir(parents=True, exist_ok=True)
        with open(self.report_path.path, "w") as f:
            f.write(text)

    @classmethod
    def load(cls, file: Any):
        """ Load a mutation vector report from a file. """
        file_path = path.MutVectorReportFilePath.parse_path(file)
        report = cls.parse_file(file_path.path)
        if (file_path.top != report.top_dir.top
                or file_path.sample != report.sample
                or file_path.ref != report.ref
                or file_path.end5 != report.end5
                or file_path.end3 != report.end3):
            raise ValueError(f"Report fields do not match path '{file_path}'")
        report.assert_valid_batches()
        return report


class VectorWriter(VectorIO):
    """
    Compute mutation vectors for all reads from one sample mapping to one
    region of one reference sequence.
    """
    def __init__(self, *,
                 bam_path: path.OneRefAlignmentInFilePath,
                 min_qual: int,
                 rerun: bool,
                 **kwargs):
        super().__init__(sample=bam_path.sample,
                         ref=bam_path.ref,
                         **kwargs)
        self.bam_path = bam_path
        self.min_qual = min_qual
        self.seq_bytes = bytes(self.region_seq)
        self.rerun = rerun
        self._vectorized = False

    def _write_report(self, began: datetime, ended: datetime):
        report = VectorReport(**{**self.prof_fields,
                                 "top_str": self.top_dir,
                                 "ref_str": self.ref_seq,
                                 "began": began,
                                 "ended": ended})
        report.save()
        return report

    def _write_batch(self,
                     read_names: list[str, ...],
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
        df.to_orc(mv_file.path, engine="pyarrow")
        return mv_file, n_records

    def _vectorize_record(self, rec: SamRecord):
        """
        Compute the mutation vector of one record from a SAM file.

        ** Arguments **
        rec (SamRecord) --> SAM record for which to compute a mutation vector

        ** Returns **
        muts (bytearray) <- mutation vector
        """
        if rec.ref != self.ref:
            raise ValueError(f"Ref name of record ({rec.ref}) failed to "
                             f"match ref name of profile ({self.ref})")
        muts = rec.vectorize(self.seq_bytes, self.end5, self.end3)
        if not any(muts):
            raise ValueError("Mutation vector was entirely blank.")
        return rec.read_name, muts

    def _vectorize_batch(self, sam_viewer: SamViewer, batch_num: int,
                         start: int, stop: int):
        """
        Generate a batch of mutation vectors and write them to an ORC file.

        ** Arguments **
        sam_viewer (SamViewer) -> viewer to the SAM file for which to generate
                                  a batch of mutation vectors
        batch_num (int) --------> non-negative integer label for the batch
        start (int) ------------> seek to this position in the SAM file, then
                                  start generating vectors; position must be
                                  the beginning of a line, otherwise integrity
                                  checks will fail
        stop (int) -------------> stop generating vectors upon reaching this
                                  position in the SAM file; position must be
                                  the beginning of a line, otherwise integrity
                                  checks will fail
        
        ** Returns **
        n_records (int) <-------- number of records read from the SAM file
                                  between positions start and stop
        checksum (str) <--------- MD5 checksum of the ORC file of vectors
        """
        if stop > start:
            with sam_viewer as sv:
                # Use the SAM viewer to generate the mutation vectors.
                # Collect them as a single, 1-dimensional bytes object.
                read_names, muts = zip(*map(self._vectorize_record,
                                            sv.get_records(start, stop)))
        else:
            read_names, muts = (), ()
        # Write the mutation vectors to a file and compute its checksum.
        mv_file, n_records = self._write_batch(read_names, muts, batch_num)
        checksum = self.digest_file(mv_file.path)
        return n_records, checksum

    def _vectorize_sam(self, max_cpus: int):
        if self._vectorized:
            raise RuntimeError(f"Vectoring was already run for {self}.")
        self._vectorized = True
        with SamViewer(top_dir=self.top_dir,
                       max_cpus=max_cpus,
                       xam_path=self.bam_path,
                       ref_name=self.ref,
                       end5=self.end5,
                       end3=self.end3,
                       spanning=self.spanning,
                       min_qual=self.min_qual) as sv:
            batch_size = max(1, DEFAULT_BATCH_SIZE // self.length)
            indexes = list(sv.get_batch_indexes(batch_size))
            starts = indexes[:-1]
            stops = indexes[1:]
            self.num_batches = len(starts)
            svs = [SamViewer(top_dir=self.top_dir,
                             max_cpus=max_cpus,
                             xam_path=sv.sam_path,
                             ref_name=self.ref,
                             end5=self.end5,
                             end3=self.end3,
                             spanning=self.spanning,
                             min_qual=self.min_qual,
                             owner=False)
                   for _ in self.batch_nums]
            args = list(zip(svs, self.batch_nums, starts, stops, strict=True))
            if max_cpus > 1 and self.num_batches > 1:
                n_procs = min(max_cpus, self.num_batches)
                with Pool(n_procs, maxtasksperchild=1) as pool:
                    results = list(pool.starmap(self._vectorize_batch, args,
                                                chunksize=1))
            else:
                results = list(starmap(self._vectorize_batch, args))
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
    
    def vectorize(self, max_cpus: int):
        if self.rerun or not self.outputs_valid:
            self.batch_dir.path.mkdir(parents=True, exist_ok=True)
            logging.info(f"{self}: computing vectors")
            began = datetime.now()
            self._vectorize_sam(max_cpus)
            ended = datetime.now()
            logging.info(f"{self}: writing report")
            self._write_report(began, ended)
            logging.info(f"{self}: finished")
        else:
            logging.warning(f"{self}: already finished. To rerun, add --rerun")
        return self.report_path


class VectorWriterSpawner(object):
    def __init__(self,
                 top_dir: str,
                 fasta: str,
                 bam_files: list[str],
                 coords: list[tuple[str, int, int]],
                 primers: list[tuple[str, DNA, DNA]],
                 fill: bool,
                 parallel: str,
                 max_cpus: int,
                 min_phred: int,
                 phred_enc: int,
                 rerun: bool):
        self.top_dir = path.TopDirPath.parse_path(top_dir)
        self.bam_paths = [path.OneRefAlignmentInFilePath.parse_path(bam)
                          for bam in bam_files]
        self.ref_path = path.RefsetSeqInFilePath.parse_path(fasta)
        self.coords = coords
        self.primers = primers
        self.fill = fill
        self.parallel = parallel
        self.max_cpus = max_cpus
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
        seqs = dict(FastaParser(self.ref_path.path).parse())
        if not seqs:
            raise ValueError(f"'{self.ref_path}' contained no sequences")
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
                                    ref=ref, end5=first, end3=last))
        for ref, fwd, rev in self.primers:
            add_region(RegionFinder(ref_seq=self.ref_seqs.get(ref),
                                    ref=ref, fwd=fwd, rev=rev))
        if self.fill:
            for ref, seq in self.ref_seqs.items():
                if ref not in regions:
                    add_region(RegionFinder(ref_seq=seq, ref=ref))
        return regions

    def _get_writers(self):
        for bam in self.bam_paths:
            for region in self.regions[bam.ref]:
                if region.ref != bam.ref:
                    raise RuntimeError(f"Refs of region ({region.ref}) "
                                       f"and BAM file ({bam.ref}) disagree.")
                yield VectorWriter(top_dir=self.top_dir,
                                   bam_path=bam,
                                   ref_seq=self.ref_seqs[bam.ref],
                                   end5=region.end5,
                                   end3=region.end3,
                                   min_qual=self.min_qual,
                                   rerun=self.rerun)

    @cached_property
    def writers(self):
        return list(self._get_writers())

    @property
    def num_writers(self):
        return len(self.writers)

    @property
    def parallel_broad(self):
        return (self.parallel == ParallelChoice.BROAD
                or (self.parallel == ParallelChoice.AUTO
                    and self.num_writers > 1))

    @property
    def parallel_deep(self):
        return (self.parallel == ParallelChoice.DEEP
                or (self.parallel == ParallelChoice.AUTO
                    and self.num_writers == 1))

    @cached_property
    def procs_per_profile(self):
        return self.max_cpus if self.parallel_deep else 1

    def generate_profiles(self):
        if self.num_writers == 0:
            raise ValueError("No BAM files and/or regions specified")
        args = tuple(zip(self.writers, repeat(self.procs_per_profile)))
        if self.parallel_broad:
            n_procs = max(min(self.max_cpus, self.num_writers), 1)
            with Pool(n_procs, maxtasksperchild=1) as pool:
                report_files = tuple(pool.starmap(VectorWriter.vectorize, args,
                                                  chunksize=1))
        else:
            report_files = tuple(starmap(VectorWriter.vectorize, args))
        return report_files


class VectorReader(VectorIO):
    @property
    def shape(self):
        return self.num_vectors, self.length

    @property
    def vectors(self):
        # FIXME: implement method to load the mutation vectors efficiently.
        return

    @classmethod
    def from_report(cls, report: VectorReport):
        return cls(top_dir=report.top_dir,
                   sample=report.sample,
                   ref_seq=report.ref_seq,
                   ref=report.ref,
                   end5=report.end5,
                   end3=report.end3)

    @classmethod
    def from_report_file(cls, report_file):
        return cls.from_report(VectorReport.load(report_file))
