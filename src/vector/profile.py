from __future__ import annotations

from collections import defaultdict, namedtuple
from datetime import datetime
from functools import cache, cached_property, partial
from hashlib import md5
import itertools
import logging
from multiprocessing import Pool
import re
import sys
import time
from typing import Any, ClassVar, Sequence

import numpy as np
import pandas as pd
from pydantic import (BaseModel, Extra, Field, NonNegativeInt, NonNegativeFloat,
                      PositiveInt, StrictBool, StrictStr)
from pydantic import validator, root_validator

from ..util import path
from ..util.cli import (opt_bamf, opt_bamd,
                        opt_library, opt_cfill,
                        opt_coords, opt_primers, opt_primer_gap,
                        opt_out_dir, opt_temp_dir,
                        opt_phred_enc, opt_min_phred,
                        opt_strict_pairs, opt_ambid, opt_batch_size,
                        opt_parallel, opt_max_procs,
                        opt_rerun, opt_resume, opt_save_temp)
from ..util.seq import (BLANK_INT, MATCH_INT, DELET_INT, INS_5_INT, INS_3_INT,
                        SUB_A_INT, SUB_C_INT, SUB_G_INT, SUB_T_INT, AMBIG_INT,
                        DNA, FastaParser)
from ..util.util import get_num_parallel
from ..vector.samread import SamReader
from ..vector.vector import SamRecord

RegionTuple = namedtuple("PrimerTuple", ["pos5", "pos3"])


def get_bytes_per_patch(batch_size: float):
    """
    Return the number of bytes per batch of a given size in mebibytes.

    Parameters
    ----------
    batch_size: float
        Size of the batch in mebibytes (MiB: 1 MiB = 2^20 bytes)

    Return
    ------
    int
        Number of bytes per batch, to the nearest integer
    """
    return round(batch_size * 1048576)  # 1048576 = 2^20


class Region(object):
    """
    Represent a region of a reference sequence between two coordinates.

    Attributes
    ----------
    ref: str
        Name of the reference sequence
    end5: int (1 ≤ end5 ≤ end3)
        Coordinate of the reference sequence at which the region's
        5' end is located (1-indexed)
    end3: int (end5 ≤ end3 ≤ len(ref_seq))
        Coordinate of the reference sequence at which the region's
        3' end is located (1-indexed; end3 itself is included)
    seq: DNA
        Sequence of the region between end5 and end3 (inclusive)
    spans: bool
        Whether the region spans the entire reference sequence

    Examples
    --------
    >>> seq = DNA(b"CATCTGGA")
    >>> name = "example"
    >>> region = Region(ref_seq=seq, ref=name, end5=1, end3=8)
    >>> assert region.seq == seq
    >>> region = Region(ref_seq=seq, ref=name, end5=1, end3=-1)
    >>> assert region.seq == seq
    >>> region = Region(ref_seq=seq, ref=name, end5=-8, end3=8)
    >>> assert region.seq == seq
    >>> region = Region(ref_seq=seq, ref=name, end5=3, end3=7)
    >>> assert region.seq == DNA(b"TCTGG")
    >>> region = Region(ref_seq=seq, ref=name, end5=-5, end3=-3)
    >>> assert region.seq == DNA(b"CTG")
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
            Entire reference sequence (not of just the region)
        ref: str
            Name of the reference sequence
        end5: int (-len(ref_seq) ≤ end5 ≤ len(ref_seq); end5 ≠ 0)
            Coordinate of the reference sequence at which the 5' end of
            the region is located. If positive, number the coordinates
            of the reference sequence 1, 2, ... starting at the 5' end
            (i.e. 1-based indexing). If negative, number the coordinates
            of the reference sequence -1, -2, ... starting at the 3' end
            (i.e. 1-based indexing from the other side), then convert to
            the corresponding (positive) 1-based index from the 5' end
        end3: int (-len(ref_seq) ≤ end3 ≤ len(ref_seq); end3 ≠ 0)
            Coordinate of the reference sequence at which the region's
            3' end is located. Follows the same coordinate numbering
            convention as end5
        """
        if end5 < 0:
            # If end5 < 0, find the corresponding positive coordinate.
            end5 += len(ref_seq) + 1
        if end3 < 0:
            # If end3 < 0, find the corresponding positive coordinate.
            end3 += len(ref_seq) + 1
        if not 1 <= end5 <= end3 <= len(ref_seq):
            raise ValueError("Must have 1 ≤ end5 ≤ end3 ≤ len(ref_seq), "
                             f"but got end5 = {end5}, end3 = {end3}, and "
                             f"len(ref_seq) = {len(ref_seq)}")
        self.ref = ref
        self.end5 = end5
        self.end3 = end3
        self.seq = ref_seq[end5 - 1: end3]
        self.spans = self.seq == ref_seq

    @property
    def path_fields(self):
        return {"ref": self.ref, "end5": self.end5, "end3": self.end3}

    @property
    def report_fields(self):
        return {"ref": self.ref, "end5": self.end5, "end3": self.end3,
                "spans": self.spans}

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
        return {"module": path.Module.VECTOR,
                "sample": self.sample,
                **super().path_fields}

    @property
    def report_fields(self):
        return {"sample": self.sample,
                **super().report_fields}

    @property
    def tag(self):
        return tuple([self.sample, *super().tag])

    def __str__(self):
        return f"{self.sample}:{super().__str__()}"


class VectorIO(MutationalProfile):
    """
    Read and write mutation mut_vectors.

    Fields
    ------
    top_dir: str
        Path to the top-level directory in which all temporary and final
        output files will be written
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.num_batches: int = 0
        self.num_vectors: int = 0
        self.checksums: list[str] = list()

    @property
    def report_fields(self):
        return {"num_batches": self.num_batches,
                "num_vectors": self.num_vectors,
                "checksums": self.checksums,
                **super().report_fields}

    def get_batch_dir(self, out_dir: str):
        return path.RegionOutDirPath(top=out_dir, **self.path_fields)

    def get_mv_batch_path(self, out_dir: str, batch: int):
        return path.MutVectorBatchFilePath(top=out_dir,
                                           **self.path_fields,
                                           batch=batch,
                                           ext=".orc")

    @property
    def batch_nums(self) -> list[int]:
        """ Return a range of all batch numbers. """
        return list(range(self.num_batches))

    def get_mv_batch_paths(self, out_dir: str):
        """ Return the path of every mutation vector batch file. """
        return [self.get_mv_batch_path(out_dir, batch)
                for batch in self.batch_nums]

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
    Read and write a report about a mutational profile, including:
    - the sample, reference, and region
    - number of mutation mut_vectors
    - number of mutation vector batch files and their checksums
    - beginning and ending time, duration, and speed of vectoring

    Examples
    --------
    >>> report = VectorReport(sample="dms2", ref="tmv-rna", end5=1, end3=20,
    ...                       seq_str=DNA(b"GTATTTTTACAACAATTACC"),
    ...                       num_vectors=10346, num_batches=2,
    ...                       checksums=["b47260fcad8be60470bee67862f187b4",
    ...                                  "098f40cfadc266ea5bc48ab2e18cdc95"],
    ...                       began=datetime.now(),
    ...                       ended=(time.sleep(1E-5), datetime.now())[-1])
    >>> report.seq_str
    'GTATTTTTACAACAATTACC'
    """

    class Config:
        allow_population_by_field_name = True
        extra = Extra.ignore

    # Fields
    sample: StrictStr = Field(alias="Sample name")
    ref: StrictStr = Field(alias="Reference name")
    end5: PositiveInt = Field(alias="5' end of region")
    end3: PositiveInt = Field(alias="3' end of region")
    seq_str: StrictStr = Field(alias="Sequence of region")
    spans: StrictBool = Field(alias="Region spans entire reference sequence")
    num_vectors: NonNegativeInt = Field(alias="Number of mut_vectors")
    num_batches: NonNegativeInt = Field(alias="Number of batches")
    checksums: list[StrictStr] = Field(alias="MD5 checksums of vector batches")
    began: datetime = Field(alias="Began vectoring")
    ended: datetime = Field(alias="Ended vectoring")
    duration: NonNegativeFloat = Field(default=float("nan"),
                                       alias="Duration of vectoring (s)")
    speed: NonNegativeFloat = Field(default=float("nan"),
                                    alias="Speed of vectoring (mut_vectors/s)")

    # Format of dates and times in the report file
    dt_fmt: ClassVar[str] = "on %Y-%m-%d at %H:%M:%S.%f"

    @validator("seq_str", pre=True)
    def convert_ref_seq_to_str(cls, seq_str: DNA):
        """ Return reference sequence (DNA) as a string.
        Must be str in order to write to and load from JSON correctly. """
        return str(seq_str)

    @property
    def ref_seq(self):
        """ Return reference sequence string (from JSON) as a DNA object.
        The methods from VectorIO expect ref_seq to be of type DNA. """
        return DNA(self.seq_str.encode())

    @root_validator(pre=False)
    def calculate_duration_and_speed(cls, values):
        """
        Calculate and return the duration and speed of vectoring:
        - duration: difference between ending and beginning time (sec)
        - speed: number of mut_vectors processed per unit time (mut_vectors/sec)
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
        # - nan (i.e. 0 / 0) if no mut_vectors were processed
        # when calculating the speed of processing mut_vectors.
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
        end5, end3 = values["end5"], values["end3"]
        length = len(values["seq_str"])
        if end5 < 1 or length != end3 - end5 + 1:
            raise ValueError(f"Invalid region '{end5}-{end3}' for reference "
                             f"sequence of length {length}")
        return values

    def find_invalid_batches(self, out_dir: str, validate_checksums: bool):
        """ Return all the batches of mutation mut_vectors that either do not exist
        or do not match their expected checksums. """
        missing = list()
        badsum = list()
        for file, checksum in zip(self.get_mv_batch_paths(out_dir),
                                  self.checksums,
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

    def save(self, out_dir: str):
        report_path = path.MutVectorReportFilePath(top=out_dir,
                                                   module=path.Module.VECTOR.value,
                                                   sample=self.sample,
                                                   ref=self.ref,
                                                   end5=self.end5,
                                                   end3=self.end3,
                                                   ext=".json")
        report_path.path.parent.mkdir(parents=True, exist_ok=True)
        with open(report_path.path, "w") as f:
            f.write(self.json(by_alias=True))
        return report_path

    @classmethod
    def load(cls, report_file: str, validate_checksums: bool = True):
        """ Load a mutation vector report from a file. """
        try:
            report = cls.parse_file(report_file)
        except (FileNotFoundError, path.PathError):
            return
        report_path = path.MutVectorReportFilePath.parse(report_file)
        if (report_path.sample != report.sample
                or report_path.ref != report.ref
                or report_path.end5 != report.end5
                or report_path.end3 != report.end3):
            logging.error(f"Report fields do not match path '{report_path}'")
            return
        missing, badsum = report.find_invalid_batches(report_path.top,
                                                      validate_checksums)
        if missing or badsum:
            return
        return report


class VectorWriter(VectorIO):
    """
    Compute mutation mut_vectors for all reads from one sample mapping to one
    region of one reference sequence.
    """

    def __init__(self, *, bam_path: path.OneRefAlignmentInFilePath, **kwargs):
        super().__init__(sample=bam_path.sample,
                         ref=bam_path.ref,
                         **kwargs)
        self.bam_path = bam_path
        self.seq_bytes = bytes(self.seq)

    def get_report_path(self, out_dir: str):
        return path.MutVectorReportFilePath(top=out_dir,
                                            **self.path_fields,
                                            ext=".json")

    def _write_report(self, out_dir: str,
                      began: datetime, ended: datetime):
        report = VectorReport(**{**self.report_fields,
                                 "seq_str": self.seq,
                                 "began": began,
                                 "ended": ended})
        report_path = report.save(out_dir)
        return report_path

    def _write_batch(self,
                     out_dir: str,
                     batch_num: int,
                     read_names: list[str],
                     mut_vectors: tuple[bytearray, ...]):
        """
        Write a batch of mutation vectors to an ORC file.

        Parameters
        ----------
        out_dir: str
            Top-level output directory
        batch_num: int (≥ 0)
            Non-negative integer label for the batch
        read_names: list[str]
            Names of the reads (also of the vectors) in the batch
        mut_vectors: NDArray
            Batch of mutation vectors in which each row is a mutation
            vector and each column is a position in the region

        Return
        ------
        tuple[path.MutVectorBatchFilePath, int]
            [0] Path to which the batch of vectors was written
            [1] Number of vectors in the batch
        """
        # Process the mutation vectors into a 2D NumPy array.
        mut_array = np.frombuffer(b"".join(mut_vectors), dtype=np.byte)
        n_records = len(mut_array) // self.length
        mut_array.shape = (n_records, self.length)
        # Data must be converted to pd.DataFrame for PyArrow to write.
        # Explicitly set copy=False to copying the mutation mut_vectors.
        mut_frame = pd.DataFrame(data=mut_array,
                                 index=read_names,
                                 columns=self.columns,
                                 copy=False)
        mv_file = self.get_mv_batch_path(out_dir, batch_num)
        mv_file.path.parent.mkdir(parents=True, exist_ok=True)
        mut_frame.to_orc(mv_file.path, index=True, engine="pyarrow")
        return mv_file, n_records

    def _vectorize_record(self, rec: SamRecord, **kwargs):
        """
        Compute the mutation vector of one record from a SAM file.

        ** Arguments **
        rec (SamRecord) --> SAM record for which to compute a mutation vector

        ** Returns **
        mut_vectors (bytearray) <- mutation vector
        """
        try:
            ref_name = rec.read1.rname.decode()
            read_name = rec.read1.qname.decode()
            if ref_name != self.ref:
                raise ValueError(
                    f"Read '{read_name}' had reference '{ref_name}' different "
                    f"from profile reference '{self.ref}'")
            muts = rec.vectorize(region_seq=self.seq_bytes,
                                 region_end5=self.end5,
                                 region_end3=self.end3,
                                 **kwargs)
            if not any(muts):
                raise ValueError(f"Vector for read '{read_name}' was blank")
        except ValueError as error:
            logging.error(f"Read '{rec.read1.qname.decode()}' failed to "
                          f"vectorize due to the following error: {error}")
            return "", bytearray()
        return read_name, muts

    def _vectorize_records(self, /,
                           reader: SamReader,
                           start: int,
                           stop: int, *,
                           strict_pairs: bool = opt_strict_pairs.default,
                           phred_enc: int = opt_phred_enc.default,
                           min_phred: int = opt_min_phred.default,
                           ambid: bool = opt_ambid.default):
        """
        Generate a batch of mutation mut_vectors and write them to a file.

        Parameters
        ----------
        reader: SamReader
            Reader for the SAM file for which to generate a batch of
            mutation mut_vectors
        start: int (≥ 0)
            Seek to this position in the SAM file, then start generating
            mut_vectors; position must be the beginning of a line, else an
            error will be raised.
        stop: int (≥ start)
            Stop generating mut_vectors upon reaching this position in the
            SAM file; position must be the beginning of a line, else an
            error will be raised.
        

        """
        if stop > start:
            with reader as reading:
                # Use the SAM reader to generate the mutation mut_vectors.
                # Collect them as a single, 1-dimensional bytes object.
                vectorize_record = partial(self._vectorize_record,
                                           min_qual=get_min_qual(min_phred,
                                                                 phred_enc),
                                           ambid=ambid)
                read_names, muts = zip(*map(vectorize_record,
                                            reading.get_records(start, stop,
                                                                strict_pairs)))
                # For every read for which creating a mutation vector
                # failed, an empty string was returned as the read name
                # and an empty bytearray as the mutation vector. The
                # empty read names must be filtered out, while the empty
                # mutation mut_vectors will not cause problems because,
                # being of length zero, they will effectively disappear
                # when all the mut_vectors are concatenated into a 1D array.
                read_names = tuple(filter(None, read_names))
        else:
            logging.warning(f"Stop position ({stop}) is not after "
                            f"start position ({start})")
            read_names, muts = (), ()
        return read_names, muts

    def _vectorize_batch(self, /,
                         out_dir: str,
                         batch_num: int,
                         reader: SamReader,
                         start: int,
                         stop: int,
                         **kwargs):
        """

        Return
        ------
        tuple[int, str]
            - Number of records read from the SAM file between positions
              start and stop
            - MD5 checksum of the ORC file of the batch of mut_vectors
        """
        # Compute the read names and mutation mut_vectors.
        read_names, muts = self._vectorize_records(reader, start, stop,
                                                   **kwargs)
        # Write them to a file and compute its checksum.
        mv_file, n_records = self._write_batch(out_dir, batch_num,
                                               list(read_names), muts)
        checksum = self.digest_file(mv_file.path)
        return n_records, checksum

    def _vectorize_sam(self, /, *,
                       out_dir: str = opt_out_dir.default,
                       temp_dir: str = opt_temp_dir.default,
                       batch_size: int = get_bytes_per_patch(
                           opt_batch_size.default),
                       n_procs: int = opt_max_procs.default,
                       save_temp: bool = opt_save_temp.default,
                       rerun: bool = opt_rerun.default,
                       resume: bool = opt_resume.default,
                       **kwargs):
        if self.outputs_valid(out_dir) and not rerun:
            logging.warning(f"Skipping vectorization of {self} because output "
                            f"files already exist. To rerun, add --rerun")
            return
        # Open the primary SAM file reader to write the subset of SAM
        # records to a temporary SAM file and determine the number and
        # start/stop indexes of each batch of records in the file.
        with SamReader(xam_path=self.bam_path,
                       end5=self.end5,
                       end3=self.end3,
                       spans=self.spans,
                       temp_dir=temp_dir,
                       n_procs=n_procs,
                       save_temp=save_temp,
                       resume=resume,
                       owner=True) as reader:
            # Use integer division to round down the number of mut_vectors
            # per batch to avoid exceeding the given batch size. But
            # also ensure that there is at least one vector per batch.
            vectors_per_batch = max(1, batch_size // self.length)
            # Compute for each batch the positions in the SAM file
            # (given by file.tell and used by file.seek) at which the
            # batch starts and stops.
            indexes = list(reader.get_batch_indexes(vectors_per_batch))
            starts = indexes[:-1]
            stops = indexes[1:]
            self.num_batches = len(starts)
            # Once the number of batches has been determined, a list of
            # new SAM file readers is created. Each is responsible for
            # converting one batch of SAM reads into one batch of
            # mutation mut_vectors. Setting owner=False prevents the SAM
            # reader from creating a new SAM file (as well as from
            # deleting it upon exit); instead, it uses the file written
            # by the primary SAM reader (which is reader.sam_path).
            readers = [SamReader(xam_path=reader.sam_path,
                                 end5=self.end5,
                                 end3=self.end3,
                                 spans=self.spans,
                                 temp_dir=temp_dir,
                                 n_procs=n_procs,
                                 save_temp=save_temp,
                                 resume=resume,
                                 owner=False)
                       for _ in self.batch_nums]
            # Create lists of the positional and keyword aguments for
            # each SamReader.
            iter_args = [(out_dir, batch, reader, start, stop)
                         for batch, reader, start, stop in
                         zip(self.batch_nums, readers, starts, stops,
                             strict=True)]
            # Then pass the arguments to _vectorize_batch.
            vectorize_batch = partial(self._vectorize_batch, **kwargs)
            if (pool_size := min(n_procs, self.num_batches)) > 1:
                with Pool(pool_size) as pool:
                    results = list(pool.starmap(vectorize_batch, iter_args))
            else:
                results = list(itertools.starmap(vectorize_batch, iter_args))
            # results is a list containing, for each batch, a tuple of
            # the number of mut_vectors in the batch and the checksum of
            # the batch. Gather them into num_vectors and checksums.
            for num_vectors, checksum in results:
                self.num_vectors += num_vectors
                self.checksums.append(checksum)
            return self.get_report_path(out_dir)

    def outputs_valid(self, out_dir: str):
        return VectorReport.load(self.get_report_path(out_dir).path) is not None

    def vectorize(self, **kwargs):
        began = datetime.now()
        report_path = self._vectorize_sam(**kwargs)
        ended = datetime.now()
        if report_path is not None:
            written = self._write_report(report_path.top, began, ended)
            if written != report_path:
                logging.critical("Intended and actual paths of report differ: "
                                 f"{report_path} ≠ {written}")
        return report_path


class VectorReader(VectorIO):
    INDEX_COL = "__index_level_0__"

    def __init__(self,
                 out_dir: str,
                 num_vectors: int,
                 num_batches: int,
                 **kwargs):
        super().__init__(**kwargs)
        self.out_dir = out_dir
        self.num_vectors = num_vectors
        self.num_batches = num_batches

    @property
    def shape(self):
        return self.num_vectors, self.length

    def get_batch(self,
                  batch: int,
                  positions: Sequence[int] | None = None):
        """
        Return the mutation mut_vectors from one batch. Optionally, select
        a subset of the columns of the mutation mut_vectors.

        Parameters
        ----------
        batch: int (≥ 0)
            Number of the batch of mutation mut_vectors
        positions: sequence[int] | None (default: None)
            If given, use only these positions from the mutation mut_vectors;
            otherwise, use all positions.

        Return
        ------
        DataFrame
            Mutation mut_vectors; each row is a vector indexed by its name,
            each column a position indexed by its positional number
        """
        if batch not in self.batch_nums:
            raise ValueError(f"Invalid batch number: {batch} "
                             f"(expected one of {list(self.batch_nums)})")
        columns = ([self.INDEX_COL] + self.pos_to_cols(positions) if positions
                   else None)
        vectors = pd.read_orc(self.get_mv_batch_path(self.out_dir, batch).path,
                              columns=columns)
        vectors.set_index(self.INDEX_COL, drop=True, inplace=True)
        vectors.columns = positions if positions else self.positions
        return vectors

    def get_all_batches(self, positions: Sequence[int] | None = None):
        """
        Yield every batch of mutation mut_vectors.

        Parameters
        ----------
        positions: sequence[int] | None (default: None)
            If given, use only these position from the mutation mut_vectors;
            otherwise, use all positions.

        Yield
        -----
        DataFrame
            Mutation mut_vectors; each row is a vector indexed by its name,
            each column a position indexed by its base and number
        """
        for batch in self.batch_nums:
            yield self.get_batch(batch, positions)

    def get_all_vectors(self, positions: Sequence[int] | None = None):
        """
        Return all mutation mut_vectors for this vector reader. Note that
        reading all mut_vectors could take more than the available memory
        and cause the program to crash. Thus, use this method only if
        all mut_vectors will fit into memory. Otherwise, use the method
        ```get_all_batches``` to process the mut_vectors in small batches.

        Parameters
        ----------
        positions: sequence[int] | None (default: None)
            If given, use only these position from the mutation mut_vectors;
            otherwise, use all positions.

        Return
        ------
        DataFrame
            Mutation mut_vectors; each row is a vector indexed by its name,
            each column a position indexed by its base and number
        """
        return pd.concat(self.get_all_batches(positions), axis=0)

    @staticmethod
    def _query_vectors(vectors: pd.DataFrame, query: int) -> pd.DataFrame:
        """
        Return a boolean array of the same shape as mut_vectors where
        element i,j is True if and only if the byte at element i,j of
        mut_vectors is both non-zero and a bitwise subset of query (a byte
        equal to query also counts as a subset).

        Parameters
        ----------
        vectors: DataFrame
            Mutation mut_vectors
        query: int
            Byte to query

        Return
        ------
        DataFrame
            Boolean type DataFrame of the same shape as mut_vectors where
            each element is True if the element at the same position in
            mut_vectors matched the query and False otherwise
        """
        # Flag as False all bytes in mut_vectors that have no bits set to 1,
        # since these bytes represent positions in mut_vectors that were not
        # covered by reads and thus should not count as query matches.
        covered = vectors.astype(bool, copy=False)
        if query == AMBIG_INT:
            # If the query byte is all 1s (i.e. decimal 255), then the
            # next step (bitwise OR) will return True for every byte,
            # so the return value will equal that of covered. For the
            # sake of speed, return covered now.
            # Note that for the sake of speed and memory, covered is
            # a view to the SAME DataFrame as mut_vectors. Thus, functions
            # that use covered should merely read it, NEVER modify it.
            return covered
        # Flag as True all bytes in mut_vectors that are subsets of the
        # query byte (including, for now, bytes in vector that are
        # 00000000). In order for a byte in mut_vectors to be a subset of
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
        mut_vectors that match the query.

        Parameters
        ----------
        query: int (0 ≤ query < 256)
            Query byte: to match, a byte in the vector must be equal to
            or a bitwise subset of the query byte, and must be non-zero.
        positions: sequence[int] | None (default: None)
            If given, use only these positions from the mutation mut_vectors;
            otherwise, use all positions.

        Return
        ------
        Series
            Number of mut_vectors matching the query at each position; index
            is the position in the vector.
        """
        # For each batch of mut_vectors, get a DataFrame of boolean values
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
    def load(cls, report_file: str, validate_checksums: bool = True):
        if not (report := VectorReport.load(report_file, validate_checksums)):
            logging.critical(f"Failed to load report from {report_file}")
            return
        return cls(out_dir=path.MutVectorReportFilePath.parse(report_file).top,
                   sample=report.sample,
                   ref_seq=report.ref_seq,
                   ref=report.ref,
                   end5=report.end5,
                   end3=report.end3,
                   num_vectors=report.num_vectors,
                   num_batches=report.num_batches)


def get_min_qual(min_phred: int, phred_enc: int):
    """
    Return the minimum quality for a base in a read to be considered
    informative, as the ASCII integer corresponding to the character
    in the FASTQ file that is the minimum valid quality.

    Return
    ------
    int
        The ASCII value corresponding to the character in the FASTQ
        file read quality string that represents the minimum quality
        for a base to be considered informative

    Examples
    --------
    For example, if the minimum Phred score (```min_phred```) that
    is accepted as informative is 20, and the Phred encoding of the
    FASTQ file (```phred_enc```) is 33 (i.e. ASCII+33), then the
    minimum quality as an ASCII integer (```min_qual```) is 20 + 33
    = 53, which is character '5'. If ```min_phred``` were 37, then
    ```min_qual``` would be 37 + 33 = 70, which is character 'F'.
    """
    return min_phred + phred_enc


def get_regions(ref_seqs: dict[str, DNA], /,
                coords: tuple[tuple[str, int, int]],
                primers: tuple[tuple[str, DNA, DNA]], *,
                primer_gap: int = opt_primer_gap.default,
                cfill: bool = opt_cfill.default):
    regions: dict[str, list[RegionFinder]] = defaultdict(list)

    def add_region(region: RegionFinder):
        if any(region == other for other in regions[region.ref]):
            logging.warning(f"Skipping duplicate region: {region.ref_coords}")
        regions[region.ref].append(region)

    for ref, first, last in coords:
        add_region(RegionFinder(ref_seq=ref_seqs.get(ref),
                                ref=ref, end5=first, end3=last,
                                primer_gap=primer_gap))
    for ref, fwd, rev in primers:
        add_region(RegionFinder(ref_seq=ref_seqs.get(ref),
                                ref=ref, fwd=fwd, rev=rev,
                                primer_gap=primer_gap))
    if cfill:
        for ref, seq in ref_seqs.items():
            if ref not in regions:
                add_region(RegionFinder(ref_seq=seq, ref=ref,
                                        primer_gap=primer_gap))
    return regions


def get_writers(refset_path: path.RefsetSeqInFilePath,
                bam_paths: list[path.OneRefAlignmentInFilePath], /,
                **kwargs):
    ref_seqs = dict(FastaParser(refset_path.path).parse())
    regions = get_regions(ref_seqs, **kwargs)
    writers: dict[tuple, VectorWriter] = dict()
    for bam in bam_paths:
        for region in regions[bam.ref]:
            if region.ref != bam.ref:
                logging.error(f"Skipping region {region} of {bam.path} "
                              "because its reference does not match that "
                              f"of the BAM file ('{bam.ref}').")
                continue
            writer = VectorWriter(bam_path=bam,
                                  ref_seq=ref_seqs[bam.ref],
                                  end5=region.end5,
                                  end3=region.end3)
            if writer.tag in writers:
                logging.warning("Skipping duplicate mutational profile: "
                                f"{writer}.")
                continue
            writers[writer.tag] = writer
    return list(writers.values())


def generate_profiles(writers: list[VectorWriter], /, *,
                      phred_enc: int = opt_phred_enc.default,
                      min_phred: int = opt_min_phred.default,
                      max_procs: int = opt_max_procs.default,
                      parallel: bool = opt_parallel.default,
                      **kwargs):
    n_profiles = len(writers)
    if n_profiles == 0:
        logging.critical("No BAM files and/or regions specified")
        return ()
    # Determine method of parallelization. Do not use hybrid mode, which
    # would try to process multiple SAM files in parallel and use more
    # than one processe for each file. Python ```multiprocessing.Pool```
    # forbids a daemon process (one for each SAM file in parallel) from
    # spawning additional processes.
    n_tasks_parallel, n_procs_per_task = get_num_parallel(n_profiles,
                                                          max_procs,
                                                          parallel)
    # List the arguments of VectorWriter.vectorize for each writer.
    iter_args = [(writer,) for writer in writers]
    vectorize = partial(VectorWriter.vectorize,
                        phred_enc=phred_enc,
                        min_phred=min_phred,
                        n_procs=n_procs_per_task,
                        **kwargs)
    # Call the vectorize method of each writer, passing args.
    if n_tasks_parallel > 1:
        with Pool(n_tasks_parallel) as pool:
            report_files = list(pool.starmap(vectorize, iter_args))
    else:
        report_files = list(itertools.starmap(vectorize, iter_args))
    return report_files


vector_trans_table = bytes.maketrans(*map(b"".join, zip(*[(
    i.to_bytes(length=1, byteorder=sys.byteorder),
    (
        b"." if i == BLANK_INT
        else b"~" if i == MATCH_INT
        else b"/" if i == DELET_INT
        else b"{" if i == (INS_5_INT | MATCH_INT)
        else b"}" if i == (INS_3_INT | MATCH_INT)
        else b"A" if i == SUB_A_INT
        else b"C" if i == SUB_C_INT
        else b"G" if i == SUB_G_INT
        else b"T" if i == SUB_T_INT
        else b"?"
    )
) for i in range(256)])))


def trans_vectors_iter(vectors: pd.DataFrame):
    for index, row in zip(vectors.index, vectors.values, strict=True):
        translated = row.tobytes(order='C').translate(vector_trans_table)
        yield b"%b\t%b\n" % (index.encode(), translated)


def trans_vectors_block(vectors: pd.DataFrame):
    return b"".join(trans_vectors_iter(vectors))
