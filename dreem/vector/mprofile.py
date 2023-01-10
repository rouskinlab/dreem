from __future__ import annotations
from collections import defaultdict
from functools import cached_property
import itertools
import os
import pathlib
import re
from tqdm import tqdm
from datetime import datetime
from hashlib import file_digest
from multiprocessing import Pool
from typing import Any, List, Optional, Tuple, Dict

import numpy as np
import pandas as pd

from dreem.util.dflt import NUM_PROCESSES
from dreem.util.fa import FastaParser
from dreem.util.seq import DNA
from dreem.vector.samview import SamViewer
from dreem.vector.vector import SamRecord


DEFAULT_BATCH_SIZE = 2**25  # 33,554,432 bytes


def ref_from_bam_file(bam_file: str) -> bytes:
    """
    Get the name of the reference from the path of a BAM file. DREEM requires
    that the BAM file be named after the reference, up to the '.bam' extension,
    and that reference names be bytes objects.

    ** Arguments **
    bam_file (str) ---> path (absolute or relative) to the BAM file

    ** Returns **
    ref_name (bytes) <- name of the reference deduced from the BAM file
    """
    ref_name = pathlib.PosixPath(bam_file).stem.encode()
    return ref_name


def sample_from_bam_file(bam_file: str) -> str:
    """
    Get the name of the sample from the path of a BAM file. DREEM requires
    that the BAM file be in a directory named after the sample.

    ** Arguments **
    bam_file (str) ----> path (absolute or relative) to the BAM file

    ** Returns **
    sample_name (str) <- name of the sample deduced from the BAM file directory
    """
    sample_name = os.path.basename(os.path.dirname(os.path.abspath(bam_file)))
    return sample_name


def samples_from_bam_files(bam_files: List[str]) -> Dict[str, List[str]]:
    """
    Given a list of BAM files, return a dictionary where each key/value pair is
    a (key) sample and (value) list of BAM files deriving from that sample.

    ** Arguments **
    bam_files (list[str]) ----------> list of path(s) (absolute or relative)
                                      to n ≥ 0 BAM file(s)

    ** Returns **
    samples (dict[str, list[str]]) <- dict where each key is a sample name and
                                      maps to a list of absolute paths of all
                                      BAM files from that sample
    """
    samples = defaultdict(list)
    for bam in bam_files:
        samples[sample_from_bam_file(bam)].append(os.path.abspath(bam))
    return dict(samples)


class Region(object):
    __slots__ = ["_first", "_last", "_ref_seq", "_ref_name"]

    def __init__(self, ref_name: bytes, first: int, last: int, ref_seq: DNA):
        """
        Initialize a Region object to represent a region of interest (between a
        first and a last position) in a reference sequence.
        
        ** Arguments **
        ref_name (bytes) ----> name of the reference
        first (int) ---------> 5'-most position of the reference that lies in
                               the region of interest (1-indexed, inclusive)
        last (int) ----------> 3'-most position of the reference that lies in
                               the region of interest (1-indexed, inclusive)
        ref_seq (DNA) -------> sequence of the entire reference (not just the
                               region of interest)
        
        ** Returns **
        None
        """
        if first <= 0:
            raise ValueError(f"first ({first}) must be >= 1")
        if last > len(ref_seq):
            raise ValueError(f"last ({last}) must be <= the length of "
                             f"ref_seq ({len(ref_seq)}) {ref_seq}.")
        if -len(ref_seq) <= last < 0:
            # This option allows using non-positive end coordinates to mean
            # distance from the 3' end (similar to Python's indexing), with
            # -1 meaning up to and including the last coordinate, -2 meaning
            # up to an including the second-to-last coordinate, and so on
            # until -len(ref_seq), which means the first coordinate.
            last += len(ref_seq)
        if last < first:
            raise ValueError(f"last ({last}) must be >= first ({first})")
        self._first = first
        self._last = last
        self._ref_seq = ref_seq
        self._ref_name = ref_name
    
    @property
    def first(self):
        return self._first
    
    @property
    def last(self):
        return self._last
    
    @property
    def ref_seq(self):
        return self._ref_seq
    
    @property
    def ref_name(self):
        return self._ref_name
    
    @property
    def spanning(self) -> bool:
        """ Return whether the region spans the entire reference sequence. """
        return self.first == 1 and self.last == len(self.ref_seq)
    
    @property
    def region_seq(self) -> DNA:
        """ Return the sequence of the region of interest. """
        return self.ref_seq[self.first - 1: self.last]

    @property
    def length(self) -> int:
        """ Return the length of the region of interest. """
        return self.last - self.first + 1

    @property
    def positions(self) -> np.ndarray:
        """ Return all positions in the region of interest as an NDArray. """
        return np.arange(self.first, self.last + 1)

    @property
    def ref_coords(self) -> Tuple[bytes, int, int]:
        """ Return the name of the reference and the first and last positions
        of the region of interest; for equality testing and hashing. """
        return self.ref_name, self.first, self.last
    
    @property
    def columns(self) -> List[str]:
        """ Return a list of the bases and positions in the region of interest,
        each of the form '{base}{position}' (e.g. ['G13', 'C14', 'A15']). """
        return [f"{chr(base)}{pos}" for base, pos
                in zip(self.region_seq, self.positions)]
    
    def __eq__(self, other: object) -> bool:
        """ Two Region objects are equal iff their ref_coords match. """
        if isinstance(other, Region):
            return self.ref_coords == other.ref_coords
        else:
            return NotImplemented        


class RegionFinder(Region):
    # For amplicon-based DMS-MaPseq, primer_gap is the number of positions
    # to leave between the end of the fwd/rev primer in the RT-PCR step
    # and the 5'/3' end of the region. It is used because bases adjacent to
    # primers may show artifacts. Setting it to 0 causes the 5'/3' ends of
    # the region to lie immediately 3'/5' of the fwd/rev primers, respectively.
    # Historically, primer_gap has been set to 2.
    primer_gap = 2

    def __init__(self, ref_name: str, ref_seq: DNA,
                 first: Optional[int] = None, last: Optional[int] = None,
                 fwd: Optional[DNA] = None, rev: Optional[DNA] = None):
        """
        Initialize a RegionFinder object to determine the first and last
        positions in a region. The positions are determined as follows:
        - first
          - If first is given as an argument, this value is used for first.
          - Else if fwd is given, first is based on the location of its 3' end.
          - Else first is set to 1.
        - last
          - If last is given as an argument, last is set to this value.
          - Else if rev is given, last is based on the location of the 5' end
            of its reverse complement.
          - Else last is set to the length of ref_seq.
        
        ** Arguments **
        ref_name (bytes) ----> name of the reference
        ref_seq (DNA) -------> sequence of the entire reference (not just the
                               region of interest)
        first (int) ---------> (optional) 5'-most position of the reference
                               that lies in the region of interest
        last (int) ----------> (optional) 3'-most position of the reference
                               that lies in the region of interest
        fwd (DNA) -----------> (optional) sequence of the forward primer used
                               to amplify the region (for amplicons only)
        rev (DNA) -----------> (optional) sequence of the reverse primer used
                               to amplify the region (for amplicons only)
        
        ** Returns **
        None
        """
        if first is None:
            # If first is to be determined from the fwd primer sequence,
            # the primer is aligned to the reference, and first is set to
            # the position (primer_gap + 1) downstream of its 3' coordinate.
            # The (+ 1) makes the region start one position 3' of the primer
            # if primer_gap == 0.
            first = (1 if fwd is None else self.locate_primer(
                ref_seq, fwd)[1] + self.primer_gap + 1)
        if last is None:
            # If last is to be determined from the rev primer sequence,
            # the reverse complement of the primer is aligned to the reference,
            # and last is set to the position (primer_gap + 1) upstream of its
            # 5' coordinate. The (- 1) makes the region start one position
            # 5' of the primer if primer_gap == 0.
            last = (len(ref_seq) if rev is None else self.locate_primer(
                ref_seq, rev.rc)[0] - self.primer_gap - 1)
        super().__init__(ref_name, first, last, ref_seq)

    @staticmethod
    def locate_primer(ref_seq: DNA, primer: DNA) -> tuple[int, int]:
        """
        Return the 5' and 3' positions (1-indexed) of a primer within a
        reference sequence. The primer must occur exactly once in the
        reference, otherwise an error is raised.

        ** Arguments **
        ref_seq (DNA) -> sequence of the reference
        primer (DNA) --> sequence of the forward primer, or of the reverse
                         complement of the reverse primer
        
        ** Returns **
        pos5 (int) ----> 5'-most position of the primer aligned to the reference
        pos3 (int) ----> 3'-most position of the primer aligned to the reference
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
        return pos5, pos3


class MutationalProfile(Region):
    __slots__ = ["_sample_name"]

    def __init__(self, sample_name: str, ref_name: bytes,
                 first: int, last: int, ref_seq: DNA):
        """
        Initialize a MutationalProfile object that represents all of the reads
        from a particular sample that overlap a region of interest.
        
        ** Arguments **
        sample_name (str) ---> name of the sample
        ref_name (bytes) ----> name of the reference
        first (int) ---------> 5'-most position of the reference that lies in
                               the region of interest (1-indexed, inclusive)
        last (int) ----------> 3'-most position of the reference that lies in
                               the region of interest (1-indexed, inclusive)
        ref_seq (DNA) -------> sequence of the entire reference (not just the
                               region of interest)
        
        ** Returns **
        None
        """
        super().__init__(ref_name, first, last, ref_seq)
        self._sample_name = sample_name
    
    @property
    def sample_name(self):
        return self._sample_name

    @property
    def short_path(self):
        return os.path.join(self.sample_name,
                            self.ref_name.decode(),
                            f"{self.first}-{self.last}")
    
    def __str__(self) -> str:
        return (f"Mutational Profile of sample '{self.sample_name}' reference"
                f" '{self.ref_name}' region {self.first}-{self.last}")


class VectorIO(MutationalProfile):
    __slots__ = ["_out_dir", "_num_batches", "_num_vectors", "_checksums"]

    digest_algo = "md5"

    def __init__(self, out_dir: str, sample_name: str, ref_name: bytes,
                 first: int, last: int, ref_seq: DNA):
        """
        Initialize a VectorIO object to read and write mutation vectors.
        
        ** Arguments **
        out_dir (str) -------> path to directory in which the output files
                               will be written ({out_dir}/{sample_name}/
                               {ref_name}/{first}-{last})
        sample_name (str) ---> name of the sample
        ref_name (bytes) ----> name of the reference
        first (int) ---------> 5'-most position of the reference that lies in
                               the region of interest (1-indexed, inclusive)
        last (int) ----------> 3'-most position of the reference that lies in
                               the region of interest (1-indexed, inclusive)
        ref_seq (DNA) -------> sequence of the entire reference (not just the
                               region of interest)

        ** Returns **
        None
        """
        super().__init__(sample_name, ref_name, first, last, ref_seq)
        self._out_dir = os.path.abspath(out_dir)
        self._num_batches = 0
        self._num_vectors = 0
        self._checksums = list()
    
    @property
    def out_dir(self):
        return self._out_dir
    
    @property
    def num_batches(self):
        return self._num_batches
    
    @property
    def num_vectors(self):
        return self._num_vectors
    
    @property
    def checksums(self):
        return self._checksums

    @property
    def full_path(self):
        return os.path.join(self.out_dir, self.short_path)

    @property
    def report_file(self):
        return f"{self.full_path}_report.txt"

    @property
    def mv_dir(self):
        return self.full_path

    def get_mv_filename(self, batch_num: int):
        return os.path.join(self.mv_dir, f"{batch_num}.orc")
    
    @property
    def batch_nums(self):
        """List all of the batch numbers."""
        return range(self.num_batches)

    @property
    def mv_files(self):
        return list(map(self.get_mv_filename, self.batch_nums))

    @classmethod
    def digest_file(cls, path: str) -> str:
        """
        Compute the MD5 checksum of a file.
        
        ** Arguments **
        path (str) ---> path of the file
        
        ** Returns **
        digest (str) <- MD5 checksum of the file
        """
        with open(path, "rb") as f:
            digest = file_digest(f, cls.digest_algo).hexdigest()
        return digest


class Report(VectorIO):
    __slots__ = ["_began", "_ended"]

    # fields is a dict that maps the name of each field to its data type
    # and defines the order of the fields in the report file.
    fields = {"Sample Name": str, "Ref Name": bytes, "First": int, "Last": int,
              "Ref Seq": DNA, "Num Batches": int, "Num Vectors": int,
              "Checksums": list, "Began": datetime, "Ended": datetime,
              "Duration": float, "Speed": float}
    
    # units is a dict that defines the units of several fields that have them.
    units = {"Speed": "vec/s", "Duration": "s",
             "Checksums": VectorIO.digest_algo}

    # format of dates and times in the report file
    datetime_fmt = "on %Y-%m-%d at %H:%M:%S.%f"

    def __init__(self, out_dir: str, sample_name: str, ref_name: bytes,
                 first: int, last: int, ref_seq: DNA, num_batches: int,
                 num_vectors: int, checksums: List[str],
                 began: datetime, ended: datetime):
        """
        Initialize a Report object to
        - record information about a mutational profile and its mutation vectors
          in a human-readable plain-text format
        - facilitate validating and loading ORC files of the mutation vectors
          during subsequent steps of DREEM
        
        ** Arguments **
        out_dir (str) ---------> path to directory in which the output files
                                 will be written ({out_dir}/{sample_name}/
                                 {ref_name}/{first}-{last})
        sample_name (str) -----> name of the sample
        ref_name (bytes) ------> name of the reference
        first (int) -----------> 5'-most position of the reference that lies in
                                 the region of interest (1-indexed, inclusive)
        last (int) ------------> 3'-most position of the reference that lies in
                                 the region of interest (1-indexed, inclusive)
        ref_seq (DNA) ---------> sequence of the entire reference (not just the
                                 region of interest)
        num_batches (int) -----> number of batches in the mutational profile
        num_vectors (int) -----> number of vectors in the mutational profile
        checksums (list[str]) -> list of checksums for mutation vector files

        ** Returns **
        None
        """
        super().__init__(out_dir, sample_name, ref_name, first, last, ref_seq)
        self._num_batches = num_batches
        self._num_vectors = num_vectors
        self._checksums = checksums
        assert ended >= began
        self._began = began
        self._ended = ended
    
    @property
    def began(self):
        return self._began
    
    @property
    def ended(self):
        return self._ended
    
    @property
    def duration(self) -> float:
        """ Return duration of the computation (in seconds) """
        dt = self.ended - self.began
        return dt.seconds + dt.microseconds / 1E6
    
    @property
    def speed(self) -> float:
        """ Return speed of the computation (in vectors per second) """
        try:
            return self.num_vectors / self.duration
        except ZeroDivisionError:
            return float("inf" if self.num_vectors else "nan")

    @classmethod
    def append_unit(cls, field: str) -> str:
        """ Append the unit to the name of a field. """
        return f"{field} ({unit})" if (unit := cls.units.get(field)) else field

    @classmethod
    def remove_unit(cls, field: str) -> str:
        """ Remove the unit from a field, leaving only its name.
        NOTE: This method only works assuming that every field that has a unit
        (namely "Speed" and "Duration") contains no whitespace. This method
        will need rewriting if this assumption changes in the future. """
        return word if cls.units.get(word := field.split(" ")[0]) else field

    @classmethod
    def attr_to_field(cls, attr: str):
        """ Convert the name of an attribute of a Report object (in snake_case)
        to the name of a field in the report file (in Capital Case). """
        return cls.append_unit(" ".join(map(str.capitalize, attr.split("_"))))

    @classmethod
    def field_to_attr(cls, field: str):
        """ Convert the name of a field in the report file (in Capital Case)
        to the name of an attribute of a Report object (in snake_case). """
        return cls.remove_unit(field).replace(' ', '_').lower()

    @classmethod
    def format_val(cls, val: Any) -> str:
        """ Return a string representation of the value of a field. """
        dtype = type(val)
        if dtype is str or dtype is int or dtype is DNA:
            return val
        if dtype is bytes:
            return val.decode()
        if dtype is float:
            return round(val, 2)
        if dtype is datetime:
            return val.strftime(cls.datetime_fmt)
        if dtype is list:
            return ", ".join(val)
        raise ValueError(dtype)

    @classmethod
    def parse_valstr(cls, field: str, valstr: str) -> Any:
        """ Return the value of a field from a string representation. """
        dtype = cls.fields[cls.remove_unit(field)]
        if dtype is str or dtype is int or dtype is float:
            return dtype(valstr)
        if dtype is bytes:
            return valstr.encode()
        if dtype is DNA:
            return DNA(valstr.encode())
        if dtype is datetime:
            return datetime.strptime(valstr, cls.datetime_fmt)
        if dtype is list:
            return valstr.split(", ")
        raise ValueError(dtype)

    def save(self):
        """ Save the information in a Report to a file. """
        # Determine the maximum number of characters in a label.
        width = max(map(len, map(self.append_unit, self.fields.keys())))
        # Create a format string that pads every label to that length.
        pattern = "{field: <" + str(width) + "}\t{val}\n"
        lines: List[str] = list()
        # Iterate through all the fields to include in the report.
        for field in self.fields.keys():
            # Get the name of the Report attribute corresponding to the field.
            attr = self.field_to_attr(field)
            # Get a string representation of the value of the field.
            valstr = self.format_val(self.__getattribute__(attr))
            # Format the line of the report with the field's name and value.
            line = pattern.format(field=self.append_unit(field), val=valstr)
            lines.append(line)
        # Write the lines to the report file.
        with open(self.report_file, "w") as f:
            f.write("".join(lines))
    
    @classmethod
    def load(cls, report_file: str) -> Report:
        """ Return a Report from a saved report file. """
        # FIXME: make this function to find out_dir more elegant
        region_dir = os.path.dirname(os.path.abspath(report_file))
        ref_dir = os.path.dirname(region_dir)
        sample_dir = os.path.dirname(ref_dir)
        module_dir = os.path.dirname(sample_dir)
        out_dir = os.path.dirname(module_dir)
        # Initialize the dict of attribute values with out_dir.
        vals: Dict[str, Any] = {"out_dir": out_dir}
        with open(report_file) as f:
            for line in f:
                # Read the field and the string representation of its value.
                field, valstr = map(str.rstrip, line.split("\t"))
                # Get the name of the attribute corresponding to the field.
                attr = cls.field_to_attr(field)
                # Parse the string representation and set the attribute value.
                vals[attr] = cls.parse_valstr(field, valstr)
        # Return a new Report object from the attributes and their values.
        return cls(**vals)


class VectorWriter(VectorIO):
    __slots__ = ["_bam_file", "_parallel_reads", "_region_seq_bytes"]

    """
    Computes mutation vectors for all reads from one sample mapping to one
    region of one reference sequence.
    """
    def __init__(self, out_dir: str, bam_file: str, ref_name: bytes,
                 first: int, last: int, ref_seq: DNA, parallel_reads: bool):
        sample = sample_from_bam_file(bam_file)
        super().__init__(out_dir, sample, ref_name, first, last, ref_seq)
        self._bam_file = bam_file
        self._parallel_reads = parallel_reads
        self._region_seq_bytes = bytes(self.region_seq)

    def _comp_vector(self, rec: SamRecord):
        """
        Compute the mutation vector of one record from a SAM file.

        ** Arguments **
        rec (SamRecord) --> SAM record for which to compute a mutation vector

        ** Returns **
        muts (bytearray) <- mutation vector
        """
        assert rec.ref_name == self.ref_name
        muts = rec.vectorize(self._region_seq_bytes, self.first, self.last)
        assert muts != bytes(len(muts))
        return muts
    
    def _write_vector_batch(self, muts: np.ndarray, batch_num: int) -> str:
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
        assert batch_num >= 0
        # Data must be converted to pd.DataFrame for PyArrow to write.
        # Explicitly set copy=False to copying the mutation vectors.
        df = pd.DataFrame(data=muts, columns=self.columns, copy=False)
        mv_file = self.get_mv_filename(batch_num)
        df.to_orc(mv_file, engine="pyarrow")
        return mv_file

    def _gen_vector_batch(self, sam_viewer: SamViewer, batch_num: int,
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
        with sam_viewer as sv:
            # Use the SAM viewer to generate the mutation vectors.
            # Collect them as a single, 1-dimensional bytes object.
            mut_bytes = b"".join(map(self._comp_vector,
                                     sv.get_records(start, stop)))
        if not mut_bytes:
            raise Warning(f"{self} contained no reads.")
        # Convert the concatenated mutation vectors to a 2D NumPy array.
        muts = np.frombuffer(mut_bytes, dtype=np.byte)
        n_records, rem = divmod(len(muts), self.length)
        assert rem == 0
        muts.resize((n_records, self.length))
        # Write the mutation vectors to a file and compute its checksum.
        mv_file = self._write_vector_batch(muts, batch_num)
        checksum = self.digest_file(mv_file)
        return n_records, checksum

    def _gen_vectors(self):
        with SamViewer(self._bam_file, self.ref_name, self.first, self.last,
                       self.spanning) as sv:
            batch_size = max(1, DEFAULT_BATCH_SIZE // self.length)
            indexes = list(sv.get_batch_indexes(batch_size))
            starts = indexes[:-1]
            stops = indexes[1:]
            self._num_batches = len(starts)
            assert self.num_batches == len(stops)
            svs = [SamViewer(sv.sam_subset, self.ref_name, self.first,
                             self.last, self.spanning, make=False, remove=False)
                   for _ in self.batch_nums]
            args = list(zip(svs, self.batch_nums, starts, stops))
            if self._parallel_reads:
                n_procs = max(1, min(NUM_PROCESSES, self.num_batches))
                with Pool(n_procs, maxtasksperchild=1) as pool:
                    results = pool.starmap(self._gen_vector_batch, args,
                                           chunksize=1)
            else:
                results = list(itertools.starmap(self._gen_vector_batch, args))
            assert len(results) == self.num_batches
            assert self.num_vectors == 0

            assert len(self.checksums) == 0
            for num_vectors, checksum in results:
                self._num_vectors += num_vectors
                self._checksums.append(checksum)

    def _write_report(self, t_start: datetime, t_end: datetime):
        Report(self.out_dir, self.sample_name, self.ref_name, self.first,
               self.last, self.ref_seq, self.num_batches, self.num_vectors,
               self.checksums, t_start, t_end).save()
    
    def gen_vectors(self):
        if not (all(os.path.isfile(f) for f in self.mv_files)
                and os.path.isfile(self.report_file)):
            os.makedirs(self.mv_dir, exist_ok=True)
            print(f"{self}: computing vectors")
            t_start = datetime.now()
            self._gen_vectors()
            t_end = datetime.now()
            print(f"{self}: writing report")
            self._write_report(t_start, t_end)
            print(f"{self}: finished")


class VectorWriterSpawner(object):
    def __init__(self,
                 out_dir: str,
                 fasta: str,
                 bam_files: List[str],
                 coords: List[Tuple[bytes, int, int]],
                 primers: List[Tuple[bytes, DNA, DNA]],
                 fill: bool,
                 parallel: str):
        self.out_dir = out_dir
        self.bam_files = bam_files
        self.ref_file = fasta
        self.coords = coords
        self.primers = primers
        self.fill = fill
        if parallel == "auto":
            parallel = ("reads" if self.num_samples == self.num_regions == 1
                        else "profiles")
        if parallel == "profiles":
            self.parallel_profiles = True
            self.parallel_reads = False
        elif parallel == "reads":
            self.parallel_profiles = False
            self.parallel_reads = True
        elif parallel == "off":
            self.parallel_profiles = False
            self.parallel_reads = False
        else:
            raise ValueError(f"Invalid value for parallel: '{parallel}'")
    
    @property
    def num_samples(self):
        return len(samples_from_bam_files(self.bam_files))
    
    @property
    def num_regions(self):
        return len(self.coords) + len(self.primers)
    
    @cached_property
    def ref_seqs(self):
        seqs = dict(FastaParser(self.ref_file).parse())
        if not seqs:
            raise ValueError(f"'{self.ref_file}' contained no sequences")
        return seqs
    
    @property
    def regions(self):
        regions: Dict[bytes, List[RegionFinder]] = defaultdict(list)

        def add_region(region: RegionFinder):
            if any(region == other for other in regions[region.ref_name]):
                raise ValueError(f"Duplicate region: {region.ref_coords}")
            regions[region.ref_name].append(region)
        
        def ref_to_seq(ref_seqs: Dict[bytes, DNA], ref: bytes):
            try:
                return ref_seqs[ref]
            except KeyError:
                raise ValueError(f"No reference named '{ref.decode()}'")

        for ref, first, last in self.coords:
            add_region(RegionFinder(ref, ref_to_seq(self.ref_seqs, ref),
                                    first=first, last=last))
        for ref, fwd, rev in self.primers:
            add_region(RegionFinder(ref, ref_to_seq(self.ref_seqs, ref),
                                    fwd=fwd, rev=rev))
        if self.fill:
            for ref, seq in self.ref_seqs.items():
                if ref not in regions:
                    add_region(RegionFinder(ref, seq))
        return regions

    @property
    def writers(self):
        no_writers = True
        for bam_file in self.bam_files:
            ref_name = ref_from_bam_file(bam_file)
            for region in self.regions[ref_name]:
                assert region.ref_name == ref_name
                no_writers = False
                yield VectorWriter(self.out_dir, bam_file, ref_name,
                                   region.first, region.last, region.ref_seq,
                                   self.parallel_reads)
        if no_writers:
            raise ValueError("No samples and/or regions were given.")

    def gen_mut_profiles(self, processes: int = 0):
        if self.parallel_profiles:
            writers = list(self.writers)
            with Pool(processes if processes
                      else min(NUM_PROCESSES, len(writers)),
                      maxtasksperchild=1) as pool:
                pool.map(VectorWriter.gen_vectors, writers,
                         chunksize=1)
        else:
            for writer in self.writers:
                writer.gen_vectors()


'''
class VectorReader(VectorIO):
    @property
    def shape(self):
        return self.num_vectors, self.length

    def run_checksums(self):
        if len(mv_files := self.mv_files) != len(self.checksums):
            raise ValueError(f"Got {len(mv_files)} files but "
                             f"{len(self.checksums)} checksums")
        for mv_file, checksum in zip(mv_files, self.checksums):
            digest = self.digest_file(mv_file)
            if digest != checksum:
                raise ValueError(f"Hex digest of {mv_file} ({digest}) "
                                 f"did not match checksum ({checksum})")

    @property
    def vectors(self):
        self.run_checksums()
        # FIXME: I suspect there are more efficient ways to load these files
        mvs = np.full(self.shape, fill_value=BLANK, dtype=bytes)
        row = 0
        for mv_file in self.mv_files:
            data = orc.read_table(mv_file).to_pandas().values
            n_vectors, n_cols = data.shape
            if n_cols != self.length:
                raise ValueError(f"Expected DataFrame with {self.length}"
                                 f" columns but got {n_cols} columns.")
            mvs[row: (row := row + n_vectors)] = data
        if row != self.num_vectors:
            raise ValueError(f"Expected DataFrame with {self.num_vectors}"
                             f" rows but got {row} rows.")
        vectors = pd.DataFrame(data=mvs, columns=self.positions)
        return VectorSet(self.sample_name, self.ref_name, self.first,
                         self.last, self.ref_seq, vectors)

    @classmethod
    def load(cls, report_file):
        rep = Report.load(report_file)
        return cls(rep.out_dir, rep.sample_name, rep.ref_name, rep.first,
                   rep.last, rep.ref_seq, rep.num_batches, rep.num_vectors,
                   rep.checksums)


class VectorSet(MutationalProfile):
    __slots__ = ["_vectors", "_cover_count", "_mm_count", "_del_count"]

    color_dict = {"A": "red", "C": "blue", "G": "orange", "T": "green"}

    def __init__(self, sample_name: str, ref_name: str, first: int, last: int,
                 ref_seq: DNA, vectors: pd.DataFrame):
        super().__init__(sample_name, ref_name, first, last, ref_seq)
        if (vectors.shape[1] != self.length
                or (vectors.columns != self.positions).any()):
            raise ValueError("Columns of vectors do not match positions.")
        self._vectors = vectors
        self._cover_count = None
        self._mm_count = None
        self._del_count = None

    @staticmethod
    def series(f):
        def wrapper(self):
            return pd.Series(f(self), index=self.positions)
        return wrapper

    @property
    def vectors(self):
        return self._vectors

    @property
    def num_vectors(self):
        return self.vectors.shape[0]

    @property
    def mismatch_frac(self):
        return self.mismatch_count / self.coverage_count

    @property
    @series
    def deletion_count(self):
        if self._del_count is None:
            self._del_count = (self.vectors.values == DELET).sum(axis=0)
        return self._del_count

    @property
    def deletion_frac(self):
        return self.deletion_count / self.coverage_count

    @property
    def mutation_count(self):
        return self.mismatch_count + self.deletion_count

    @property
    def mutation_frac(self):
        return self.mutation_count / self.coverage_count

    @property
    def colors(self):
        return [self.color_dict[chr(base)] for base in self.region_seq]

    def plot_muts(self):
        fig, ax = plt.subplots()
        data = self.mutation_frac
        ax.bar(data.index, data, color=self.colors)
        plt.show()
'''
