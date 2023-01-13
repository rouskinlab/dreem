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
from dreem.util.seq import FastaParser
from dreem.util.path import BasePath, FastaInPath, XamInPath, MutVectorBatchPath, MutVectorBatchSegment, MutVectorReportPath, MutVectorReportSegment, MOD_VEC, TEMP_DIR, OUTPUT_DIR, RegionSegment, RegionPath
from dreem.util.seq import DNA
from dreem.vector.samview import SamViewer
from dreem.vector.vector import SamRecord


DEFAULT_BATCH_SIZE = 33_554_432  # 2^25 bytes ≈ 33.6 Mb


class Region(object):
    __slots__ = ["first", "last", "ref_seq", "ref_name"]

    def __init__(self, ref_name: str, first: int, last: int, ref_seq: DNA):
        """
        Initialize a Region object to represent a region of interest (between a
        first and a last position) in a reference sequence.
        
        ** Arguments **
        ref_name (str) ----> name of the reference
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
        self.first = first
        self.last = last
        self.ref_seq = ref_seq
        self.ref_name = ref_name
    
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
    def ref_coords(self) -> Tuple[str, int, int]:
        """ Return the name of the reference and the first and last positions
        of the region of interest; for equality testing and hashing. """
        return self.ref_name, self.first, self.last
    
    @property
    def columns(self) -> List[str]:
        """ Return a list of the bases and positions in the region of interest,
        each of the form '{base}{position}' (e.g. ['G13', 'C14', 'A15']). """
        # FIXME: add a column for the read names
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
            first = (1 if fwd is None else self.locate_primer(
                ref_seq, fwd)[1] + self.primer_offset)
        if last is None:
            # If last is to be determined from the rev primer sequence,
            # the reverse complement of the primer is aligned to the reference,
            # and last is set to the position (primer_gap + 1) upstream of its
            # 5' coordinate.
            last = (len(ref_seq) if rev is None else self.locate_primer(
                ref_seq, rev.rc)[0] - self.primer_offset)
        super().__init__(ref_name, first, last, ref_seq)

    @property
    def primer_offset(self):
        return self.primer_gap + 1

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
    __slots__ = ["sample_name"]

    def __init__(self, sample_name: str, ref_name: str,
                 first: int, last: int, ref_seq: DNA):
        """
        Initialize a MutationalProfile object that represents all of the reads
        from a particular sample that overlap a region of interest.
        
        ** Arguments **
        sample_name (str) ---> name of the sample
        ref_name (str) ----> name of the reference
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
        self.sample_name = sample_name
    
    def __str__(self) -> str:
        return (f"Mutational Profile of sample '{self.sample_name}' reference"
                f" '{self.ref_name}' region {self.first}-{self.last}")


class VectorIO(MutationalProfile):
    __slots__ = ["base_path", "num_batches", "num_vectors", "checksums"]

    digest_algo = "md5"

    def __init__(self, base_path: BasePath, sample_name: str, ref_name: str,
                 first: int, last: int, ref_seq: DNA):
        """
        Initialize a VectorIO object to read and write mutation vectors.
        
        ** Arguments **
        base_path (str) -------> path to directory in which the output files
                               will be written ({base_path}/{sample_name}/
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
        self.base_path = base_path
        self.num_batches = 0
        self.num_vectors = 0
        self.checksums = list()
    
    @property
    def args(self):
        return (self.base_path.path,
                OUTPUT_DIR,
                MOD_VEC,
                self.sample_name,
                self.ref_name,
                RegionSegment.format(self.first, self.last))

    @property
    def report_path(self):
        mv_report = MutVectorReportSegment.format(self.first, self.last)
        args = self.args[:-1] + (mv_report,)
        return MutVectorReportPath(*args)

    @property
    def batch_dir(self):
        return RegionPath(*self.args)

    def get_mv_batch_path(self, batch_num: int):
        mv_batch = MutVectorBatchSegment.format(batch_num)
        args = self.args + (mv_batch,)
        return MutVectorBatchPath(*args)
    
    @property
    def batch_nums(self):
        """ List all of the batch numbers. """
        return range(self.num_batches)

    @property
    def mv_batch_paths(self):
        return list(map(self.get_mv_batch_path, self.batch_nums))

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
    __slots__ = ["began", "ended"]

    # fields is a dict that maps the name of each field to its data type
    # and defines the order of the fields in the report file.
    fields = {"Sample Name": str, "Ref Name": str, "First": int, "Last": int,
              "Ref Seq": DNA, "Num Batches": int, "Num Vectors": int,
              "Checksums": list, "Began": datetime, "Ended": datetime,
              "Duration": float, "Speed": float}
    
    # units is a dict that defines the units of several fields that have them.
    units = {"Speed": "vec/s", "Duration": "s",
             "Checksums": VectorIO.digest_algo}

    # format of dates and times in the report file
    datetime_fmt = "on %Y-%m-%d at %H:%M:%S.%f"

    def __init__(self, base_path: str, sample_name: str, ref_name: str,
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
        base_path (str) ---------> path to directory in which the output files
                                 will be written ({base_path}/{sample_name}/
                                 {ref_name}/{first}-{last})
        sample_name (str) -----> name of the sample
        ref_name (str) ------> name of the reference
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
        super().__init__(base_path, sample_name, ref_name, first, last, ref_seq)
        self.num_batches = num_batches
        self.num_vectors = num_vectors
        self.checksums = checksums
        assert ended >= began
        self.began = began
        self.ended = ended
    
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
        with open(self.report_path.path, "w") as f:
            f.write("".join(lines))
    
    @classmethod
    def load(cls, report_file: str) -> Report:
        """ Return a Report from a saved report file. """
        # FIXME: make this function to find base_path more elegant
        region_dir = os.path.dirname(os.path.abspath(report_file))
        ref_dir = os.path.dirname(region_dir)
        sample_dir = os.path.dirname(ref_dir)
        module_dir = os.path.dirname(sample_dir)
        base_path = os.path.dirname(module_dir)
        # Initialize the dict of attribute values with base_path.
        vals: Dict[str, Any] = {"base_path": base_path}
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
    __slots__ = ["bam_path", "parallel_reads", "region_seqb"]

    """
    Computes mutation vectors for all reads from one sample mapping to one
    region of one reference sequence.
    """
    def __init__(self, base_path: BasePath, bam_path: XamInPath, ref_name: str,
                 first: int, last: int, ref_seq: DNA, parallel_reads: bool):
        sample = bam_path.sample.name
        super().__init__(base_path, sample, ref_name, first, last, ref_seq)
        self.bam_path = bam_path
        self.parallel_reads = parallel_reads
        self.region_seqb = bytes(self.region_seq)

    def _write_batch(self, read_names: Tuple[str], muts: Tuple[bytearray],
                     batch_num: int) -> Tuple[pathlib.Path, int]:
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
        # Convert the concatenated mutation vectors to a 2D NumPy array.
        muts_array = np.frombuffer(b"".join(muts), dtype=np.byte)
        n_records = len(muts_array) // self.length
        muts_array.resize((n_records, self.length))
        # Data must be converted to pd.DataFrame for PyArrow to write.
        # Explicitly set copy=False to copying the mutation vectors.
        df = pd.DataFrame(data=muts_array, index=read_names,
                          columns=self.columns, copy=False)
        mv_file = self.get_mv_batch_path(batch_num).path
        df.to_orc(mv_file, engine="pyarrow")
        return mv_file, n_records

    def _write_report(self, t_start: datetime, t_end: datetime):
        Report(self.base_path, self.sample_name, self.ref_name, self.first,
               self.last, self.ref_seq, self.num_batches, self.num_vectors,
               self.checksums, t_start, t_end).save()

    def _vectorize_record(self, rec: SamRecord):
        """
        Compute the mutation vector of one record from a SAM file.

        ** Arguments **
        rec (SamRecord) --> SAM record for which to compute a mutation vector

        ** Returns **
        muts (bytearray) <- mutation vector
        """
        assert rec.ref_name == self.ref_name
        muts = rec.vectorize(self.region_seqb, self.first, self.last)
        assert any(muts)
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
            raise Warning(f"{self} contained no reads.")
        # Write the mutation vectors to a file and compute its checksum.
        mv_file, n_records = self._write_batch(read_names, muts, batch_num)
        checksum = self.digest_file(mv_file)
        return n_records, checksum

    def _vectorize_sam(self):
        with SamViewer(self.base_path, self.bam_path, self.ref_name,
                       self.first, self.last, self.spanning) as sv:
            batch_size = max(1, DEFAULT_BATCH_SIZE // self.length)
            indexes = list(sv.get_batch_indexes(batch_size))
            starts = indexes[:-1]
            stops = indexes[1:]
            self.num_batches = len(starts)
            assert self.num_batches == len(stops)
            svs = [SamViewer(self.base_path, sv.sam_path, self.ref_name,
                             self.first, self.last, self.spanning, owner=False)
                   for _ in self.batch_nums]
            args = list(zip(svs, self.batch_nums, starts, stops))
            if self.parallel_reads:
                n_procs = max(1, min(NUM_PROCESSES, self.num_batches))
                with Pool(n_procs, maxtasksperchild=1) as pool:
                    results = pool.starmap(self._vectorize_batch, args,
                                           chunksize=1)
            else:
                results = list(itertools.starmap(self._vectorize_batch, args))
            assert len(results) == self.num_batches
            assert self.num_vectors == 0
            assert len(self.checksums) == 0
            for num_vectors, checksum in results:
                self.num_vectors += num_vectors
                self.checksums.append(checksum)
    
    def vectorize(self):
        if not (all(f.path.is_file() for f in self.mv_batch_paths)
                and self.report_path.path.is_file()):
            self.batch_dir.path.mkdir(parents=True, exist_ok=True)
            print(f"{self}: computing vectors")
            t_start = datetime.now()
            self._vectorize_sam()
            t_end = datetime.now()
            print(f"{self}: writing report")
            self._write_report(t_start, t_end)
            print(f"{self}: finished")


class VectorWriterSpawner(object):
    def __init__(self,
                 base_dir: str,
                 fasta: str,
                 bam_files: List[str],
                 coords: List[Tuple[str, int, int]],
                 primers: List[Tuple[str, DNA, DNA]],
                 fill: bool,
                 parallel: str):
        self.base_path = BasePath.parse(base_dir)
        self.bam_paths = list(map(XamInPath.parse, bam_files))
        self.ref_path = FastaInPath.parse(fasta)
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
    def bams_per_sample(self):
        samples: Dict[str, List[XamInPath]] = defaultdict(list)
        for bam in self.bam_paths:
            samples[bam.sample].append(bam)
        return samples
    
    @property
    def samples(self):
        return set(self.bams_per_sample.keys())
    
    @property
    def num_samples(self):
        return len(self.samples)
    
    @property
    def num_regions(self):
        return len(self.coords) + len(self.primers)
    
    @cached_property
    def ref_seqs(self):
        seqs = dict(FastaParser(self.ref_path.path).parse())
        if not seqs:
            raise ValueError(f"'{self.ref_path}' contained no sequences")
        return seqs
    
    @cached_property
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
        for bam in self.bam_paths:
            print(bam, bam.xam._segment)
            ref_name = bam.xam.name
            for region in self.regions[ref_name]:
                assert region.ref_name == ref_name
                yield VectorWriter(self.base_path, bam, ref_name,
                                   region.first, region.last, region.ref_seq,
                                   self.parallel_reads)

    def profile(self, processes: int = 0):
        writers = list(self.writers)
        if not writers:
            raise ValueError("No samples and/or regions were given.")
        if self.parallel_profiles:
            with Pool(processes if processes
                      else min(NUM_PROCESSES, len(writers)),
                      maxtasksperchild=1) as pool:
                pool.map(VectorWriter.vectorize, writers,
                         chunksize=1)
        else:
            for writer in writers:
                writer.vectorize()


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
        return cls(rep.base_path, rep.sample_name, rep.ref_name, rep.first,
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
