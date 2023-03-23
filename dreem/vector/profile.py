from __future__ import annotations

from collections import defaultdict, namedtuple
from datetime import datetime, timedelta
from functools import cache, cached_property, partial
from inspect import isclass
import itertools
import json
import logging
from multiprocessing import Pool
import re
import sys
import time
from typing import Any, Sequence

import numpy as np
import pandas as pd

from ..util import path
from ..util.seq import (BLANK, MATCH, DELET, INS_5, INS_3,
                        SUB_A, SUB_C, SUB_G, SUB_T, EVERY,
                        BASES, DNA, FastaParser)
from ..util.util import digest_file, get_num_parallel
from ..vector.samread import SamReader
from ..vector.vector import vectorize_line, vectorize_pair, VectorError

SectionTuple = namedtuple("PrimerTuple", ["pos5", "pos3"])


def mib_to_bytes(batch_size: float):
    """
    Return the number of bytes per batch of a given size in mebibytes.

    Parameters
    ----------
    batch_size: float
        Size of the batch in mebibytes (1 MiB = 2^20 bytes)

    Return
    ------
    int
        Number of bytes per batch, to the nearest integer
    """
    return round(batch_size * 1048576)  # 1048576 = 2^20


class Section(object):
    """
    Represent a section of a reference sequence between two coordinates.

    Attributes
    ----------
    ref: str
        Name of the reference sequence
    end5: int (1 ≤ end5 ≤ end3)
        Coordinate of the reference sequence at which the section's
        5' end is located (1-indexed)
    end3: int (end5 ≤ end3 ≤ len(ref_seq))
        Coordinate of the reference sequence at which the section's
        3' end is located (1-indexed; end3 itself is included)
    seq: DNA
        Sequence of the section between end5 and end3 (inclusive)
    isfullref: bool
        Whether the section sequence is the full reference sequence

    Examples
    --------
    >>> seq = DNA(b"CATCTGGA")
    >>> name = "example"
    >>> section = Section(ref_seq=seq, ref=name, end5=1, end3=8)
    >>> assert section.seq == seq
    >>> section = Section(ref_seq=seq, ref=name, end5=1, end3=-1)
    >>> assert section.seq == seq
    >>> section = Section(ref_seq=seq, ref=name, end5=-8, end3=8)
    >>> assert section.seq == seq
    >>> section = Section(ref_seq=seq, ref=name, end5=3, end3=7)
    >>> assert section.seq == DNA(b"TCTGG")
    >>> section = Section(ref_seq=seq, ref=name, end5=-5, end3=-3)
    >>> assert section.seq == DNA(b"CTG")
    >>> try:
    ...     section = Section(ref_seq=seq, ref=name, end5=-9, end3=5)
    ...     assert False, "Failed to catch end5 < -len(ref_seq)"
    ... except ValueError:
    ...     pass
    >>> try:
    ...     section = Section(ref_seq=seq, ref=name, end5=6, end3=5)
    ...     assert False, "Failed to catch end3 < end5"
    ... except ValueError:
    ...     pass
    >>> try:
    ...     section = Section(ref_seq=seq, ref=name, end5=1, end3=9)
    ...     assert False, "Failed to catch end3 > len(ref_seq)"
    ... except ValueError:
    ...     pass
    """

    def __init__(self, /, *,
                 ref: str,
                 end5: int,
                 end3: int,
                 ref_seq: DNA | None = None,
                 sect_seq: DNA | None = None,
                 isfullref: bool = False):
        """
        Parameters
        ----------
        ref: str
            Name of the reference sequence
        end5: int (-len(ref_seq) ≤ end5 ≤ len(ref_seq); end5 ≠ 0)
            Coordinate of the reference sequence at which the 5' end of
            the section is located. If positive, number the coordinates
            of the reference sequence 1, 2, ... starting at the 5' end
            (i.e. 1-based indexing). If negative, number the coordinates
            of the reference sequence -1, -2, ... starting at the 3' end
            (i.e. 1-based indexing from the other side), then convert to
            the corresponding (positive) 1-based index from the 5' end
        end3: int (-len(ref_seq) ≤ end3 ≤ len(ref_seq); end3 ≠ 0)
            Coordinate of the reference sequence at which the section's
            3' end is located. Follows the same coordinate numbering
            convention as end5
        ref_seq: DNA | None
            Sequence of the entire reference (must provide either
            ```ref_seq``` or ```sect_seq```, but not both)
        sect_seq: DNA | None
            Sequence of the section only (must provide either
            ```ref_seq``` or ```sect_seq```, but not both)
        """
        self.ref = ref
        if ref_seq and sect_seq:
            raise ValueError("Cannot give both ref_seq and sect_seq")
        if ref_seq:
            if end5 < 0:
                # Compute the corresponding positive coordinate.
                end5 += len(ref_seq) + 1
            if end3 < 0:
                # Compute the corresponding positive coordinate.
                end3 += len(ref_seq) + 1
            self.end5 = end5
            self.end3 = end3
            if not 1 <= end5 <= end3 <= len(ref_seq):
                raise ValueError("Must have 1 ≤ end5 ≤ end3 ≤ len(ref_seq), "
                                 f"but got end5 = {end5}, end3 = {end3}, and "
                                 f"len(ref_seq) = {len(ref_seq)}")
            self.seq = ref_seq[end5 - 1: end3]
            self.isfullref = self.seq == ref_seq
        elif sect_seq:
            self.end5 = end5
            self.end3 = end3
            if not 1 <= end5 <= end3:
                raise ValueError("Must have 1 ≤ end5 ≤ end3 ≤ len(ref_seq), "
                                 f"but got end5 = {end5}, end3 = {end3}")
            if not self.length == len(sect_seq):
                raise ValueError(f"Calculated length of {self.length} from "
                                 f"end5 = {end5} and end3 = {end3}, but got "
                                 f"sect_seq of length {len(sect_seq)}")
            self.seq = sect_seq
            self.isfullref = isfullref
        else:
            raise ValueError("Must give either ref_seq or sect_seq")

    @property
    def path_fields(self):
        return {"ref": self.ref, "end5": self.end5, "end3": self.end3}

    @property
    def report_fields(self):
        return {"ref": self.ref, "end5": self.end5, "end3": self.end3,
                "isfullref": self.isfullref}

    @property
    def length(self):
        """ Return the length of the section of interest. """
        return self.end3 - self.end5 + 1

    @cached_property
    def positions(self):
        """ Return all positions in the section of interest. """
        return np.arange(self.end5, self.end3 + 1)

    def subseq(self, positions: Sequence[int] | None):
        """ Return a subset of the sequence at the given positions, or
        the entire sequence if positions is None. """
        if positions is None:
            # Return the entire sequence if no positions are selected.
            return self.seq
        n_pos = len(positions)
        if n_pos == 0:
            raise ValueError("Positions is an empty sequence")
        pos5, pos3 = positions[0], positions[-1]
        if n_pos != pos3 - pos5 + 1:
            raise ValueError(
                "Positions must be a sequence of monotonically increasing "
                f"consecutive integers, but got {positions}")
        if pos5 < self.end5:
            raise ValueError(f"5' end ({pos5}) out of bounds for {self}")
        if pos3 > self.end3:
            raise ValueError(f"3' end ({pos3}) out of bounds for {self}")
        return self.seq[pos5 - self.end5: pos3 - self.end5 + 1]

    @property
    def ref_coords(self):
        """ Return the name of the reference and the 5' and 3' positions
        of the section of interest; for hashing and equality test_input. """
        return self.ref, self.end5, self.end3

    @staticmethod
    def seq_pos_to_cols(seq: bytes, positions: Sequence[int]):
        """ Convert sequence and positions to column names. Each column
        name is a string of the base followed by the position. """
        # Use chr(base) instead of base.decode() because base is an int.
        return [f"{chr(base)}{pos}" for base, pos in zip(seq, positions,
                                                         strict=True)]

    @staticmethod
    def cols_to_seq_pos(columns: list[str]):
        """ Convert column names to sequence and positions. Each column
        name is a string of the base followed by the position. """
        # Regex pattern "^([ACGT])([0-9]+)$" finds the base and position
        # in each column name.
        pattern = re.compile(f"^([{BASES.decode()}])([0-9]+)$")
        # Match each column name using the pattern.
        matches = list(map(pattern.match, columns))
        try:
            # Obtain the two groups (base and position) from each match
            # and unzip them into two tuples.
            bases, pos_strs = zip(*map(re.Match.groups, matches))
        except TypeError:
            # TypeError is raised if any match is None, which happens if
            # a column fails to match the pattern.
            invalid_cols = [col for col, match in zip(columns, matches,
                                                      strict=True)
                            if match is None]
            raise ValueError(f"Invalid columns: {invalid_cols}")
        # Join the tuple of bases into a DNA sequence.
        seq = DNA("".join(bases).encode())
        # Cast the tuple of strings of positions into an integer array.
        positions = np.array(list(map(int, pos_strs)))
        return seq, positions

    @cached_property
    def columns(self):
        """ Return the column names of the section. """
        return self.seq_pos_to_cols(self.seq, self.positions)

    @property
    def tag(self):
        """ Return a hashable identifier for the section. """
        return self.ref, self.end5, self.end3

    @property
    def section(self):
        return f"{self.end5}-{self.end3}"

    def __str__(self):
        return f"{self.ref}:{self.section}"


class SectionFinder(Section):
    """
    The 5' and 3' ends of a section can be given explicitly as integers, but if
    the sample is of an amplicon (i.e. generated by RT-PCR using site-specific
    primers), then it is often more convenient to enter the sequences of the
    PCR primers and have the software determine the coordinates. SectionFinder
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

    def __init__(self, /, *, ref_seq: DNA | None, ref: str, primer_gap: int,
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
            primer to exclude from the section. Coordinates within 1 - 2
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
            raise ValueError(f"No sequence for reference named '{ref}'. Check "
                             "that you gave the right reference sequence file "
                             "and spelled the name of the reference correctly.")
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
    def locate(ref_seq: DNA, primer: DNA) -> SectionTuple:
        """
        Return the 5' and 3' positions (1-indexed) of a primer within a
        reference sequence. The primer must occur exactly once in the
        reference, otherwise an error is raised.

        Parameters
        ----------
        ref_seq: DNA
            Sequence of the entire reference (not just the section of interest)
        primer: DNA
            Sequence of the forward PCR primer or of the reverse complement of
            the reverse PCR primer
        
        Returns
        -------
        SectionTuple
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
        return SectionTuple(pos5, pos3)


class MutationalProfile(Section):
    """
    Represent all reads from one sample that overlap a section of interest in a
    reference sequence.

    Fields
    ------
    sample: str
        Name of the sample
    """

    def __init__(self, /, *, sample: str, **kwargs):
        super().__init__(**kwargs)
        self.sample = sample

    @property
    def path_fields(self):
        return {"module": path.Module.VECTOR.value,
                "sample": self.sample,
                **super().path_fields}

    @property
    def report_fields(self):
        return {"sample": self.sample,
                **super().report_fields}

    @property
    def tag(self):
        return tuple([self.sample, *super().tag])

    def get_batch_dir(self, out_dir: str):
        """ Directory in which all batches of mutation vectors are
        written. """
        return path.SectionOutDirPath(top=out_dir, **self.path_fields)

    def get_mv_batch_path(self, out_dir: str, batch: int):
        """ File in which one batch of mutation vectors is written. """
        return path.MutVectorBatchFilePath(top=out_dir,
                                           **self.path_fields,
                                           batch=batch,
                                           ext=".orc")

    def __str__(self):
        return f"{self.sample}@{super().__str__()}"


class VectorWriter(MutationalProfile):
    """
    Compute mutation vectors for all reads from one sample mapping to one
    section of one reference sequence.
    """

    def __init__(self, /, bam_file: Any, **kwargs):
        bam_path = path.OneRefAlignmentInFilePath.parse(bam_file)
        super().__init__(sample=bam_path.sample,
                         ref=bam_path.ref,
                         **kwargs)
        self.bam_path = bam_path
        self.byteseq = bytes(self.seq)

    def get_report_path(self, out_dir: str):
        return path.MutVectorReportFilePath(top=out_dir,
                                            **self.path_fields,
                                            ext=".json")

    def outputs_valid(self, /, out_dir: str):
        """ Return whether the report file exists and, if so, whether
        all batch files of mutation vectors listed in the report exist
        and match the checksums recorded in the report. """
        try:
            VectorReport.load(self.get_report_path(out_dir).path)
        except (FileNotFoundError, ValueError):
            return False
        else:
            return True

    def _write_report(self, /, *, out_dir: str, **kwargs):
        report = VectorReport(seq=self.seq, **self.report_fields, **kwargs)
        report_path = report.save(out_dir)
        return report_path

    def _write_batch(self, /,
                     out_dir: str,
                     batch_num: int,
                     read_names: list[str],
                     vectors: tuple[bytearray, ...],
                     n_passed: int):
        """ Write a batch of mutation vectors to an ORC file. """
        # Process the mutation vectors into a 2D NumPy array.
        array = np.frombuffer(b"".join(vectors), dtype=np.byte)
        try:
            array.shape = (n_passed, self.length)
        except ValueError:
            return None
        # Data must be converted to pd.DataFrame for PyArrow to write.
        # Set copy=False to prevent copying the mutation vectors.
        mut_frame = pd.DataFrame(data=array,
                                 index=read_names,
                                 columns=self.columns,
                                 copy=False)
        mv_file = self.get_mv_batch_path(out_dir, batch_num)
        mv_file.path.parent.mkdir(parents=True, exist_ok=True)
        mut_frame.to_orc(mv_file.path, index=True, engine="pyarrow")
        return mv_file

    def _vectorize_record(self, /,
                          read_name: bytes, line1: bytes, line2: bytes, *,
                          min_qual: int, ambid: bool):
        """ Compute the mutation vector of a record in a SAM file. """
        try:
            muts = bytearray(self.length)
            if line2:
                # Passing byteseq causes
                vectorize_pair(line1, line2, muts, self.byteseq, self.length,
                               self.end5, self.ref, min_qual, ambid)
            else:
                vectorize_line(line1, muts, self.byteseq, self.length,
                               self.end5, self.ref, min_qual, ambid)
            if not any(muts):
                raise VectorError(f"Vector was blank")
        except VectorError as error:
            logging.error(
                f"Read '{read_name.decode()}' failed to vectorize: {error}")
            return "", b""
        return read_name, muts

    def _vectorize_records(self, /, *,
                           reader: SamReader,
                           start: int,
                           stop: int,
                           strict_pairs: bool,
                           phred_enc: int,
                           min_phred: int,
                           ambid: bool):
        """ Generate a batch of mutation vectors, and return them along
        with the name of every vector. """
        if stop > start:
            with reader as reading:
                # Use the SAM reader to generate the mutation vectors.
                # Collect them as a single, 1-dimensional bytes object.
                vectorize_record = partial(self._vectorize_record,
                                           min_qual=get_min_qual(min_phred,
                                                                 phred_enc),
                                           ambid=ambid)
                iter_records = reading.get_records(start, stop, strict_pairs)
                read_names, muts = zip(*itertools.starmap(vectorize_record,
                                                          iter_records))
                # For every read for which creating a mutation vector
                # failed, an empty string was returned as the read name
                # and an empty bytearray as the mutation vector. The
                # empty read names must be filtered out, while the empty
                # mutation vectors will not cause problems because,
                # being of length zero, they will effectively disappear
                # when all the vectors are concatenated into a 1D array.
                read_names = list(filter(None, read_names))
        else:
            logging.warning(f"Stop position ({stop}) is not after "
                            f"start position ({start})")
            read_names, muts = [], ()
        return read_names, muts

    def _vectorize_batch(self, /,
                         batch_num: int,
                         reader: SamReader,
                         start: int,
                         stop: int, *,
                         out_dir: str,
                         **kwargs):
        """ Compute mutation vectors for every SAM record in one batch,
        write the vectors to a batch file, and return its MD5 checksum
        and the number of vectors. """
        # Compute the read names and mutation vectors.
        read_names, muts = self._vectorize_records(reader=reader,
                                                   start=start,
                                                   stop=stop,
                                                   **kwargs)
        # Compute the number of reads that passed and failed.
        n_reads = len(muts)
        n_pass = len(read_names)
        n_fail = n_reads - n_pass
        if n_pass == 0:
            logging.critical(f"Failed to assemble batch {batch_num} of {self}")
            return n_pass, n_fail, "FAIL"
        # Write the names and vectors to a file.
        b_file = self._write_batch(out_dir, batch_num, read_names, muts, n_pass)
        if b_file is None:
            logging.critical(f"Failed to assemble batch {batch_num} of {self}")
            return n_pass, n_fail, "FAIL"
        # Compute the MD5 checksum of the file.
        checksum = digest_file(b_file.path)
        return n_pass, n_fail, checksum

    def _vectorize_bam(self, /, *,
                       out_dir: str,
                       temp_dir: str,
                       batch_size: int,
                       n_procs: int,
                       save_temp: bool,
                       resume: bool,
                       **kwargs):
        """ Compute a mutation vector for every record in a BAM file,
        split among one or more batches. For each batch, write the
        vectors to one batch file, and compute its checksum. """
        # Open the primary SAM file reader to write the subset of SAM
        # records to a temporary SAM file and determine the number and
        # start/stop indexes of each batch of records in the file.
        # The SAM file will remain open until exiting the with block.
        with SamReader(xam_path=self.bam_path,
                       end5=self.end5,
                       end3=self.end3,
                       isfullref=self.isfullref,
                       temp_dir=temp_dir,
                       n_procs=n_procs,
                       save_temp=save_temp,
                       resume=resume,
                       owner=True) as reader:
            # Use integer division to round down the number of vectors
            # per batch to avoid exceeding the given batch size. But
            # also ensure that there is at least one vector per batch.
            vectors_per_batch = max(1, mib_to_bytes(batch_size) // self.length)
            # Compute for each batch the positions in the SAM file
            # (given by file.tell and used by file.seek) at which the
            # batch starts and stops.
            indexes = list(reader.get_batch_indexes(vectors_per_batch))
            starts = indexes[:-1]
            stops = indexes[1:]
            n_batches = len(starts)
            batch_nums = list(range(n_batches))
            # Once the number of batches has been determined, a list of
            # new SAM file readers is created. Each is responsible for
            # converting one batch of SAM reads into one batch of
            # mutation vectors. Setting owner=False prevents the SAM
            # reader from creating a new SAM file (as well as from
            # deleting it upon exit); instead, it uses the file written
            # by the primary SAM reader (which is reader.sam_path).
            readers = [SamReader(xam_path=reader.sam_path,
                                 end5=self.end5,
                                 end3=self.end3,
                                 isfullref=self.isfullref,
                                 temp_dir=temp_dir,
                                 n_procs=n_procs,
                                 save_temp=save_temp,
                                 resume=resume,
                                 owner=False)
                       for _ in batch_nums]
            # Create lists of the positional aguments for each reader.
            iter_args = list(zip(batch_nums, readers, starts, stops,
                                 strict=True))
            # Use the same keyword arguments (out_dir and other kwargs)
            # for every batch being vectorized.
            vectorize_batch = partial(self._vectorize_batch,
                                      out_dir=out_dir,
                                      **kwargs)
            if (pool_size := min(n_procs, n_batches)) > 1:
                # Process batches of records simultaneously in parallel.
                with Pool(pool_size) as pool:
                    results = list(pool.starmap(vectorize_batch, iter_args))
            else:
                # Process batches of records one at a time in series.
                results = list(itertools.starmap(vectorize_batch, iter_args))
        # The list of results contains, for each batch, a tuple of the
        # number of mutation vectors in the batch and the MD5 checksum
        # of the batch. Compute the total number of vectors and list all
        # the checksums.
        n_pass = sum(p for p, _, _ in results)
        n_fail = sum(f for _, f, _ in results)
        checksums = [c for _, _, c in results]
        return n_pass, n_fail, n_batches, checksums

    def vectorize(self, /, *, rerun: bool, out_dir: str, **kwargs):
        """ Compute a mutation vector for every record in a BAM file,
        write the vectors into one or more batch files, compute their
        checksums, and write a report summarizing the results. """
        report_path = self.get_report_path(out_dir)
        if self.outputs_valid(out_dir) and not rerun:
            logging.warning(f"Skipping vectorization of {self} because output "
                            f"files already exist. To rerun, add --rerun")
        else:
            # Compute the mutation vectors, write them to batch files,
            # and generate a report.
            began = datetime.now()

            n_p, n_f, n_b, chk = self._vectorize_bam(out_dir=out_dir, **kwargs)
            ended = datetime.now()
            written = self._write_report(out_dir=out_dir,
                                         n_vectors=n_p,
                                         n_readerr=n_f,
                                         n_batches=n_b,
                                         checksums=chk,
                                         began=began,
                                         ended=ended)
            if written != report_path:
                logging.critical(
                    "Intended and actual paths of report differ: "
                    f"{report_path} ≠ {written}")
        return report_path


class WrittenVectors(MutationalProfile):
    """ Represents a collection of mutation vectors that have already
    been written to one or more files. """

    def __init__(self, /, *,
                 n_batches: int,
                 n_vectors: int,
                 checksums: list[str],
                 **kwargs):
        super().__init__(**kwargs)
        self.n_batches = n_batches
        self.n_vectors = n_vectors
        self.checksums = checksums

    @property
    def batch_nums(self) -> list[int]:
        """ Return a list of all batch numbers. """
        return list(range(self.n_batches))

    def get_mv_batch_paths(self, out_dir: str):
        """ Return the path of every mutation vector batch file. """
        return [self.get_mv_batch_path(out_dir, batch)
                for batch in self.batch_nums]


class VectorReport(object):
    """
    Read and write a report about a mutational profile, including:
    - the sample, reference, and section
    - number of mutation vectors
    - number of mutation vector batch files and their checksums
    - beginning and ending time, duration, and speed of vectoring

    Examples
    --------
    >>> report = VectorReport(sample="dms2", ref="tmv-rna", end5=1, end3=20,
    ...                       seq=DNA(b"GTATTTTTACAACAATTACC"),
    ...                       n_vectors=10346, n_batches=2,
    ...                       checksums=["b47260fcad8be60470bee67862f187b4",
    ...                                  "098f40cfadc266ea5bc48ab2e18cdc95"],
    ...                       began=datetime.now(),
    ...                       ended=(time.sleep(1E-5), datetime.now())[-1])
    >>> str(report.seq)
    'GTATTTTTACAACAATTACC'
    """

    class AbstractField(object):
        dtype: type
        key: str
        alias: str

        def __init__(self, value: Any):
            if not isinstance(value, self.dtype):
                raise TypeError(f"{self.__class__.__name__} expected a value "
                                f"of type '{self.dtype.__name__}' but got "
                                f"type '{type(value).__name__}'")
            self.value: Any = value

        @classmethod
        def parse(cls, value: str):
            # For the base class, just cast value to the expected type.
            return cls(cls.dtype(value))

        def __str__(self):
            return str(self.value)

    class AbstractStrField(AbstractField):
        """ String field """
        dtype = str

    class AbstractNumField(AbstractField):
        """ Numeric field with optional minimal and maximal values """
        dtype: int | float
        minimum: int | float | None = None
        maximum: int | float | None = None

        def __init__(self, value: int | float):
            super().__init__(value)
            # Verify bounds. Note that if value is NaN, the < and >
            # comparisons always return False, so the checks will pass,
            # which is the intended behavior. Do not check for
            # if not value >= self.minimum and
            # if not value <= self.maximum
            # because then the checks will fail for NaN values.
            if self.minimum is not None:
                if value < self.minimum:
                    raise ValueError(f"{self.alias} must be ≥ {self.minimum}, "
                                     f"but got {value}")
            if self.maximum is not None:
                if value > self.maximum:
                    raise ValueError(f"{self.alias} must be ≤ {self.maximum}, "
                                     f"but got {value}")

    class AbstractPosIntField(AbstractNumField):
        """ Positive integer """
        dtype, minimum = int, 1

    class AbstractNonNegIntField(AbstractNumField):
        """ Non-negative integer """
        dtype, minimum = int, 0

    class AbstractNonNegFloatField(AbstractNumField):
        """ Non-negative real number """
        dtype, minimum = float, 0.0

    class SampleField(AbstractStrField):
        key, alias = "sample", "Sample Name"

    class RefField(AbstractStrField):
        key, alias = "ref", "Reference Name"

    class End5Field(AbstractPosIntField):
        key, alias = "end5", "Section 5' End"

    class End3Field(AbstractPosIntField):
        key, alias = "end3", "Section 3' End"

    class LengthField(AbstractPosIntField):
        key, alias = "length", "Section Length"

    class AbstractDnaField(AbstractField):
        dtype = DNA

        @classmethod
        def parse(cls, value: str):
            # Need to encode the value to bytes before casting.
            return cls(cls.dtype(value.encode()))

        def __str__(self):
            # Need to decode the value from DNA to str.
            return self.value.decode()

    class SeqField(AbstractDnaField):
        key, alias = "seq", "Section Sequence"

    class BoolField(AbstractField):
        dtype = bool

        @classmethod
        def parse(cls, value: str):
            lower = value.strip().lower()
            if lower == "true":
                return cls(True)
            if lower == "false":
                return cls(False)
            raise ValueError(f"Cannot parse '{value}' as {cls.dtype.__name__}")

    class IsFullField(BoolField):
        key, alias = "isfullref", "Section is Full Reference"

    class NumVectorsField(AbstractNonNegIntField):
        key, alias = "n_vectors", "Reads Vectorized"

    class NumReadErrorsField(AbstractNonNegIntField):
        key, alias = "n_readerr", "Reads with Errors"

    class FracVectorizedField(AbstractNonNegFloatField):
        key, alias = "f_success", "Fraction Vectorized"

    class NumBatchesField(AbstractNonNegIntField):
        key, alias = "n_batches", "Batches"

    class ListStrField(AbstractField):
        dtype, delim = list, ", "

        @classmethod
        def parse(cls, value: str):
            return cls(value.split(cls.delim))

        def __str__(self):
            return self.delim.join(self.value)

    class ChecksumsField(ListStrField):
        key, alias = "checksums", "MD5 Checksums"

    class AbstractDateTimeField(AbstractField):
        dtype, dtfmt = datetime, "%Y-%m-%d %H:%M:%S.%f"

        @classmethod
        def parse(cls, value: str):
            return cls(cls.dtype.strptime(value, cls.dtfmt))

        def __str__(self):
            return self.value.strftime(self.dtfmt)

    class TimeBeganField(AbstractDateTimeField):
        key, alias = "began", "Began"

    class TimeEndedField(AbstractDateTimeField):
        key, alias = "ended", "Ended"

    class TimeTakenField(AbstractNonNegFloatField):
        key, alias = "taken", "Time taken (s)"

    class SpeedField(AbstractNonNegFloatField):
        key, alias = "speed", "Speed (1/s)"
    
    fields = [SampleField, RefField, End5Field, End3Field, LengthField,
              SeqField, IsFullField, NumVectorsField, NumReadErrorsField,
              FracVectorizedField, NumBatchesField, ChecksumsField,
              TimeBeganField, TimeEndedField, TimeTakenField, SpeedField]

    fields_by_key = {field.key: field for field in fields}

    fields_by_alias = {field.alias: field for field in fields}

    @staticmethod
    def compute_length(rep: VectorReport):
        return rep[rep.End3Field] - rep[rep.End5Field] + 1

    @staticmethod
    def compute_f_success(rep: VectorReport):
        try:
            return round((rep[rep.NumVectorsField]
                          / (rep[rep.NumVectorsField]
                             + rep[rep.NumReadErrorsField])
                          ), 5)
        except ZeroDivisionError:
            return float("nan")

    @staticmethod
    def compute_time_taken(rep: VectorReport):
        dtime: timedelta = rep[rep.TimeEndedField] - rep[rep.TimeBeganField]
        return round(dtime.seconds + dtime.microseconds / 1000000, 3)

    @staticmethod
    def compute_speed(rep: VectorReport):
        try:
            return round(rep[rep.NumVectorsField] / rep[rep.TimeTakenField], 1)
        except ZeroDivisionError:
            return float("inf" if rep[rep.NumVectorsField] else "nan")

    compute_fields = {LengthField.key: compute_length,
                      FracVectorizedField.key: compute_f_success,
                      TimeTakenField.key: compute_time_taken,
                      SpeedField.key: compute_speed}

    def __init__(self, /, **kwargs):
        # Fill fields using keyword arguments.
        missing = list()
        for key, field in self.fields_by_key.items():
            if (value := kwargs.get(key)) is None:
                if key not in self.compute_fields:
                    missing.append(key)
            else:
                self[field] = value
        # Raise an error if fields that cannot be computed are missing.
        if missing:
            raise TypeError("Missing the following required keyword arguments: "
                            + ", ".join(missing))
        # Compute any missing fields.
        for key, func in self.compute_fields.items():
            if key not in self:
                self[key] = func(self)

    def __getitem__(self, item: str | type[AbstractField]):
        if isclass(item) and issubclass(item, self.AbstractField):
            item = item.key
        return self.__getattribute__(item)

    def __setitem__(self, field: str | type[AbstractField], value: Any):
        if isinstance(field, str):
            field = self.fields_by_key[field]
        self.__setattr__(field.key, field(value).value)

    def __contains__(self, item: str | type[AbstractField]):
        if isclass(item) and issubclass(item, self.AbstractField):
            item = item.key
        return hasattr(self, item)

    def get_values_by_key(self, fields: Sequence[type[AbstractField]]):
        return {field.key: self[field] for field in fields}

    @cached_property
    def reader_dict(self):
        rdict = self.get_values_by_key([self.SampleField, self.RefField,
                                        self.End5Field, self.End3Field,
                                        self.IsFullField, self.NumVectorsField,
                                        self.SeqField, self.NumBatchesField,
                                        self.ChecksumsField])
        rdict["sect_seq"] = rdict.pop("seq")
        return rdict

    @cached_property
    def path_dict(self):
        return dict(module=path.Module.VECTOR,
                    **self.get_values_by_key([self.SampleField, self.RefField,
                                              self.End5Field, self.End3Field]))

    def find_invalid_batches(self, out_dir: str, validate_checksums: bool):
        """ Return all the batches of mutation vectors that either do not exist
        or do not match their expected checksums. """
        missing = list()
        badsum = list()
        for batch, checksum in enumerate(self[self.ChecksumsField]):
            file = path.MutVectorBatchFilePath(top=out_dir,
                                               **self.path_dict,
                                               batch=batch,
                                               ext=".orc")
            file_path = file.path
            if file_path.is_file():
                if validate_checksums and digest_file(file_path) != checksum:
                    # The batch file exists but does not match the checksum.
                    badsum.append(file)
            else:
                # The batch file does not exist.
                missing.append(file)
        return missing, badsum

    def save(self, out_dir: str):
        dict_strs = {alias: str(field(self[field]))
                     for alias, field in self.fields_by_alias.items()}
        report_path = path.MutVectorReportFilePath(top=out_dir,
                                                   ext=".json",
                                                   **self.path_dict)
        report_path.path.parent.mkdir(parents=True, exist_ok=True)
        with open(report_path.path, "w") as f:
            json.dump(dict_strs, f)
        return report_path

    @classmethod
    def load(cls,
             report_file: str,
             validate_checksums: bool = True) -> VectorReport:
        """ Load a mutation vector report from a file. """
        with open(report_file, "r") as f:
            dict_strs = json.load(f)
        dict_vals = dict()
        for alias, value in dict_strs.items():
            field = cls.fields_by_alias[alias]
            dict_vals[field.key] = field.parse(value).value
        report = cls(**dict_vals)
        report_path = path.MutVectorReportFilePath.parse(report_file)
        if (report_path.sample != report[cls.SampleField]
                or report_path.ref != report[cls.RefField]
                or report_path.end5 != report[cls.End5Field]
                or report_path.end3 != report[cls.End3Field]):
            raise ValueError(f"Report fields do not match path: {report_path}")
        missing, badsum = report.find_invalid_batches(report_path.top,
                                                      validate_checksums)
        if missing:
            raise FileNotFoundError(f"Batches: {', '.join(map(str, missing))}")
        if badsum:
            raise ValueError(f"Bad MD5 sums: {', '.join(map(str, badsum))}")
        return report


class VectorReader(WrittenVectors):
    INDEX_COL = "__index_level_0__"

    def __init__(self, *,
                 out_dir: str,
                 **kwargs):
        super().__init__(**kwargs)
        self.out_dir = out_dir

    @property
    def shape(self):
        return self.n_vectors, self.length

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
        if positions is None:
            columns = None
        else:
            subseq = self.subseq(positions)
            columns = [self.INDEX_COL] + self.seq_pos_to_cols(subseq, positions)
        # Read the vectors from the ORC file using PyArrow as backend.
        vectors = pd.read_orc(self.get_mv_batch_path(self.out_dir, batch).path,
                              columns=columns)
        # Remove the column of read names and set it as the index.
        vectors.set_index(self.INDEX_COL, drop=True, inplace=True)
        # Set the remaining columns to the positions.
        vectors.columns = positions if positions else self.positions
        # The vectors are stored as signed 8-bit integers (np.int8) and
        # must be cast to unsigned 8-bit integers (np.uint8) so that the
        # bitwise operations work. This step must be doneafter removing
        # the column of read names (which cannot be cast to np.uint8).
        return vectors.astype(np.uint8, copy=False)

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
        # If there are no vectors, then return an empty DataFrame.
        if self.n_batches == 0:
            return pd.DataFrame(columns=(self.positions if positions is None
                                         else positions), dtype=int)
        # Load and concatenate every vector batch into one DataFrame.
        return pd.concat(self.get_all_batches(positions), axis=0)

    @classmethod
    def _query_vectors(cls, /,
                       vectors: pd.DataFrame,
                       query: int, *,
                       subsets: bool = False,
                       supersets: bool = False) -> pd.DataFrame:
        """
        Return a boolean array of the same shape as vectors where
        element i,j is True if and only if the byte at element i,j of
        vectors matches the given query byte. By default, a byte in
        vectors matches if it equals the query byte and is not blank
        (i.e. 00000000). Matches can be extended to include bitwise
        subsets and supersets of query by setting the corresponding
        parameters to True. Blank bytes in vectors never match.

        Parameters
        ----------
        vectors: DataFrame
            Mutation vectors
        query: int
            Byte to query; must be in range [0, 255]
        subsets: bool
            Whether to count non-blank bitwise subsets of the query
        supersets: bool
            Whether to count non-blank bitwise supersets of the query

        Return
        ------
        DataFrame
            Boolean type DataFrame of the same shape as vectors where
            each element is True if the element at the same position in
            vectors matched the query and False otherwise
        """
        if not isinstance(query, int):
            raise TypeError(
                f"Expected query of type int, but got {type(query).__name__}")
        if not BLANK <= query <= EVERY:
            raise ValueError(
                f"Expected query in range {BLANK} - {EVERY}, but got {query}")
        if supersets and subsets:
            # Count both supersets and subsets.
            return (cls._query_vectors(vectors, query, supersets=True)
                    | cls._query_vectors(vectors, query, subsets=True))
        if supersets:
            # Non-blank vector bytes that are matches and supersets of
            # the query byte count.
            if query == BLANK:
                # If query is BLANK (00000000), then every non-blank
                # byte is a superset.
                return vectors.astype(bool, copy=True)
            if query == EVERY:
                # If the query is EVERY (11111111), then no bitwise
                # supersets exist. Since subsets do not count, only
                # exact matches count.
                return cls._query_vectors(vectors, query)
            # No shortcut method can be used, so the supersets must be
            # computed explicitly. A vector byte is a match or superset
            # of the query byte iff both of the following are true:
            # - The bitwise intersection of the vector and query bytes
            #   equals the query byte, meaning that every bit set to 1
            #   in the query byte is also set to 1 in the vector byte,
            #   and thus the vector byte is a superset of the query.
            #   Equivalently, the union equals the vector byte.
            # - The vector byte is not blank. But since the query byte
            #   is not blank, if a vector byte satisfies the first
            #   condition, then it must also satisfy this one, so this
            #   latter condition need not be checked.
            return np.equal(np.bitwise_and(vectors, query), query)
        if query == BLANK:
            # If supersets do not count, then only matches and subsets
            # count. But if the query is BLANK (00000000), then there
            # are no subsets, and matches do not count because blank
            # vector bytes never count. Thus, no byte counts.
            return pd.DataFrame(False, dtype=bool,
                                index=vectors.index,
                                columns=vectors.columns)
        if subsets:
            # Non-blank vector bytes that are matches and subsets of the
            # query byte count.
            if query == EVERY:
                # If query is EVERY (11111111), then every non-blank
                # byte is a subset.
                return vectors.astype(bool, copy=True)
            if (query & (query - 1)) == 0:
                # If query is a power of 2, then it has exactly one bit
                # set to 1. Thus, the only possible subset of the query
                # is the blank byte, which never counts. Since supersets
                # do not count either, only exact matches count.
                return cls._query_vectors(vectors, query)
            # No shortcut method can be used, so the subsets must be
            # computed explicitly. A vector byte is a match or subset of
            # the query byte iff both of the following are true:
            # - The bitwise union of the vector and query bytes equals
            #   the query byte, meaning that there are no bits set to 1
            #   in the vector byte that are not 1 in the query byte,
            #   and thus the vector byte is a subset of the query.
            #   Equivalently, the intersection equals the vector byte.
            # - The vector byte is not blank.
            return (vectors.astype(bool, copy=False)
                    & np.equal(np.bitwise_or(vectors, query), query))
        # If neither subsets nor supersets count and query is not BLANK,
        # then count vector bytes that match the query byte exactly.
        return np.equal(vectors, query)

    @cache
    def count_muts_by_pos(self, /,
                          query: int, *,
                          subsets: bool = False,
                          supersets: bool = False,
                          positions: Sequence[int] | None = None) -> pd.Series:
        """
        Return the number of mutations that match the query at each
        position in the mutational profile.

        Parameters
        ----------
        query: int
            Byte to query; must be in range [0, 255]
        subsets: bool
            Whether to count non-blank bitwise subsets of the query
        supersets: bool
            Whether to count non-blank bitwise supersets of the query
        positions: sequence[int] | None
            Use only these positions from the mutation vectors; if None,
            then use all positions.

        Return
        ------
        Series
            Number of mutations matching the query at each position
        """
        # Initialize empty Series to count mutations at each position.
        counts = pd.Series(np.zeros(self.length, dtype=int),
                           index=(self.positions if positions is None
                                  else positions))
        # Iterate over all batches of vectors.
        for vectors in self.get_all_batches(positions):
            # Add the number of mutations at each position in the batch
            # to the cumulative count of mutations at each position.
            counts += self._query_vectors(vectors,
                                          query,
                                          subsets=subsets,
                                          supersets=supersets).sum(axis=0)
        return counts

    @cache
    def count_muts_by_vec(self, /,
                          query: int, *,
                          subsets: bool = False,
                          supersets: bool = False,
                          positions: Sequence[int] | None = None) -> pd.Series:
        """
        Return the number of mutations that match the query for each
        vector in the mutational profile.

        Parameters
        ----------
        query: int
            Byte to query; must be in range [0, 255]
        subsets: bool
            Whether to count non-blank bitwise subsets of the query
        supersets: bool
            Whether to count non-blank bitwise supersets of the query
        positions: sequence[int] | None
            Use only these positions from the mutation vectors; if None,
            then use all positions.

        Return
        ------
        Series
            Number of mutations matching the query in each vector
        """
        # If there are no vectors, then return an empty Series.
        if self.n_batches == 0:
            return pd.Series([], dtype=int)
        # Initialize empty list to count the mutations in each vector.
        counts = list()
        # Iterate over all batches of vectors.
        for vectors in self.get_all_batches(positions):
            # Count the number of mutations in each vector in the batch
            # and append them to the list of counts.
            counts.append(self._query_vectors(vectors,
                                              query,
                                              subsets=subsets,
                                              supersets=supersets).sum(axis=1))
        # Concatenate and return the number of mutations in each vector
        # among all batches.
        return pd.concat(counts, axis=0)


    @classmethod
    def load(cls, report_file: str, validate_checksums: bool = True):
        if not (report := VectorReport.load(report_file, validate_checksums)):
            logging.critical(f"Failed to load report from {report_file}")
            return
        return cls(out_dir=path.MutVectorReportFilePath.parse(report_file).top,
                   **report.reader_dict)


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


def get_sections(ref_seqs: dict[str, DNA], *,
                coords: tuple[tuple[str, int, int]],
                primers: tuple[tuple[str, DNA, DNA]],
                primer_gap: int,
                cfill: bool):
    """ Return all the sections corresponding to the given coordinates
    and/or primers in the given reference sequences. """
    sections: dict[str, list[SectionFinder]] = defaultdict(list)

    def add_section(section: SectionFinder):
        if any(section == other for other in sections[section.ref]):
            logging.warning(f"Skipping duplicate section: {section.ref_coords}")
        sections[section.ref].append(section)

    for ref, first, last in coords:
        add_section(SectionFinder(ref_seq=ref_seqs.get(ref), ref=ref,
                                  end5=first, end3=last, primer_gap=primer_gap))
    for ref, fwd, rev in primers:
        add_section(SectionFinder(ref_seq=ref_seqs.get(ref), ref=ref,
                                  fwd=fwd, rev=rev, primer_gap=primer_gap))
    if cfill:
        for ref, seq in ref_seqs.items():
            if ref not in sections:
                add_section(SectionFinder(ref_seq=seq, ref=ref,
                                          primer_gap=primer_gap))
    return sections


def get_writers(fasta: str,
                bam_files: list[str],
                **kwargs):
    ref_seqs = dict(FastaParser(fasta).parse())
    sections = get_sections(ref_seqs, **kwargs)
    writers: dict[tuple, VectorWriter] = dict()
    for bam_file in bam_files:
        bam_path = path.OneRefAlignmentInFilePath.parse(bam_file)
        for section in sections[bam_path.ref]:
            if section.ref != bam_path.ref:
                logging.error(f"Skipping section {section} of {bam_path.path} "
                              "because its reference does not match that "
                              f"of the BAM file ('{bam_path.ref}').")
                continue
            writer = VectorWriter(bam_file=bam_path,
                                  ref_seq=ref_seqs[bam_path.ref],
                                  end5=section.end5,
                                  end3=section.end3)
            if writer.tag in writers:
                logging.warning("Skipping duplicate mutational profile: "
                                f"{writer}.")
                continue
            writers[writer.tag] = writer
    return list(writers.values())


def generate_profiles(writers: list[VectorWriter], *,
                      phred_enc: int,
                      min_phred: int,
                      max_procs: int,
                      parallel: bool,
                      **kwargs) -> tuple[str, ...]:
    """ Generate mutational profiles of one or more vector writers. """
    n_profiles = len(writers)
    if n_profiles == 0:
        logging.warning("No BAM files and/or sections specified")
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
            report_files = tuple(pool.starmap(vectorize, iter_args))
    else:
        report_files = tuple(itertools.starmap(vectorize, iter_args))
    return tuple(map(str, report_files))


vector_trans_table = bytes.maketrans(*map(b"".join, zip(*[(
    i.to_bytes(length=1, byteorder=sys.byteorder),
    (
        b"." if i == BLANK
        else b"~" if i == MATCH
        else b"/" if i == DELET
        else b"{" if i == (INS_5 | MATCH)
        else b"}" if i == (INS_3 | MATCH)
        else b"A" if i == SUB_A
        else b"C" if i == SUB_C
        else b"G" if i == SUB_G
        else b"T" if i == SUB_T
        else b"?"
    )
) for i in range(256)])))


def trans_vectors_iter(vectors: pd.DataFrame):
    for index, row in zip(vectors.index, vectors.values, strict=True):
        translated = row.tobytes(order='C').translate(vector_trans_table)
        yield b"%b\t%b\n" % (index.encode(), translated)


def trans_vectors_block(vectors: pd.DataFrame):
    return b"".join(trans_vectors_iter(vectors))
