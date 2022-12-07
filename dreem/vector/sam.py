from io import FileIO
import itertools
import math
import re
from typing import List, Optional, Set, Tuple

from core import *


CIG_MAT = b"="
CIG_SUB = b"X"
CIG_DEL = b"D"
CIG_INS = b"I"
CIG_SCL = b"S"
CIG_PATTERN = re.compile(b"".join(
    [rb"(\d+)([", CIG_MAT, CIG_SUB, CIG_DEL, CIG_INS, CIG_SCL, b"])"]
))
SAM_HEADER = b"@"
NOLIM = -1
INT_A, INT_C, INT_G, INT_T, INT_N = BASES + BASEN


def parse_cigar(cigar_string: bytes):
    if not cigar_string:
        raise ValueError("CIGAR string was empty.")
    length_matched = 0
    for match in CIG_PATTERN.finditer(cigar_string):
        length_bytes, operation = match.groups()
        length_matched += len(length_bytes) + len(operation)
        if (length_int := int(length_bytes)) < 1:
            raise ValueError("length of CIGAR operation must be >= 1")
        yield operation, length_int
    if length_matched != len(cigar_string):
        raise ValueError(f"Invalid CIGAR string: '{cigar_string.decode()}'")


class SamFlag(object):
    __slots__ = ["paired", "proper", "unmap", "munmap", "rev", "mrev",
                 "first", "second", "secondary", "qcfail", "dup", "supp"]

    max_flag = 2**len(__slots__) - 1
    pattern = "".join(["{:0<", str(len(__slots__)), "}"])

    def __init__(self, flag: int):
        if not 0 <= flag <= self.max_flag:
            raise ValueError(f"Invalid flag: '{flag}'")
        (self.paired, self.proper, self.unmap, self.munmap, self.rev,
         self.mrev, self.first, self.second, self.secondary, self.qcfail,
         self.dup, self.supp) = (x == "1" for x in
                                 self.pattern.format(bin(flag)[:1:-1]))


class SamRead(object):
    __slots__ = ["qname", "flag", "rname", "pos", "mapq", "cigar",
                 "tlen", "seq", "qual"]
    
    min_fields = 11

    def __init__(self, line: bytes):
        fields = line.rstrip().split(b"\t")
        if len(fields) < self.min_fields:
            raise ValueError(f"Invalid SAM line:\n{line}")
        self.qname = fields[0]
        self.flag = SamFlag(int(fields[1]))
        self.rname = fields[2]
        self.pos = int(fields[3])
        self.mapq = int(fields[4])
        self.cigar = fields[5]
        #self.rnext = fields[6]
        #self.pnext = int(fields[7])
        self.tlen = int(fields[8])
        self.seq = AmbigDNA(fields[9])
        self.qual = fields[10]
        if len(self) != len(self.qual):
            raise ValueError(f"Lengths of seq ({len(self)}) and qual "
                             f"string {len(self.qual)} did not match.")
        
    def __len__(self):
        return len(self.seq)


def get_consensus_mut(byte1: int, byte2: int):
    return intersect if (intersect := byte1 & byte2) else byte1 | byte2


def get_muts_pair(region_seq: bytes, region_first: int, region_last: int,
                  read1: SamRead, read2: SamRead):
    muts1 = get_muts_read(region_seq, region_first, region_last, read1)
    muts2 = get_muts_read(region_seq, region_first, region_last, read2)
    return bytearray(map(get_consensus_mut, muts1, muts2))


def op_consumes_ref(op: bytes):
    return op == CIG_MAT or op == CIG_SUB or op == CIG_DEL


def op_consumes_read(op: bytes):
    return op == CIG_MAT or op == CIG_SUB or op == CIG_INS or op == CIG_SCL


def encode_base(base: int):
    if base == INT_T:
        return SUB_T
    if base == INT_G:
        return SUB_G
    if base == INT_C:
        return SUB_C
    if base == INT_A:
        return SUB_A
    raise ValueError(f"Invalid base: {base.to_bytes().decode()}")


def encode_compare(ref_base: int, read_base: int):
    return (MATCH[0] if ref_base == read_base
            else (ANY_N ^ encode_base(ref_base) if read_base == INT_N
                  else encode_base(read_base)))


def get_muts_read(region_seq: bytes, region_first: int, region_last: int,
                  read: SamRead):
    """
    :param region_seq: str, reference sequence (must contain T, not U)
    :param region_first: int, the first coordinate (w.r.t. the reference
        sequence) of the region of interest (1-indexed)
    :param region_last: int, the last coordinate (w.r.t. the reference
        sequence) of the region of interest, inclusive (1-indexed)
    :param read: SamRead, the read for which to compute mutations
    :return:
    """
    region_length = region_last - region_first + 1
    assert region_length == len(region_seq)
    # current position in the read
    # 0-indexed from beginning of read
    curr_read_pos = 0
    # position at which the current CIGAR operation starts
    # 0-indexed from beginning of region
    curr_op_start = read.pos - region_first
    # position at which the current CIGAR operation ends (0-indexed)
    # does not include the last position of the operation (like Python slicing)
    # 0-indexed from beginning of region
    curr_op_end = curr_op_start
    # Initialize the mutation vector. Pad the beginning with missing bytes
    # if the read starts after the first position in the region.
    muts = bytearray(BLANK * min(curr_op_start, region_length))
    # Keep track of all positions with deletions and insertions.
    # 0-indexed from beginning of region
    del_ref_coords: List[int] = list()  # deletions in reference coordinates
    del_read_coords: List[int] = list()  # deletions in read coordinates
    ins_ref_coords: List[int] = list()  # insertions in reference coordinates
    ins_read_coords: List[int] = list()  # insertions in read coordinates
    # Read the CIGAR string one operation at a time.
    for cigar_op, op_length in parse_cigar(read.cigar):
        if op_consumes_ref(cigar_op):
            # Advance the end of the operation if it consumes the reference.
            curr_op_end += op_length
        if curr_op_end > 0 and curr_op_start < region_length:
            # Run this block once the operation has entered the region.
            if curr_op_start < 0:
                # If the current operation starts before the region,
                # then advance it to the beginning of the region.
                op_length += curr_op_start
                curr_op_start = 0
            if curr_op_end > region_length:
                # If the current operation ends after the region,
                # then truncate it to the end of the region.
                op_length -= curr_op_end - region_length
                curr_op_end = region_length
            # Perform an action based on the CIGAR operation and its length.
            if cigar_op == CIG_MAT:
                # Condition: read matches the reference at this position
                # The quality must meet the minimum quality to count as a match
                # Low quality bases were already masked as Ns before alignment
                muts.extend(MATCH * op_length)
            elif cigar_op == CIG_SUB:
                # Condition: read has a substitution w.r.t. reference at this
                # position or has a low-quality base that was masked as an N.
                muts.extend(map(encode_compare,
                                region_seq[len(muts): len(muts) + op_length],
                                read.seq[curr_read_pos:
                                         curr_read_pos + op_length]))
            elif cigar_op == CIG_DEL:
                # Condition: read contains a deletion w.r.t. the reference.
                del_ref_coords.extend(range(len(muts), len(muts) + op_length))
                del_read_coords.extend((curr_read_pos,) * op_length)
                muts.extend(DELET * op_length)
            elif cigar_op == CIG_INS:
                # Condition: read contains an insertion w.r.t. the reference.
                # Position added to insertions is of the base 3' of the insert.
                ins_ref_coords.append((len(muts),) * op_length)
                ins_read_coords.append(range(curr_read_pos,
                                             curr_read_pos + op_length))
                # Insertions do not consume the reference, so do not add any
                # information to muts yet. That information is added later.
            elif cigar_op == CIG_SCL:
                # Condition: read contains a soft clipping
                pass
            else:
                raise ValueError(
                    f"Invalid CIGAR operation: '{cigar_op.decode()}'")
        if op_consumes_read(cigar_op):
            # Advance the read position if the operation consumes the read.
            curr_read_pos += op_length
        # Advance the start position to the end of the current operation.
        curr_op_start = curr_op_end
    # Pad the end of the mutation vector with any non-covered positions.
    muts.extend(BLANK * (region_length - len(muts)))
    assert len(muts) == region_length
    # Ensure the CIGAR string matched the length of the read.
    if curr_read_pos != len(read):
        raise ValueError(
            f"CIGAR string '{read.cigar.decode()}' consumed {curr_read_pos} "
            f"bases from read, but read is {len(read)} bases long.")
    # Flag all positions that are ambiguous due to indels.
    #flag_ambigs(muts, region_seq, read.seq,
    #            del_ref_coords, del_read_coords,
    #            ins_ref_coords, ins_read_coords)
    return muts


class SamRecord(object):
    __slots__ = ["read1", "read2"]

    def __init__(self, read1: SamRead, read2: Optional[SamRead] = None):
        self.read1 = read1
        if read2 is not None:
            if self.paired:
                if not read2.flag.paired:
                    raise ValueError("read1 is paired but read2 is not")
                if read1.qname != read2.qname:
                    raise ValueError(
                        "Paired reads had inconsistent query names: "
                        f"'{read1.qname}' and '{read2.qname}'")
                if read1.rname != read2.rname:
                    raise ValueError(
                        "Paired reads had inconsistent reference names: "
                        f"'{read1.rname}' and '{read2.rname}'")
                if abs(read1.tlen) != abs(read2.tlen):
                    raise ValueError(
                        "Paired reads had inconsistent template lengths: "
                        f"{read1.tlen} and {read2.tlen}")
                if read1.flag.second or not read1.flag.first:
                    raise ValueError("read1 is not flagged as first read")
                if read2.flag.first or not read2.flag.second:
                    raise ValueError("read2 is not flagged as second read")
                if read1.flag.rev == read2.flag.rev:
                    raise ValueError(
                        "read1 and read2 are in the same orientation")
            else:
                raise ValueError("read1 is unpaired, but read2 was given")
        self.read2 = read2

    @property
    def ref_name(self):
        return self.read1.rname

    @property
    def paired(self):
        return self.read1.flag.paired

    def comp_muts(self, region_seq, region_start, region_end):
        if self.read2 is None:
            return get_muts_read(region_seq, region_start, region_end,
                                 self.read1)
        else:
            return get_muts_pair(region_seq, region_start, region_end,
                                 self.read1, self.read2)


class SamViewer(object):
    __slots__ = ["bam_path", "ref_name", "start", "end", "spanning",
                 "_sam_path", "_sam_file", "_paired", "_n_lines", "_n_records",
                 "_n_records_est"]

    def __init__(self, bam_path: str, ref_name: bytes, start: int, end: int,
                 spanning: bool):
        self.bam_path = bam_path
        self.ref_name = ref_name
        self.start = start
        self.end = end
        self.spanning = spanning
        self._sam_path: Optional[str] = None
        self._sam_file: Optional[FileIO] = None
        self._paired: Optional[bool] = None
        self._n_lines: Optional[int] = None
        self._n_records: Optional[int] = None
        self._n_records_est: Optional[int] = None

    def __enter__(self):
        # Convert the BAM file to a temporary SAM file
        bam_dir, bam_base = os.path.split(self.bam_path)
        bam_name, ext = os.path.splitext(bam_base)
        self._sam_path = name_temp_file(dirname=bam_dir, prefix=f"{bam_name}_",
                                        suffix=".sam")
        cmd = (f"samtools view -h '{self.bam_path}' "
               f"{self.ref_name.decode()}:{self.start}-{self.end} "
               f"| samtools sort -n -o '{self._sam_path}'")
        os.system(cmd)
        self._sam_file = open(self._sam_path, "rb")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._sam_file.close()
        try_remove(self._sam_path)
    
    @staticmethod
    def _reset_seek(func):
        def wrapper(self, *args, **kwargs):
            prev_pos = self._sam_file.tell()
            func(self, *args, **kwargs)
            self._sam_file.seek(prev_pos)
        return wrapper

    def _seek_beginning(self):
        self._sam_file.seek(0)
    
    def _seek_record_1(self):
        while (line := self._sam_file.readline()).startswith(SAM_HEADER):
            pass
        self._sam_file.seek(-len(line), 1)
    
    def _seek_record_by_start(self, start: Optional[int] = None):
        if start is None:
            self._seek_record_1()
        else:
            self._sam_file.seek(start)
    
    def _assert_record_by_stop(self, stop: Optional[int] = None):
        if stop is not None and stop != self._sam_file.tell():
            raise ValueError(f"Requested stopping at {stop} "
                f"but stopped at {self._sam_file.tell()} instead.")
    
    @staticmethod
    def _range_of_records(get_records):
        def wrapper(self, start, stop):
            self._seek_record_by_start(start)
            for record in get_records(self):
                yield record
            self._assert_record_by_stop(stop)
        return wrapper

    @property
    def paired(self):
        if self._paired is None:
            prev_pos = self._seek_record_1()
            first_line = self._sam_file.readline()
            if first_line:
                self._paired = SamRead(first_line).flag.paired
            else:
                self._paired = False
            self._sam_file.seek(prev_pos)
        return self._paired

    @property
    def n_records(self):
        """
        Count the number of records in the SAM file.
        :return:
        """
        if self._n_records is None:
            self._n_records = sum(1 for _ in self.get_records())
        return self._n_records

    @property
    def est_records(self):
        """
        Estimate the number of records in the SAM file. If spanning is True
        or all records are single-end reads, then this function returns the
        exact number of records (the same as n_records). Otherwise, it may
        underestimate the number of records because paired-end reads where
        only one mate overlaps the region of interest will appear only once
        (while paired-end reads are expected to appear twice). However, this
        function is much faster than n_records because the latter needs to
        construct a SamRecord object for every record in the file, while
        est_records merely counts the number of lines.
        :return:
        """
        if self._n_records_est is None:
            if self._n_records is None:
                prev_pos = self._seek_record_1()
                n = sum(1 for _ in self._sam_file)
                self._n_records_est = math.ceil(n / 2) if self.paired else n
                self._sam_file.seek(prev_pos)
            else:
                self._n_records_est = self._n_records
        return self._n_records_est
    
    @_reset_seek
    @_range_of_records
    def _get_records_single(self):
        while read := SamRead(self._sam_file.readline()):
            if read.flag.paired:
                raise ValueError("Found paired-end read in single-end file")
            yield SamRecord(read)

    @_reset_seek
    @_range_of_records
    def _get_records_paired_flexible(self):
        prev_read: Optional[SamRead] = None
        while line := self._sam_file.readline():
            if not (read := SamRead(line)).flag.paired:
                raise ValueError("Found single-end read in paired-end file")
            if prev_read:
                # The previous read has not yet been yielded
                if prev_read.qname == read.qname:
                    # The current read is the mate of the previous read
                    if prev_read.flag.first:
                        yield SamRecord(prev_read, read)
                    else:
                        yield SamRecord(read, prev_read)
                    prev_read = None
                else:
                    # The previous read is paired, but its mate is not in
                    # the SAM file. This situation can happen when only the
                    # reads mapping to a specific region are exported to a
                    # SAM file and only one mate overlaps that region.
                    yield SamRecord(prev_read)
                    # Save the current read so that if its mate is the next
                    # read, it will be returned as a pair.
                    prev_read = read
            else:
                # Save the current read so that if its mate is the next
                # read, it will be returned as a pair.
                prev_read = read
        if prev_read:
            # In case the last record
            yield SamRecord(prev_read)
    
    @_reset_seek
    @_range_of_records
    def _get_records_paired_strict(self):
        while line := self._sam_file.readline():
            if not (read := SamRead(line)).flag.paired:
                raise ValueError("Found single-end read in paired-end file")
            yield SamRecord(read, SamRead(self._sam_file.readline()))
    
    def get_records(self, start: Optional[int] = None,
                    stop: Optional[int] = None):
        if self.paired:
            if self.spanning:
                records = self._get_records_paired_strict
            else:
                records = self._get_records_paired_flexible
        else:
            records = self._get_records_single
        return records(start, stop)
    
    def get_batch_indexes(self, batch_size: int):
        if batch_size <= 0:
            raise ValueError("batch_size must be a positive integer")
        records = self.get_records()
        while True:
            index = self._sam_file.tell()
            try:
                next(records)
            except StopIteration:
                return
            else:
                yield index
                for _, _ in zip(range(batch_size - 1), records):
                    pass
