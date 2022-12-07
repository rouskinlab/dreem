from io import FileIO
import math
from typing import Optional

from util import *
from vector.vector import *


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
