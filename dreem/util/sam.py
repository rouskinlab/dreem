from __future__ import annotations
from io import FileIO
import math
import os
from tqdm import tqdm
from typing import Optional

from dreem.util.util import name_temp_file, SAMTOOLS_CMD, run_cmd, try_remove
from dreem.vector.vector import *


class SamViewer(object):
    __slots__ = ["_bam_path", "_ref_name", "_first", "_last", "_spanning",
                 "_sam_path", "_sam_file", "_paired", "_n_lines", "_n_records",
                 "_n_records_est"]

    def __init__(self, bam_path: str, ref_name: bytes, first: int, last: int,
                 spanning: bool):
        self._bam_path = bam_path
        self._ref_name = ref_name
        self._first = first
        self._last = last
        self._spanning = spanning
        self._sam_path: Optional[str] = None
        self._sam_file: Optional[FileIO] = None
        self._paired: Optional[bool] = None
    
    @property
    def bam_path(self):
        return self._bam_path
    
    @property
    def bam_dir(self):
        return os.path.dirname(self.bam_path)
    
    @property
    def bam_name(self):
        return os.path.splitext(os.path.basename(self.bam_path))[0]
    
    @property
    def ref_name(self):
        return self._ref_name
    
    @property
    def first(self):
        return self._first
    
    @property
    def last(self):
        return self._last
    
    @property
    def spanning(self):
        return self._spanning
    
    @property
    def sam_path(self):
        return self._sam_path
    
    @property
    def sam_file(self):
        return self._sam_file
    
    @property
    def ref_coords(self):
        return f"{self.ref_name.decode()}:{self.first}-{self.last}"
    
    def __enter__(self):
        # Convert the BAM file to a temporary SAM file
        
        # FIXME !!!
        # Use functions defined in fastq.py to do these conversions elegantly
        # Write the SAM files to the temporary directory instead of the same
        # directory with the BAM file (in case that directory is read-only)

        cmd = [SAMTOOLS_CMD, "index", self.bam_path]
        run_cmd(cmd)

        sam_view = name_temp_file(dirname=self.bam_dir,
                                  prefix=f"{self.bam_name}_nosort_",
                                  suffix=".sam")
        cmd = [SAMTOOLS_CMD, "view", "-h", "-o", sam_view,
               self.bam_path, self.ref_coords]
        run_cmd(cmd)

        self._sam_path = name_temp_file(dirname=self.bam_dir,
                                        prefix=f"{self.bam_name}_",
                                        suffix=".sam")
        cmd = [SAMTOOLS_CMD, "sort", "-n", "-o", self.sam_path, sam_view]
        run_cmd(cmd)
        try_remove(sam_view)
        self._sam_file = open(self.sam_path, "rb")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.sam_file.close()
        try_remove(self.sam_path)
    
    @staticmethod
    def _reset_seek(func):
        def wrapper(self: SamViewer, *args, **kwargs):
            prev_pos = self.sam_file.tell()
            ret = func(self, *args, **kwargs)
            self.sam_file.seek(prev_pos)
            return ret
        return wrapper

    def _seek_beginning(self):
        self.sam_file.seek(0)
    
    def _seek_record_1(self):
        prev_pos = self.sam_file.tell()
        while (line := self.sam_file.readline()).startswith(SAM_HEADER):
            pass
        self.sam_file.seek(-len(line), 1)
        return prev_pos
    
    def _seek_record_by_start(self, start: Optional[int] = None):
        if start is None:
            self._seek_record_1()
        else:
            self.sam_file.seek(start)
    
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
