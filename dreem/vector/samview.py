from __future__ import annotations
from functools import cached_property, wraps
from io import BufferedReader
import os
from typing import Optional

from dreem.util.util import name_temp_file, SAMTOOLS_CMD, run_cmd, try_remove
from dreem.vector.vector import *


def _requires_open(func: function):
    @wraps(func)
    def wrapper(self: SamViewer, *args, **kwargs):
        if self._sam_file is None:
            raise ValueError(f"Function '{func.__name__}' requires "
                             "the SAM file to have been opened.")
        return func(self, *args, **kwargs)
    return wrapper


def _reset_seek(func: function):
    @wraps(func)
    @_requires_open
    def wrapper(self: SamViewer, *args, **kwargs):
        prev_pos = self._sam_file.tell()
        ret = func(self, *args, **kwargs)
        self._sam_file.seek(prev_pos)
        return ret
    return wrapper


def _range_of_records(func: function):
    @wraps(func)
    @_reset_seek
    def wrapper(self: SamViewer, start: int, stop: int):
        self._sam_file.seek(start)
        records = iter(func(self))
        if stop is None:
            return records
        else:
            while True:
                if self._sam_file.tell() < stop:
                    yield next(records)
                else:
                    assert self._sam_file.tell() == stop
                    break
    return wrapper


class SamViewer(object):
    def __init__(self, xam_path: str, ref_name: bytes, first: int, last: int,
                 spanning: bool, make: bool = True, remove: bool = True):
        self._xam_in = xam_path
        self._ref_name = ref_name
        self._first = first
        self._last = last
        self._spanning = spanning
        self._make = make
        self._remove = remove
        self._sam_temp: str = ""
        self._sam_file: Optional[BufferedReader] = None
        self._paired: Optional[bool] = None
    
    @property
    def input_dir(self):
        return os.path.dirname(self._xam_in)
    
    @property
    def input_name(self):
        return os.path.splitext(os.path.basename(self._xam_in))[0]
    
    @property
    def ref_coords(self):
        return f"{self._ref_name.decode()}:{self._first}-{self._last}"
    
    @property
    def sam_subset(self):
        if not self._sam_temp:
            raise ValueError("Not currently viewing SAM subset file")
        return self._sam_temp
    
    def __enter__(self):
        # Convert the BAM file to a temporary SAM file
        
        # FIXME !!!
        # Use functions defined in fastq.py to do these conversions elegantly
        # Write the SAM files to the temporary directory instead of the same
        # directory with the BAM file (in case that directory is read-only)

        if self._make:
            cmd = [SAMTOOLS_CMD, "index", self._xam_in]
            run_cmd(cmd)

            sam_view = name_temp_file(dirname=self.input_dir,
                                    prefix=f"{self.input_name}_nosort_",
                                    suffix=".sam")
            cmd = [SAMTOOLS_CMD, "view", "-h", "-o", sam_view,
                self._xam_in, self.ref_coords]
            run_cmd(cmd)

            self._sam_temp = name_temp_file(dirname=self.input_dir,
                                            prefix=f"{self.input_name}_",
                                            suffix=".sam")
            cmd = [SAMTOOLS_CMD, "sort", "-n", "-o", self._sam_temp, sam_view]
            run_cmd(cmd)
            try_remove(sam_view)
            self._sam_file = open(self._sam_temp, "rb")
        else:
            self._sam_file = open(self._xam_in, "rb")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._sam_file.close()
        self._sam_file = None
        if self._remove:
            try_remove(self._sam_temp)
        self._sam_temp = ""

    @_requires_open
    def _seek_beginning(self):
        self._sam_file.seek(0)
    
    @cached_property
    @_reset_seek
    def _rec1_pos(self):
        self._seek_beginning()
        while (line := self._sam_file.readline()).startswith(SAM_HEADER):
            pass
        return self._sam_file.tell() - len(line)
    
    @_requires_open
    def _seek_rec1(self):
        self._sam_file.seek(self._rec1_pos)
        
    @cached_property
    @_reset_seek
    def paired(self):
        if self._paired is None:
            self._seek_rec1()
            first_line = self._sam_file.readline()
            self._paired = (SamRead(first_line).flag.paired if first_line
                            else False)
        return self._paired
    
    @_range_of_records
    def _get_records_single(self):
        while read := SamRead(self._sam_file.readline()):
            assert not read.flag.paired
            yield SamRecord(read)

    @_range_of_records
    def _get_records_paired_flexible(self):
        prev_read: Optional[SamRead] = None
        while line := self._sam_file.readline():
            read = SamRead(line)
            assert read.flag.paired
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
    
    @_range_of_records
    def _get_records_paired_strict(self):
        while line := self._sam_file.readline():
            yield SamRecord(SamRead(line), SamRead(self._sam_file.readline()))
    
    def get_records(self, start: int, stop: int):
        if self.paired:
            if self._spanning:
                records = self._get_records_paired_strict
            else:
                records = self._get_records_paired_flexible
        else:
            records = self._get_records_single
        return records(start, stop)
    
    @_reset_seek
    def get_batch_indexes(self, batch_size: int):
        if batch_size <= 0:
            raise ValueError("batch_size must be a positive integer")
        n_skip = (self.paired + 1) * (batch_size - 1)
        self._seek_rec1()
        while True:
            yield self._sam_file.tell()
            try:
                line = next(self._sam_file)
            except StopIteration:
                break
            else:
                for _, line in zip(range(n_skip), self._sam_file):
                    pass
                if self.paired:
                    # Compare the current and next query names.
                    try:
                        qname = line.split()[0]
                    except IndexError:
                        pass
                    else:
                        line_next = self._sam_file.readline()
                        try:
                            qname_next = line_next.split()[0]
                        except IndexError:
                            pass
                        else:
                            if qname != qname_next:
                                # If the current and next query names differ
                                # (the lines do not come from two paired mates),
                                # then backtrack to the beginning of line_next.
                                self._sam_file.seek(-len(line_next), 1)
