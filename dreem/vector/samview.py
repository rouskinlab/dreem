from __future__ import annotations
from functools import cached_property, wraps
from io import BufferedReader
from typing import Optional

from dreem.util.reads import XamBase, BamVectorSelector, SamVectorSorter
from dreem.util.path import BasePath, XamInPath, XamTempPath
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
    def __init__(self, base_path: BasePath, xam_path: XamInPath,
                 ref_name: str, first: int, last: int, spanning: bool,
                 owner: bool = True):
        self.base_path = base_path
        self.xam_path = xam_path
        self.ref_name = ref_name
        self.first = first
        self.last = last
        self.spanning = spanning
        self.owner = owner
        self._sam_path: Optional[XamInPath | XamTempPath] = None
        self._sam_file: Optional[BufferedReader] = None
    
    def __enter__(self):
        # Convert the BAM file to a temporary SAM file
        if self.owner:
            if self.spanning:
                selector = None
                xam_path = self.xam_path
            else:
                xam_base = XamBase(self.base_path, self.xam_path)
                xam_index = xam_base.xam_index
                if not xam_index.path.is_file():
                    assert xam_base.create_index().path == xam_index.path
                selector = BamVectorSelector(self.base_path, self.xam_path)
                xam_path = selector.run(self.ref_name, self.first, self.last)
            sorter = SamVectorSorter(self.base_path, xam_path)
            self._sam_path = sorter.run(name=True)
            if selector:
                selector.clean()
        else:
            self._sam_path = self.xam_path
        self._sam_file = open(self.sam_path.path, "rb")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._sam_file.close()
        self._sam_file = None
        if self.owner:
            self.sam_path.path.unlink()
            self._sam_path = None
    
    @property
    def sam_path(self):
        if self._sam_path:
            return self._sam_path
        raise RuntimeError(f"{self} has no open SAM path.")
    
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
        self._seek_rec1()
        first_line = self._sam_file.readline()
        return SamRead(first_line).flag.paired if first_line else False
    
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
            if self.spanning:
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
