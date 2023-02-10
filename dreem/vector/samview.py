from __future__ import annotations
from functools import cached_property, wraps
from io import BufferedReader
from typing import Callable, Optional

from dreem.util.reads import BamVectorSelector, SamVectorSorter
from dreem.util.path import TopDirPath, OneRefAlignmentInFilePath, OneRefAlignmentTempFilePath
from dreem.vector.vector import *


def _reset_seek(func: Callable):
    @wraps(func)
    def wrapper(samview: SamViewer, *args, **kwargs):
        prev_pos = samview.sam_file.tell()
        result = func(samview, *args, **kwargs)
        samview.sam_file.seek(prev_pos)
        return result
    return wrapper


def _range_of_records(func: Callable):
    @wraps(func)
    @_reset_seek
    def wrapper(samview: SamViewer, start: int, stop: int):
        samview.sam_file.seek(start)
        records = iter(func(samview))
        if stop is None:
            return records
        else:
            while True:
                if samview.sam_file.tell() < stop:
                    yield next(records)
                else:
                    assert samview.sam_file.tell() == stop
                    break
    return wrapper


class SamViewer(object):
    def __init__(self,
                 top_dir: TopDirPath,
                 max_cpus: int,
                 xam_path: OneRefAlignmentInFilePath,
                 ref_name: str,
                 end5: int,
                 end3: int,
                 spanning: bool,
                 min_qual: int,
                 owner: bool = True):
        self.top_dir = top_dir
        self.max_cpus = max_cpus
        self.xam_path = xam_path
        self.ref_name = ref_name
        self.end5 = end5
        self.end3 = end3
        self.spanning = spanning
        self.min_qual = min_qual
        self.owner = owner
        self._sam_path: OneRefAlignmentTempFilePath | None = None
        self._sam_file: BufferedReader | None = None
    
    def __enter__(self):
        # Convert the BAM file to a temporary SAM file
        if self.owner:
            if self.spanning:
                selector = None
                xam_selected = self.xam_path
            else:
                selector = BamVectorSelector(self.top_dir,
                                             self.max_cpus,
                                             self.xam_path,
                                             self.ref_name,
                                             self.end5,
                                             self.end3)
                xam_selected = selector.run()
            sorter = SamVectorSorter(self.top_dir, self.max_cpus, xam_selected)
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
        if not self._sam_path:
            raise RuntimeError(f"{self} has no open SAM path.")
        return self._sam_path

    @property
    def sam_file(self):
        if not self._sam_file:
            raise RuntimeError(f"{self} has no open SAM file.")
        return self._sam_file

    def _seek_beginning(self):
        self.sam_file.seek(0)
    
    @cached_property
    @_reset_seek
    def _rec1_pos(self):
        self._seek_beginning()
        while (line := self.sam_file.readline()).startswith(SAM_HEADER):
            pass
        return self.sam_file.tell() - len(line)

    def _seek_rec1(self):
        self.sam_file.seek(self._rec1_pos)
        
    @cached_property
    @_reset_seek
    def paired(self):
        self._seek_rec1()
        first_line = self.sam_file.readline()
        return (SamRead(first_line, self.min_qual).flag.paired if first_line
                else False)
    
    @_range_of_records
    @staticmethod
    def _get_records_single(samview: SamViewer):
        while read := SamRead(samview.sam_file.readline(), samview.min_qual):
            if read.flag.paired:
                raise ValueError("Found paired-end read in single-end SAM file.")
            yield SamRecord(read)

    @_range_of_records
    @staticmethod
    def _get_records_paired_flexible(samview: SamViewer):
        prev_read: Optional[SamRead] = None
        while line := samview.sam_file.readline():
            read = SamRead(line, samview.min_qual)
            assert read.flag.paired
            if prev_read:
                # The previous read has not yet been yielded
                if prev_read.qname == read.qname:
                    # The current read is the mate of the previous read
                    if prev_read.flag.end5:
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
    @staticmethod
    def _get_records_paired_strict(samview: SamViewer):
        while line := samview.sam_file.readline():
            yield SamRecord(SamRead(line,
                                    samview.min_qual),
                            SamRead(samview.sam_file.readline(),
                                    samview.min_qual))
    
    def get_records(self, start: int, stop: int):
        if self.paired:
            if self.spanning:
                record_generator = self.__class__._get_records_paired_strict
            else:
                record_generator = self.__class__._get_records_paired_flexible
        else:
            record_generator = self.__class__._get_records_single
        return record_generator(self, start, stop)
    
    @_reset_seek
    def get_batch_indexes(self, batch_size: int):
        if batch_size <= 0:
            raise ValueError("batch_size must be a positive integer")
        n_skip = (self.paired + 1) * (batch_size - 1)
        self._seek_rec1()
        while True:
            yield self.sam_file.tell()
            try:
                line = next(self.sam_file)
            except StopIteration:
                break
            else:
                for _, line in zip(range(n_skip), self.sam_file):
                    pass
                if self.paired:
                    # Compare the current and next query names.
                    try:
                        qname = line.split()[0]
                    except IndexError:
                        pass
                    else:
                        line_next = self.sam_file.readline()
                        try:
                            qname_next = line_next.split()[0]
                        except IndexError:
                            pass
                        else:
                            if qname != qname_next:
                                # If the current and next query names differ
                                # (the lines do not come from two paired mates),
                                # then backtrack to the beginning of line_next.
                                self.sam_file.seek(-len(line_next), 1)
