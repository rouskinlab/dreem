from __future__ import annotations
from functools import cached_property, wraps
from io import BufferedReader
from typing import Callable, Optional

from dreem.align.reads import BamVectorSelector, SamVectorSorter
from dreem.util.path import TopDirPath, OneRefAlignmentInFilePath, OneRefAlignmentStepFilePath
from dreem.vector.vector import *


def _reset_seek(func: Callable):
    @wraps(func)
    def wrapper(reader: SamReader, *args, **kwargs):
        prev_pos = reader.sam_file.tell()
        result = func(reader, *args, **kwargs)
        reader.sam_file.seek(prev_pos)
        return result
    return wrapper


def _range_of_records(func: Callable):
    @wraps(func)
    @_reset_seek
    def wrapper(reader: SamReader, start: int, stop: int | None):
        reader.sam_file.seek(start)
        records = iter(func(reader))
        if stop is None:
            return records
        else:
            while True:
                if reader.sam_file.tell() < stop:
                    yield next(records)
                else:
                    assert reader.sam_file.tell() == stop
                    break
    return wrapper


@_range_of_records
def _get_records_single(reader: SamReader):
    while read := SamRead(reader.sam_file.readline(), reader.min_qual):
        if read.flag.paired:
            raise ValueError("Found paired-end read in single-end SAM file.")
        yield SamRecord(read)


@_range_of_records
def _get_records_paired_subspan(reader: SamReader):
    prev_read: Optional[SamRead] = None
    while line := reader.sam_file.readline():
        read = SamRead(line, reader.min_qual)
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
                # The previous read is paired, but its mate is not
                # in the SAM file. This situation can happen when
                # only the reads mapping to a specific region are
                # exported to a SAM file and only one mate overlaps
                # that region. In this case, it is valid for only
                # one of the two paired mates to be given. Thus,
                # strict mode (which normally ensures that read 2
                # is given if read 1 is paired) is turned off. The
                # given mate is processed as if a single-end read.
                yield SamRecord(prev_read, strict=False)
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
def _get_records_paired_spanning(reader: SamReader):
    while line := reader.sam_file.readline():
        yield SamRecord(SamRead(line, reader.min_qual),
                        SamRead(reader.sam_file.readline(), reader.min_qual))


class SamReader(object):
    def __init__(self,
                 temp_dir: TopDirPath,
                 save_temp: bool,
                 resume: bool,
                 n_procs: int,
                 xam_path: OneRefAlignmentInFilePath,
                 ref_name: str,
                 end5: int,
                 end3: int,
                 spanning: bool,
                 min_qual: int,
                 owner: bool = True):
        self.temp_dir = temp_dir
        self.save_temp = save_temp
        self.resume = resume
        self.max_procs = n_procs
        self.xam_path = xam_path
        self.ref_name = ref_name
        self.end5 = end5
        self.end3 = end3
        self.spanning = spanning
        self.min_qual = min_qual
        self.owner = owner
        self._sam_sorter: SamVectorSorter | None = None
        self._sam_path: OneRefAlignmentStepFilePath | None = None
        self._sam_file: BufferedReader | None = None
    
    def __enter__(self):
        # Convert the BAM file to a temporary SAM file
        if self.owner:
            selector = BamVectorSelector(top_dir=self.temp_dir,
                                         save_temp=self.save_temp,
                                         resume=self.resume,
                                         num_cpus=self.max_procs,
                                         xam=self.xam_path,
                                         ref=self.ref_name,
                                         end5=self.end5,
                                         end3=self.end3)
            if self.spanning:
                xam_selected = self.xam_path
            elif selector.output.path.is_file() and self.resume:
                xam_selected = selector.output
            else:
                xam_selected = selector.run()
            sorter = SamVectorSorter(top_dir=self.temp_dir,
                                     save_temp=self.save_temp,
                                     resume=self.resume,
                                     num_cpus=self.max_procs,
                                     xam=xam_selected)
            if sorter.output.path.is_file() and self.resume:
                self._sam_path = sorter.output
            else:
                self._sam_path = sorter.run(name=True)
            if not self.save_temp:
                selector.clean()
        else:
            self._sam_path = self.xam_path
        self._sam_file = open(self.sam_path.path, "rb")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self._sam_file.close()
        self._sam_file = None
        if self._sam_sorter is not None and not self.save_temp:
            self._sam_sorter.clean()
        if self.owner and not self.save_temp:
            self.sam_path.path.unlink(missing_ok=True)
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
    
    def get_records(self, start: int, stop: int):
        # Each record_generator method is obtained by self.__class__
        # instead of just self because the method is static
        if self.paired:
            if self.spanning:
                return _get_records_paired_spanning(self, start, stop)
            return _get_records_paired_subspan(self, start, stop)
        return _get_records_single(self, start, stop)
    
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
