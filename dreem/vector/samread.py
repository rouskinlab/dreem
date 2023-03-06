from __future__ import annotations
from functools import cached_property, wraps
from io import BufferedReader
import logging
from typing import Callable, Optional

from ..align.reads import BamVectorSelector, SamVectorSorter, SAM_HEADER
from ..util.path import OneRefAlignmentInFilePath, OneRefAlignmentStepFilePath
from ..vector.vector import *


def _reset_seek(func: Callable):
    @wraps(func)
    def wrapper(reader: SamReader, *args, **kwargs):
        prev_pos = reader.sam_file.tell()
        result = func(reader, *args, **kwargs)
        reader.sam_file.seek(prev_pos)
        return result
    return wrapper


def _range_of_records(get_records_func: Callable):
    @wraps(get_records_func)
    @_reset_seek
    def wrapper(reader: SamReader, start: int, stop: int | None):
        reader.sam_file.seek(start)
        records = get_records_func(reader)
        if stop is None:
            yield from records
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
    while read := SamRead(reader.sam_file.readline()):
        if read.flag.paired:
            logging.error(f"Got paired-end read {read} in single-end SAM file")
        else:
            yield SamRecord(read)


@_range_of_records
def _get_records_paired_lenient(reader: SamReader):
    prev_read: Optional[SamRead] = None
    while line := reader.sam_file.readline():
        read = SamRead(line)
        if not read.flag.paired:
            logging.error(f"Got single-end read {read} in paired-end SAM file")
            continue
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
def _get_records_paired_strict(reader: SamReader):
    while line := reader.sam_file.readline():
        yield SamRecord(SamRead(line),
                        SamRead(reader.sam_file.readline()),
                        strict=True)


class SamReader(object):
    """ Interface to a SAM file. Extract reads from a specific region
    of a BAM file into a SAM file, then return a SamRecord object for
    each record in the file. """

    def __init__(self,
                 xam_path: OneRefAlignmentInFilePath,
                 temp_dir: str,
                 end5: int,
                 end3: int,
                 spans: bool,
                 n_procs: int,
                 save_temp: bool,
                 resume: bool,
                 owner: bool):
        self.xam_path = xam_path
        self.temp_dir = temp_dir
        self.end5 = end5
        self.end3 = end3
        self.spans = spans
        self.n_procs = n_procs
        self.save_temp = save_temp
        self.resume = resume
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
                                         n_procs=self.n_procs,
                                         input_path=self.xam_path,
                                         ref=self.ref,
                                         end5=self.end5,
                                         end3=self.end3)
            if selector.output.path.is_file() and self.resume:
                xam_selected = selector.output
            else:
                xam_selected = selector.run()
            sorter = SamVectorSorter(top_dir=self.temp_dir,
                                     save_temp=self.save_temp,
                                     resume=self.resume,
                                     n_procs=self.n_procs,
                                     input_path=xam_selected)
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
    def sample(self):
        return self.xam_path.sample

    @property
    def ref(self):
        return self.xam_path.ref

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
        return SamRead(first_line).flag.paired if first_line else False
    
    def get_records(self, start: int, stop: int, strict_pairs: bool):
        # Each record_generator method is obtained by self.__class__
        # instead of just self because the method is static
        if self.paired:
            if strict_pairs:
                if self.spans:
                    return _get_records_paired_strict(self, start, stop)
                logging.warning(f"Disabling strict pairs for {self} because "
                                f"it does not span the reference sequence")
            return _get_records_paired_lenient(self, start, stop)
        return _get_records_single(self, start, stop)
    
    @_reset_seek
    def get_batch_indexes(self, vectors_per_batch: int):
        if vectors_per_batch <= 0:
            raise ValueError("vectors_per_batch must be a positive integer")
        n_skip = (self.paired + 1) * (vectors_per_batch - 1)
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

    def __str__(self):
        return f"{self.sample}@{self.ref}:{self.end5}-{self.end3}"
