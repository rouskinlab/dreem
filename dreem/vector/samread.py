from __future__ import annotations
from functools import cached_property, wraps
from io import BufferedReader
import logging
import pathlib
from typing import Callable

from ..align.reads import (SAM_DELIMITER, SAM_HEADER, sort_xam, view_xam)
from ..util import path


logger = logging.getLogger(__name__)


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
                    if reader.sam_file.tell() > stop:
                        raise ValueError(f"SAM reader {reader} stopped at "
                                         f"{reader.sam_file.tell()} but "
                                         f"expected {stop}")
                    break
    return wrapper


@_range_of_records
def _iter_records_single(reader: SamReader):
    """ Yield the read name and line for every read in the file. """
    while line := reader.sam_file.readline():
        yield line.split(SAM_DELIMITER, 1)[0], line, b""


@_range_of_records
def _iter_records_paired_lenient(reader: SamReader):
    prev_line: bytes = b""
    prev_name: bytes = b""
    while line := reader.sam_file.readline():
        if prev_line:
            # Read name is the first field of the line.
            name = line.split(SAM_DELIMITER, 1)[0]
            # The previous read has not yet been yielded.
            if prev_name == name:
                # The current read is the mate of the previous read.
                yield prev_name, prev_line, line
                prev_line = b""
                prev_name = b""
            else:
                # The previous read is paired, but its mate is not
                # in the SAM file. This situation can happen when
                # only the reads mapping to a specific section are
                # exported to a SAM file and only one mate overlaps
                # that section. In this case, it is valid for only
                # one of the two paired mates to be given. Thus,
                # strict mode (which normally ensures that read 2
                # is given if read 1 is paired) is turned off. The
                # given mate is processed as if a single-end read.
                yield prev_name, prev_line, b""
                # Save the current read so that if its mate is the next
                # read, it will be returned as a pair.
                prev_line = line
                prev_name = name
        else:
            # Save the current read so that if its mate is the next
            # read, it will be returned as a pair.
            prev_line = line
            prev_name = line.split(SAM_DELIMITER, 1)[0]
    if prev_line:
        # In case the last read has not yet been yielded, do so.
        yield prev_name, prev_line, b""


@_range_of_records
def _iter_records_paired_strict(reader: SamReader):
    """ Yield the common name and both lines for every pair of reads in
    the file. """
    while line := reader.sam_file.readline():
        yield line.split(SAM_DELIMITER, 1)[0], line, reader.sam_file.readline()


class SamReader(object):
    """ Interface to a SAM file. Extract reads from a specific section
    of a BAM file into a SAM file, then return a SamRecord object for
    each record in the file. """

    def __init__(self, /,
                 xam_inp: (path.OneRefAlignmentInFilePath |
                           path.SectionAlignmentTempFilePath), *,
                 temp_dir: pathlib.Path,
                 end5: int,
                 end3: int,
                 isfullref: bool,
                 n_procs: int,
                 save_temp: bool):
        self.xam_inp = xam_inp
        self.temp_dir = temp_dir
        self.end5 = end5
        self.end3 = end3
        self.isfullref = isfullref
        self.n_procs = n_procs
        self.save_temp = save_temp
        self.bam_split: pathlib.Path | None = None
        self.sam_path: path.SectionAlignmentTempFilePath | None = None
        self.sam_file: BufferedReader | None = None
    
    def __enter__(self):
        # Convert the BAM file to a temporary SAM file
        if self.is_owner:
            # Output only the reads aligning to the section of interest.
            # This step is unnecessary if the section is the full ref.
            if not self.isfullref:
                self.bam_split = path.SectionAlignmentTempFilePath(
                    top=str(self.temp_dir),
                    module=path.Module.VECTOR,
                    step=path.Step.VECTOR_SELECT,
                    sample=self.xam_inp.sample,
                    ref=self.xam_inp.ref,
                    end5=self.end5,
                    end3=self.end3,
                    ext=path.BAM_EXT).path
                view_xam(self.xam_inp.path, self.bam_split,
                         ref=self.xam_inp.ref, end5=self.end5, end3=self.end3)
            # Sort the alignment records by name.
            # Note: this step is unnecessary for single-end reads, but
            # there is currently no way to check whether the reads are
            # single- or paired-end before this step.
            self.sam_path = path.SectionAlignmentTempFilePath(
                top=str(self.temp_dir),
                module=path.Module.VECTOR,
                step=path.Step.VECTOR_SORT,
                sample=self.xam_inp.sample,
                ref=self.xam_inp.ref,
                end5=self.end5,
                end3=self.end3,
                ext=path.SAM_EXT)
            sort_xam(self.xam_inp if self.bam_split is None else self.bam_split,
                     self.sam_path.path,
                     name=True)
            # Delete the temporary BAM file, if it exists.
            if not self.save_temp and self.bam_split is not None:
                self.bam_split.unlink(missing_ok=True)
                self.bam_split = None
        else:
            # If not the owner of the SAM file, then borrow it pre-made.
            self.sam_path = self.xam_inp
        # Open the sorted, section-specific SAM file for reading.
        self.sam_file = open(self.sam_path.path, "rb")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.sam_file.close()
        self.sam_file = None
        if self.is_owner and not self.save_temp:
            # Delete the SAM file, if its owner.
            self.sam_path.path.unlink(missing_ok=True)
            self.sam_path = None

    @property
    def is_owner(self):
        return isinstance(self.xam_inp, path.OneRefAlignmentInFilePath)

    def _seek_beginning(self):
        """ Seek to the beginning of the SAM file. """
        self.sam_file.seek(0)
    
    @cached_property
    @_reset_seek
    def _rec1_pos(self):
        """ Return the position of the first record in the SAM file. """
        self._seek_beginning()
        while (line := self.sam_file.readline()).startswith(SAM_HEADER):
            pass
        return self.sam_file.tell() - len(line)

    def _seek_rec1(self):
        """ Seek to the first record in the SAM file. """
        self.sam_file.seek(self._rec1_pos)

    @cached_property
    @_reset_seek
    def paired(self):
        """ Return whether the reads in the SAM file are paired-end. """
        self._seek_rec1()
        first_line = self.sam_file.readline()
        try:
            flag = first_line.split(SAM_DELIMITER, 2)[1]
            paired = bool(int(flag) & 1)
            logger.debug(f"SAM file {self.sam_path} has "
                         f"{'paired' if paired else 'single'}-ended reads")
            return paired
        except (IndexError, ValueError):
            raise ValueError("Failed to determine pairing from first record: "
                             + first_line.decode())

    @property
    def mates_per_record(self):
        """ Mates per record: 1 if single-end, 2 if paired-end. """
        return self.paired + 1
    
    def iter_records(self, start: int, stop: int, strict_pairs: bool):
        """ Return an iterator of all records in the SAM file. """
        if self.paired:
            if strict_pairs:
                if self.isfullref:
                    return _iter_records_paired_strict(self, start, stop)
                logging.warning(f"Disabling strict pairs for {self} because "
                                f"it is not the full reference sequence")
            return _iter_records_paired_lenient(self, start, stop)
        return _iter_records_single(self, start, stop)
    
    @_reset_seek
    def iter_batch_indexes(self, records_per_batch: int):
        """ Yield the start and end positions of every batch in the SAM
        file, where each batch should have about ```records_per_batch```
        records. Assume that for nearly all records in paired-end SAM
        files, both mates are present. In the extreme case that only one
        mate is present for every paired-end record, there can be up to
        ```2 * records_per_batch``` records in a batch. """
        if records_per_batch <= 0:
            raise ValueError("records_per_batch must be a positive integer")
        # Number of lines to skip between batches: the number of records
        # per batch minus one (to account for the one line that is read
        # at the beginning of each batch, which ensures that every batch
        # has at least one line) times the number of mates per record.
        n_skip = (records_per_batch - 1) * (self.paired + 1)
        # Start at the beginning of the first record.
        self._seek_rec1()
        # Store sam_file as a local variable to speed up the iteration.
        sam_file = self.sam_file
        # Record the start position of the first batch.
        start_of_batch = sam_file.tell()
        # Yield batches until the SAM file is exhausted. If there are no
        # records in the file, then this loop will exit immediately.
        while line := sam_file.readline():
            # The start position of the batch is the beginning of the
            # line that was just read. Since the position in the file is
            # now the end of that line, the line's length is subtracted.
            if start_of_batch != sam_file.tell() - len(line):
                # The above should always equal; just a sanity check.
                raise ValueError(f"Inconsistent batch start: {start_of_batch} "
                                 f"≠ {sam_file.tell() - len(line)}")
            # Read either the prescribed number of lines or to the end
            # of the file, whichever limit is reached first.
            for _, line in zip(range(n_skip), sam_file, strict=False):
                pass
            if self.paired:
                # If the reads are paired-end, then assume that the vast
                # majority of records have two lines. If every record
                # has two lines, then the current position in the file
                # is the end of the first line / beginning of the second
                # line of a record, because an odd number of lines have
                # been read, including via line := sam_file.readline().
                # Thus, check whether the read name of the line that was
                # just read matches the read name of the next line.
                try:
                    name_prev = line.split(SAM_DELIMITER, 1)[0]
                except IndexError:
                    raise ValueError(f"Cannot read SAM line: '{line.decode()}'")

                else:
                    line_next = self.sam_file.readline()
                    try:
                        qname_next = line_next.split(SAM_DELIMITER, 1)[0]
                    except IndexError:
                        pass
                    else:
                        if qname != qname_next:
                            # If the current and next query names
                            # differ, then the lines do not come
                            # from two paired mates, so backtrack
                            # to the beginning of line_next.
                            self.sam_file.seek(-len(line_next), 1)

    def list_batch_indexes(self, vectors_per_batch: int):
        indexes = list(self.iter_batch_indexes(vectors_per_batch))


    def __str__(self):
        return str(self.sam_path)
