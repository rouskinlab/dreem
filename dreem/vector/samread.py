from __future__ import annotations
from functools import cached_property, wraps
from io import BufferedReader
from logging import getLogger
from pathlib import Path
from typing import Callable

from ..align.xams import SAM_DELIMITER, SAM_HEADER, sort_xam, view_xam
from ..util import path

logger = getLogger(__name__)


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
    def wrapper(reader: SamReader, start: int, stop: int):
        logger.debug(f"Reading records from {reader}, starting at {start} and "
                     f"stopping at {stop}")
        reader.sam_file.seek(start)
        records = get_records_func(reader)
        n_records = 0
        while True:
            if reader.sam_file.tell() < stop:
                n_records += 1
                yield next(records)
            else:
                if reader.sam_file.tell() > stop:
                    raise ValueError(f"SAM reader {reader} stopped at "
                                     f"{reader.sam_file.tell()} but "
                                     f"expected {stop}")
                break
        logger.debug(f"Read {n_records} records from {reader}, starting at "
                     f"{start} and stopping at {stop}")
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
                 temp_dir: Path,
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
        self.n_padd = n_procs - 1
        self.save_temp = save_temp
        self.bam_split: Path | None = None
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
                view_xam(self.xam_inp.path, self.bam_split, n_padd=self.n_padd,
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
                     self.sam_path.path, name=True, n_padd=self.n_padd)
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
                             + first_line.decode()+f" {self.sam_path}")

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
                logger.warning(f"Disabling strict pairs for {self} because "
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
        logger.debug(f"Computing batches for {self} (paired={self.paired}), "
                     f"targeting {records_per_batch} records per batch")
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
        # Yield batches until the SAM file is exhausted. If there are no
        # records in the file, then this loop will exit immediately.
        batch = 0
        while line := sam_file.readline():
            # The start position of the batch is the beginning of the
            # line that was just read. Since the position in the file is
            # now the end of that line, the line's length is subtracted.
            start_of_batch = sam_file.tell() - len(line)
            # Read either the prescribed number of lines or to the end
            # of the file, whichever limit is reached first.
            for _, line in zip(range(n_skip), sam_file, strict=False):
                pass
            if self.paired:
                # Get the end position of the current line.
                end_curr = sam_file.tell()
                # Try to get the name of the read in the current line.
                try:
                    name_curr = line.split(SAM_DELIMITER, 1)[0]
                except IndexError:
                    # Getting the read name failed, so it is impossible
                    # to tell whether this read is the mate of the next
                    # read. End this batch at the current read so that
                    # the next iteration can check whether the next read
                    # and the read after the next read are mates.
                    logger.error(f"Cannot find read name in '{line.decode()}'")
                    logger.debug(f"Batch {batch}: {start_of_batch}-{end_curr}")
                    yield start_of_batch, end_curr
                    continue
                # Try to read the next line.
                try:
                    line = next(sam_file)
                except StopIteration:
                    # The end of the file has been reached, so the end
                    # of the current line is the end of the batch.
                    logger.debug(f"Batch {batch}: {start_of_batch}-{end_curr}")
                    yield start_of_batch, end_curr
                    continue
                # Get the end position of the next line.
                end_next = sam_file.tell()
                # Try to get the name of the read in the next line.
                try:
                    name_next = line.split(SAM_DELIMITER, 1)[0]
                except IndexError:
                    # Getting the read name failed, so it is impossible
                    # to tell whether this next read is the mate of the
                    # current read. End this batch at the next read so
                    # that this batch will have an even number of lines,
                    # which, when in doubt, is best for paired reads.
                    logger.error(f"Cannot find read name in '{line.decode()}'")
                    logger.debug(f"Batch {batch}: {start_of_batch}-{end_next}")
                    yield start_of_batch, end_next
                    continue
                # Determine the end of the batch based on the names of
                # the current and the next reads.
                if name_curr == name_next:
                    # If the read names match, then they are mates and
                    # should be part of the same batch. So end the batch
                    # at the end of the next read.
                    logger.debug(f"Batch {batch}: {start_of_batch}-{end_next}")
                    yield start_of_batch, end_next
                else:
                    # Otherwise, the reads are not mates. If the batch
                    # were to end at the end of the next read, then it
                    # would be impossible for the next read to be in the
                    # same batch as the read after the next read, even
                    # if those two reads were mates. To allow those two
                    # reads to be compared and potentially be placed in
                    # the same batch, the current batch must end after
                    # the current read.
                    logger.debug(f"Batch {batch}: {start_of_batch}-{end_curr}")
                    yield start_of_batch, end_curr
                    # Because each iteration of the loop must start at
                    # the position at which the last iteration ended,
                    # move back to that position in the file.
                    sam_file.seek(end_curr)
            else:
                # If the file does not contain paired-end reads, then
                # the end of the batch is the end of the current read.
                end_curr = sam_file.tell()
                logger.debug(f"Batch {batch}: {start_of_batch}-{end_curr}")
                yield start_of_batch, end_curr
            # Increment the batch number.
            batch += 1

    def list_batch_indexes(self, records_per_batch: int):
        """ Return lists of the start and stop indexes of each batch in
        the SAM file. """
        indexes = SamReader.iter_batch_indexes(self, records_per_batch)
        starts, stops = map(list, zip(*indexes, strict=True))
        return starts, stops

    def __str__(self):
        return str(self.sam_path)
