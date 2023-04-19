from __future__ import annotations

from concurrent.futures import Future, ProcessPoolExecutor
from datetime import datetime
from functools import cached_property, partial
from itertools import starmap as itsmap
from logging import getLogger
from multiprocessing import Pool
from pathlib import Path
from sys import byteorder

from ..align.xams import view_xam
from ..util import path
from ..util.parallel import get_num_parallel
from ..util.seq import DNA, parse_fasta, NOCOV
from ..util.files import digest_file
from .batch import BATCH_NUM_START, write_batch, mib_to_bytes
from .cigarray import vectorize_line, vectorize_pair, VectorError
from .report import VectorReport, _get_report_path
from .sam import iter_batch_indexes, iter_records

logger = getLogger(__name__)


def _vectorize_record(read_name: bytes, line1: bytes, line2: bytes, *,
                      blank: bytearray, ref_seq: bytes, ref: str,
                      min_qual: int, ambid: bool):
    """ Compute the mutation vector of a record in a SAM file. """
    # Initialize a blank mutation vector.
    muts = blank.copy()
    # Fill the mutation vector with data from the SAM line(s).
    if line2:
        vectorize_pair(line1, line2, muts, ref_seq, len(ref_seq),
                       ref, min_qual, ambid)
    else:
        # Using seq instead of byteseq crashes vectoring.
        vectorize_line(line1, muts, ref_seq, len(ref_seq),
                       ref, min_qual, ambid)
    # Check whether the mutation vector is still blank.
    if muts == blank:
        raise VectorError(f"Mutation vector is blank")
    return read_name, muts


def _vectorize_batch(batch: int,
                     start: int,
                     stop: int,
                     *,
                     temp_sam: Path,
                     out_dir: Path,
                     sample: str,
                     ref: str,
                     ref_seq: bytes,
                     min_qual: int,
                     ambid: bool):
    """ Compute mutation vectors for every SAM record in one batch,
    write the vectors to a batch file, and return its MD5 checksum
    and the number of vectors. """
    n_reads = 0
    blank = bytearray(NOCOV.to_bytes(1, byteorder)) * len(ref_seq)
    try:
        logger.debug(f"Began vectorizing batch {batch} of {temp_sam} "
                     f"({start} - {stop})")

        # Wrap self._vectorize_record with keyword arguments and a
        # try-except block so that if one record fails to vectorize,
        # it does not crash all the others.

        def vectorize_record(read_name: bytes, line1: bytes, line2: bytes):
            try:
                return _vectorize_record(read_name, line1, line2,
                                         blank=blank, ref_seq=ref_seq, ref=ref,
                                         min_qual=min_qual, ambid=ambid)
            except Exception as err:
                logger.error(
                    f"Failed to vectorize read '{read_name.decode()}': {err}")
                return b"", bytearray()

        with open(temp_sam, "rb") as sam_file:
            # Vectorize every record in the batch.
            records = iter_records(sam_file, start, stop)
            read_names, muts = zip(*itsmap(vectorize_record, records))
            # For every read for which creating a mutation vector
            # failed, an empty string was returned as the read name
            # and an empty bytearray as the mutation vector. The
            # empty read names must be filtered out, while the empty
            # mutation vectors will not cause problems because,
            # being of length zero, they will effectively disappear
            # when all the vectors are concatenated into a 1D array.
            read_names = list(filter(None, read_names))
        # Compute the number of reads that passed and failed.
        n_reads = len(muts)
        n_pass = len(read_names)
        n_fail = n_reads - n_pass
        if not n_pass:
            logger.warning(f"Batch {batch} of {temp_sam} yielded 0 vectors")
        # Write the names and vectors to a file.
        batch_file = write_batch(batch, muts, read_names,
                                 sample=sample, ref=ref,
                                 seq=ref_seq, out_dir=out_dir)
        # Compute the MD5 checksum of the file.
        checksum = digest_file(batch_file)
        logger.debug(f"Ended vectorizing batch {batch} of {temp_sam} "
                     f"({start} - {stop}): {n_pass} pass, {n_fail} fail")
        return n_pass, n_fail, checksum
    except Exception as error:
        logger.critical(
            f"Failed to generate batch {batch} of {temp_sam}: {error}")
        return 0, n_reads, ""


class VectorWriter(object):
    """
    Compute mutation vectors for all reads from one sample mapping to
    one reference sequence.
    """

    def __init__(self, /, bam_file: Path, seq: DNA):
        self.bam = bam_file
        self.seq = bytes(seq)

    @cached_property
    def sample_ref(self):
        fields = path.parse(self.bam, path.SampSeg, path.XamSeg)
        return fields[path.SAMP], fields[path.REF]

    @property
    def sample(self):
        return self.sample_ref[0]

    @property
    def ref(self):
        return self.sample_ref[1]

    def outputs_valid(self, /, out_dir: Path):
        """ Return whether the report file exists and, if so, whether
        all batch files of mutation vectors listed in the report exist
        and match the checksums recorded in the report. """
        try:
            VectorReport.load(_get_report_path(out_dir, self.sample, self.ref))
        except (FileNotFoundError, ValueError):
            return False
        else:
            return True

    def _write_report(self, /, *, out_dir: Path, **kwargs):
        logger.info(f"Began writing report of {self}")
        report = VectorReport(seq=DNA(self.seq),
                              sample=self.sample,
                              ref=self.ref,
                              **kwargs)
        report_path = report.save(out_dir)
        logger.info(f"Ended writing report of {self} to {report_path}")
        return report_path

    def _vectorize_bam(self, /, *,
                       out_dir: Path,
                       temp_dir: Path,
                       save_temp: bool,
                       batch_size: int,
                       n_procs: int,
                       phred_enc: int,
                       min_phred: int,
                       ambid: bool):
        """ Compute a mutation vector for every record in a BAM file,
        split among one or more batches. For each batch, write the
        vectors to one batch file, and compute its checksum. """
        # Open the primary SAM file reader to write the subset of SAM
        # records to a temporary SAM file and determine the number and
        # start/stop indexes of each batch of records in the file.
        # The SAM file will remain open until exiting the with block.
        logger.info(f"Began vectorizing {self}")
        # Determine the path of the temporary SAM file.
        temp_sam = path.build(path.StepSeg, path.SampSeg, path.XamSeg,
                              top=temp_dir, step=path.STEPS_VECT[0],
                              sample=self.sample, ref=self.ref,
                              ext=path.SAM_EXT)
        # Create the temporary SAM file.
        view_xam(self.bam, temp_sam, n_procs=n_procs)
        sam_file = open(temp_sam, "rb")
        try:
            # Compute number of records per batch.
            recs_per_batch = max(1, mib_to_bytes(batch_size) // len(self.seq))
            # Create a generator for start and stop indexes of batches.
            indexes = enumerate(iter_batch_indexes(sam_file, recs_per_batch),
                                start=BATCH_NUM_START)
            # Wrap _vectorize_batch with keyword arguments.
            vb_wrapper = partial(_vectorize_batch, temp_sam=temp_sam,
                                 out_dir=out_dir, sample=self.sample,
                                 ref=self.ref, ref_seq=self.seq, ambid=ambid,
                                 min_qual=get_min_qual(min_phred, phred_enc))
            if n_procs > 1:
                # Process batches of records simultaneously in parallel.
                logger.debug(f"Initializing pool of {n_procs} processes")
                with ProcessPoolExecutor(max_workers=n_procs) as pool:
                    logger.debug(f"Opened process pool of at most {n_procs} "
                                 f"processes for {self}")
                    futures: list[Future] = list()
                    for batch, (start, stop) in indexes:
                        futures.append(pool.submit(vb_wrapper,
                                                   batch, start, stop))
                        logger.debug(f"Submitted batch {batch} to process pool "
                                     f"for {self}")
                    logger.debug(f"Waiting until all {len(futures)} processes "
                                 f"for {self} finish")
                    results = [future.result() for future in futures]
                logger.debug(f"Closed process pool for {self}")
            else:
                # Process batches of records one at a time in series.
                results = [vb_wrapper(batch, start, stop)
                           for batch, (start, stop) in indexes]
            # The list of results contains, for each batch, a tuple of the
            # number of mutation vectors in the batch and the MD5 checksum
            # of the batch. Compute the total number of vectors and list all
            # the checksums.
            n_pass = 0
            n_fail = 0
            checksums: list[str] = list()
            for p, f, c in results:
                n_pass += p
                n_fail += f
                checksums.append(c)
            logger.info(
                f"Ended vectorizing {self}: {n_pass} pass, {n_fail} fail")
            return n_pass, n_fail, checksums
        finally:
            sam_file.close()
            if not save_temp:
                # Delete the temporary SAM file before exiting.
                temp_sam.unlink(missing_ok=True)

    def vectorize(self, /, *, rerun: bool, out_dir: Path, **kwargs):
        """ Compute a mutation vector for every record in a BAM file,
        write the vectors into one or more batch files, compute their
        checksums, and write a report summarizing the results. """
        report_path = _get_report_path(out_dir, self.sample, self.ref)
        try:
            VectorReport.load(report_path)
        except (FileNotFoundError, ValueError):
            # Vectorization has not yet been run.
            pass
        else:
            # Vectorization has already been run.
            if not rerun:
                logger.warning(f"Skipping vectorization of {self} because "
                               f"all output files already exist")
                return report_path
        # Compute the mutation vectors, write them to batch files, and
        # generate a report.
        # Vectorize the BAM file and time how long it takes.
        began = datetime.now()
        n_pass, n_fail, checksums = self._vectorize_bam(out_dir=out_dir,
                                                        **kwargs)
        ended = datetime.now()
        # Write a report of the vectorization.
        written = self._write_report(out_dir=out_dir,
                                     n_vectors=n_pass,
                                     n_readerr=n_fail,
                                     checksums=checksums,
                                     began=began,
                                     ended=ended)
        if written != report_path:
            raise ValueError("Intended and actual paths of report differ: "
                             f"{report_path} ≠ {written}")
        return report_path

    def __str__(self):
        return f"the vectorization of {self.bam}"


def get_min_qual(min_phred: int, phred_enc: int):
    """
    Return the minimum quality for a base in a read to be considered
    informative, as the ASCII integer corresponding to the character
    in the FASTQ file that is the minimum valid quality.

    Return
    ------
    int
        The ASCII value corresponding to the character in the FASTQ
        file read quality string that represents the minimum quality
        for a base to be considered informative

    Examples
    --------
    For example, if the minimum Phred score (```min_phred```) that
    is accepted as informative is 20, and the Phred encoding of the
    FASTQ file (```phred_enc```) is 33 (i.e. ASCII+33), then the
    minimum quality as an ASCII integer (```min_qual```) is 20 + 33
    = 53, which is character '5'. If ```min_phred``` were 37, then
    ```min_qual``` would be 37 + 33 = 70, which is character 'F'.
    """
    return min_phred + phred_enc


def get_writers(fasta: Path, bam_files: list[Path]):
    logger.info("Began creating vector writers")
    ref_seqs = dict(parse_fasta(fasta))
    writers: list[VectorWriter] = list()
    for bam_file in set(bam_files):
        try:
            # Parse the fields of the input BAM file.
            ref = path.parse(bam_file, path.SampSeg, path.XamSeg)[path.REF]
            try:
                seq = ref_seqs[ref]
            except KeyError:
                logger.critical(
                    f"Reference '{ref}' for {bam_file} does not exist")
                continue
            writers.append(VectorWriter(bam_file, seq))
            logger.debug(f"Created vector writer for {bam_file}")
        except Exception as error:
            logger.critical(
                f"Failed to create vector writer for {bam_file}: {error}")
    logger.info(f"Ended creating {len(writers)} vector writers")
    return writers


def write(writers: list[VectorWriter], *,
          phred_enc: int,
          min_phred: int,
          max_procs: int,
          parallel: bool,
          **kwargs) -> tuple[str, ...]:
    """ Generate mutational profiles of one or more vector writers. """
    logger.info("Began generating mutational profiles")
    # Determine method of parallelization. Do not use hybrid mode, which
    # would try to process multiple SAM files in parallel and use more
    # than one processe for each file. Python ```multiprocessing.Pool```
    # forbids a daemon process (one for each SAM file in parallel) from
    # spawning additional processes.
    n_tasks_parallel, n_procs_per_task = get_num_parallel(len(writers),
                                                          max_procs,
                                                          parallel)

    # Wrap VectorWriter.vectorize in a try-except block to catch all
    # exceptions, so that if one VectorWriter crashes, it does not
    # crash all the others.

    def vectorize(writer: VectorWriter):
        try:
            return writer.vectorize(phred_enc=phred_enc,
                                    min_phred=min_phred,
                                    n_procs=n_procs_per_task,
                                    **kwargs)
        except Exception as error:
            # Alert that vectoring failed and return no report path.
            logger.critical(f"Failed to vectorize {writer}: {error}")
            return None

    # Call the vectorize method of each writer, passing args.
    if n_tasks_parallel > 1:
        logger.debug(f"Initializing pool of {n_tasks_parallel} processes")
        with Pool(n_tasks_parallel) as pool:
            logger.debug(f"Opened pool of {n_tasks_parallel} processes")
            report_files = tuple(pool.map(vectorize, writers))
        logger.debug(f"Closed pool of {n_tasks_parallel} processes")
    else:
        report_files = tuple(map(vectorize, writers))
    # Filter out any None values (indicating failure), convert report
    # paths to a tuple of strings, and return.
    reports = tuple(map(str, filter(None, report_files)))
    logger.info(f"Ended generating mutational profiles: {len(reports)} pass, "
                f"{len(report_files) - len(reports)} fail")
    return reports
