from enum import Enum
import itertools
import logging
import re
from functools import cached_property
from typing import BinaryIO, Dict, Iterable

from ..util import cli, docdef, path
from ..util.excmd import (run_cmd, BOWTIE2_CMD, BOWTIE2_BUILD_CMD,
                          CUTADAPT_CMD, FASTQC_CMD, SAMTOOLS_CMD)
from ..util.seq import FastaParser

# SAM file format specifications
SAM_HEADER = b"@"
SAM_ALIGN_SCORE = b"AS:i:"
SAM_EXTRA_SCORE = b"XS:i:"

# Bowtie2 parameters
MATCH_BONUS = "1"
MISMATCH_PENALTY = "1,1"
N_PENALTY = "0"
REF_GAP_PENALTY = "0,1"
READ_GAP_PENALTY = "0,1"


def index_bam_file(bam: path.BasePath):
    bam_index = bam.edit(ext=path.BAI_EXT)
    cmd = [SAMTOOLS_CMD, "index", bam]
    run_cmd(cmd)
    if not bam_index.path.is_file():
        logging.critical(f"Failed to index file: {bam}")


class FastqUnit(object):
    """
    Unified interface for handling the following sets of sequencing reads:

    - One FASTQ file of single-end reads from one sample
    - One FASTQ file of interleaved, paired-end reads from one sample
    - Two FASTQ files of mate 1 and mate 2 paired-end reads from one sample
    - One FASTQ file of single-end reads originating from one reference sequence
      in one sample
    - One FASTQ file of interleaved, paired-end reads originating from one
      reference sequence in one sample
    - Two FASTQ files of mate 1 and mate 2 paired-end reads originating from one
      reference sequence in one sample

    """

    MAX_PHRED_ENC = 127  # 2^7 - 1

    class KeyName(Enum):
        SINGLE = "fastqs"  # single-end reads
        INTER = "fastqi"  # interleaved paired-end reads
        MATE1 = "fastq1"  # mate 1 paired-end reads
        MATE2 = "fastq2"  # mate 2 paired-end reads

    class Bowtie2Flag(Enum):
        SINGLE = "-U"
        INTER = "--interleaved"
        MATE1 = "-1"
        MATE2 = "-2"

    class KeySampleType(Enum):
        SINGLE = path.AbstractSampleReadsFilePath
        INTER = path.AbstractSampleReadsFilePath
        MATE1 = path.AbstractSampleReads1FilePath
        MATE2 = path.AbstractSampleReads2FilePath

    class KeySampleInType(Enum):
        SINGLE = path.SampleReadsInFilePath
        INTER = path.SampleReadsInFilePath
        MATE1 = path.SampleReads1InFilePath
        MATE2 = path.SampleReads2InFilePath

    class KeyOneRefType(Enum):
        SINGLE = path.AbstractOneRefReadsFilePath
        INTER = path.AbstractOneRefReadsFilePath
        MATE1 = path.AbstractOneRefReads1FilePath
        MATE2 = path.AbstractOneRefReads2FilePath

    class KeySampleDirName(Enum):
        SINGLE = "fastqs_dir"  # single-end reads
        INTER = "fastqi_dir"  # interleaved paired-end reads
        MATE12 = "fastq12_dir"  # mate 1 and mate 2 paired-end reads

    class KeySampleDirType(Enum):
        SINGLE = path.AbstractSampleDirPath
        INTER = path.AbstractSampleDirPath
        MATE12 = path.AbstractSampleDirPath

    def __init__(self, /, *,
                 fastqs: (path.AbstractSampleReadsFilePath |
                          path.AbstractOneRefReadsFilePath |
                          None) = None,
                 fastqi: (path.AbstractSampleReadsFilePath |
                          path.AbstractOneRefReadsFilePath |
                          None) = None,
                 fastq1: (path.AbstractSampleReads1FilePath |
                          path.AbstractOneRefReads1FilePath |
                          None) = None,
                 fastq2: (path.AbstractSampleReads2FilePath |
                          path.AbstractOneRefReads2FilePath |
                          None) = None,
                 phred_enc: int):
        self._phred_enc = phred_enc
        self._init_inputs = dict(zip(self._get_keyname_strs(),
                                     (fastqs, fastqi, fastq1, fastq2),
                                     strict=True))
        self._validate_inputs()

    def _validate_inputs(self):
        _ = self.phred_arg
        _ = self.inputs
        _ = self.sample
        _ = self.ref

    def dict(self):
        return dict(sample=self.sample, ref=self.ref)

    @property
    def phred_enc(self):
        if 0 <= self._phred_enc <= self.MAX_PHRED_ENC:
            return self._phred_enc
        raise ValueError(self._phred_enc)

    @property
    def phred_arg(self):
        return f"--phred{self.phred_enc}"

    @classmethod
    def _get_keyname_strs(cls):
        return tuple(str(key.value) for key in cls.KeyName)

    @classmethod
    def _test_single(cls, keys: tuple[str, ...]):
        """ Return whether the arguments match a single-end FASTQ.
        No validation is performed. """
        return keys == (cls.KeyName.SINGLE.value,)

    @classmethod
    def _test_interleaved(cls, keys: tuple[str, ...]):
        """ Return whether the arguments match an interleaved FASTQ.
        No validation is performed. """
        return keys == (cls.KeyName.INTER.value,)

    @classmethod
    def _test_separate_mates(cls, keys: tuple[str, ...]):
        """ Return whether the arguments match two separate paired-end
        FASTQs. No validation is performed. """
        return keys == (cls.KeyName.MATE1.value, cls.KeyName.MATE2.value)

    @cached_property
    def _input_keys(self):
        """ Validate that a correct set of keyword arguments was given, and
        then return them. """
        # Find all the keys of self._inputs that got a value (not None).
        keys = tuple(key for key in self._get_keyname_strs()
                     if self._init_inputs[key] is not None)
        # The keys must match one of the following patterns to be valid.
        if not any((self._test_single(keys),
                    self._test_interleaved(keys),
                    self._test_separate_mates(keys))):
            raise ValueError(f"Got invalid FASTQ keywords: {keys}")
        # Now the keys are validated, and can validate other properties.
        return keys

    @cached_property
    def _inputs(self) -> Dict[str, (path.AbstractSampleReadsFilePath |
                                    path.AbstractOneRefReadsFilePath |
                                    path.AbstractSampleReads1FilePath |
                                    path.AbstractOneRefReads1FilePath |
                                    path.AbstractSampleReads2FilePath |
                                    path.AbstractOneRefReads2FilePath)]:
        """ Return a dictionary of keyword arguments and path instances of the
        input FASTQ files. Keys are validated, but values and types are not. """
        return {key: self._init_inputs[key] for key in self._input_keys}

    @property
    def single(self):
        """ Return whether the arguments match a single-end FASTQ.
        Validated. """
        return self._test_single(self._input_keys)

    @property
    def interleaved(self):
        """ Return whether the arguments match an interleaved FASTQ.
        Validated. """
        return self._test_interleaved(self._input_keys)

    @property
    def separate_mates(self):
        """ Return whether the arguments match two separate paired-end FASTQs.
        Validated. """
        return self._test_separate_mates(self._input_keys)

    @property
    def paired(self):
        """ Return whether the reads in the FASTQ unit are paired-end.
        Validated. """
        return not self.single

    @classmethod
    def _test_by_sample(cls, keys: tuple[str, ...], inputs: Iterable):
        """ Return whether the given path types match the path types expected
        for FASTQ files representing an entire sample. """
        return all(isinstance(i, cls.KeySampleType[cls.KeyName(key).name].value)
                   for key, i in zip(keys, inputs, strict=True))

    @classmethod
    def _test_by_ref(cls, keys: tuple[str, ...], inputs: Iterable):
        """ Return whether the given path types match the path types expected
        for FASTQ files representing one reference from a sample. """
        return all(isinstance(i, cls.KeyOneRefType[cls.KeyName(key).name].value)
                   for key, i in zip(keys, inputs, strict=True))

    @classmethod
    def _by_sample_or_ref(cls, keys: tuple[str, ...], inputs: Iterable):
        """ Validate that the types of paths given as input arguments
        match the expected types. """
        # The types must match one of the expected set of path types
        # given the keyword arguments -- either the sample- or ref-
        # based paths.
        by_sample = cls._test_by_sample(keys, inputs)
        by_ref = cls._test_by_ref(keys, inputs)
        if by_sample == by_ref:
            raise TypeError(f"Got invalid input FASTQ types: {inputs}")
        # Now the types are validated and can validate other properties.
        return by_sample, by_ref

    @cached_property
    def by_ref(self):
        _, by_ref = self._by_sample_or_ref(self._input_keys,
                                           tuple(self._inputs.values()))
        return by_ref

    @property
    def by_sample(self):
        return not self.by_ref

    @property
    def inputs(self):
        """ Return a dictionary of keyword arguments and path instances
        of the input FASTQ files. Validated. """
        # Confirm that the inputs either match the expected types for
        # sample- or ref-based inputs.
        _ = self.by_ref
        # Now the inputs are validated.
        return self._inputs

    @property
    def paths(self):
        return tuple(inp.path for inp in self.inputs.values())

    @property
    def parent(self) -> path.BasePath:
        parents = [inp.parent for inp in self.inputs.values()]
        if not parents:
            raise TypeError("Not parent directory")
        if any(parent != parents[0] for parent in parents[1:]):
            raise ValueError("More than one parent directory")
        return parents[0]

    @property
    def str_paths(self):
        return tuple(map(str, self.inputs.values()))

    @cached_property
    def sample(self):
        samples: tuple[str] = tuple({fq.sample for fq in self.inputs.values()})
        if len(samples) != 1:
            raise ValueError(f"Required exactly one sample, but got {samples}")
        return samples[0]

    @cached_property
    def ref(self):
        if not self.by_ref:
            return None
        refs: tuple[str] = tuple({fq.ref for fq in self.inputs.values()})
        if len(refs) != 1:
            raise ValueError(f"Required exactly one reference, but got {refs}")
        return refs[0]

    def compatible_fasta(self, fasta: (path.AbstractRefsetSeqFilePath |
                                       path.OneRefSeqStepFilePath)):
        """ Return whether a given FASTA file is compatible with the
        FastqUnit. They are compatible if and only if
        - the FASTA has a set of references and the FASTQ has reads from
          an entire sample
        - the FASTA has one reference and the FASTQ has reads from only
          the same reference
        """
        if isinstance(fasta, path.AbstractRefsetSeqFilePath):
            return self.by_sample
        if isinstance(fasta, path.OneRefSeqStepFilePath):
            return self.ref == fasta.ref
        raise TypeError(fasta)

    @property
    def cutadapt_input_args(self):
        return self.paths

    @property
    def _bowtie2_flags(self):
        return tuple(str(self.Bowtie2Flag[self.KeyName(key).name].value)
                     for key in self._input_keys)

    @property
    def bowtie2_inputs(self):
        return tuple(itertools.chain(*map(list, zip(self._bowtie2_flags,
                                                    self.paths,
                                                    strict=True))))

    def edit(self, **fields):
        """
        Return a new FastqUnit by modifying all the FASTQ paths, using
        the given fields.

        Parameters
        ----------
        **fields: Any
            Keyword arguments to determine the new paths

        Returns
        -------
        FastqUnit
            New FastqUnit with the same keys and modified paths
        """
        return FastqUnit(**{key: fq.edit(**{**fields, path.EXT_KEY: fq.ext})
                            for key, fq in self.inputs.items()},
                         phred_enc=self.phred_enc)

    @classmethod
    def _get_demult_files(cls, phred_enc: int, dir_path: path.BasePath,
                          key: str):
        """
        Yield a FastqUnit for each demultiplexed FASTQ file in a
        directory.

        Parameters
        ----------
        phred_enc: int
            Phred score encoding
        dir_path: SampleInDirPath
            Directory containing the FASTQ files
        key: str
            Keyword argument for every FASTQ file ('fastqs' or 'fastqi')

        Yield
        -----
        FastqUnit
            One for each FASTQ file in the directory
        """
        dp = dir_path.path
        parse_type = cls.KeyOneRefType[cls.KeyName(key).name].value
        for fname in dp.iterdir():
            yield cls(**{key: parse_type.parse(dp.joinpath(fname))},
                      phred_enc=phred_enc)

    @classmethod
    def _get_demult_pairs(cls, phred_enc: int, dir_path: path.BasePath):
        """
        Yield a FastqUnit for each pair of mated, demultiplexed FASTQ files
        in a directory.

        Parameters
        ----------
        phred_enc: int
            Phred score encoding
        dir_path: SampleInDirPath
            Directory containing the FASTQ files

        Return
        ------
        Iterable[FastqUnit]
            One for each FASTQ pair in the directory
        """
        dp = dir_path.path
        # Create empty dictionaries to store the FASTQ files for mates 1 and 2,
        # keyed by the name of the reference sequence for each FASTQ file.
        mates1: dict[str, path.AbstractOneRefReads1FilePath] = dict()
        mates2: dict[str, path.AbstractOneRefReads2FilePath] = dict()
        # Read the name of every file in the directory.
        for fname in dp.iterdir():
            # Get the full path to the file.
            fpath = str(dp.joinpath(fname))
            # Determine if the file is a mate 1 or mate 2 FASTQ file.
            is1 = any(fpath.endswith(ext) for ext in path.FQ1_EXTS)
            is2 = any(fpath.endswith(ext) for ext in path.FQ2_EXTS)
            if is1 and is2:
                # There should be no way for this error to happen, but catching
                # it just in case.
                raise ValueError(f"FASTQ path matched both mates: '{fpath}'")
            if is1:
                # The file name matched an extension for a mate 1 FASTQ file.
                fq = path.AbstractOneRefReads1FilePath.parse(fpath)
                if fq.ref in mates1:
                    raise ValueError(f"Got >1 FASTQ file for ref '{fq.ref}'")
                # Add the path to the dict of mate 1 files, keyed by reference.
                mates1[fq.ref] = fq
            elif is2:
                # The file name matched an extension for a mate 2 FASTQ file.
                fq = path.AbstractOneRefReads2FilePath.parse(fpath)
                if fq.ref in mates2:
                    raise ValueError(f"Got >1 FASTQ file for ref '{fq.ref}'")
                # Add the path to the dict of mate 2 files, keyed by reference.
                mates2[fq.ref] = fq
            else:
                # If a file name does not match the expected FASTQ name format,
                # log a warning but keep going, since the presence or absence
                # of one FASTQ file will not affect the others, and this file
                # might just be an extraneous file, such as .DS_Store on macOS.
                logging.warning(f"File name is not a valid FASTQ: '{fpath}'")
        # Iterate through all references from mate 1 and mate 2 files.
        for ref in set(mates1) | set(mates2):
            fastq1 = mates1.get(ref)
            fastq2 = mates2.get(ref)
            if fastq1 is not None and fastq2 is not None:
                # Yield a new FastqUnit from the paired FASTQ files.
                yield cls(fastq1=fastq1, fastq2=fastq2, phred_enc=phred_enc)
            elif fastq1 is None:
                logging.error(f"Missing mate 1 for reference '{ref}'")
            else:
                logging.error(f"Missing mate 2 for reference '{ref}'")

    @classmethod
    def _get_sample_files(cls, phred_enc: int, path_strs: tuple[str], key: str):
        file_type = cls.KeySampleInType[cls.KeyName(key).name].value
        for path_str in path_strs:
            yield FastqUnit(**{key: file_type.parse(path_str)},
                            phred_enc=phred_enc)

    @classmethod
    def _get_sample_pairs(cls, phred_enc: int,
                          path1_strs: tuple[str],
                          path2_strs: tuple[str]):
        keys = (cls.KeyName.MATE1, cls.KeyName.MATE2)
        for path_strs in zip(path1_strs, path2_strs, strict=True):
            paths = {
                key.value: cls.KeySampleInType[key.name].value.parse(path_str)
                for key, path_str in zip(keys, path_strs, strict=True)
            }
            yield FastqUnit(**paths, phred_enc=phred_enc)

    @classmethod
    def _from_strs(cls, /, *, phred_enc: int, **fastq_args: tuple[str]):
        """ Yield a FastqUnit for each FASTQ file (or each pair of mated
        FASTQ files) whose paths are given as strings. """
        path1_strs = None
        for fq_key, path_strs in fastq_args.items():
            try:
                # Check if fq_key is a key for a directory or file.
                dir_type = cls.KeySampleDirType[
                    cls.KeySampleDirName(fq_key).name].value
            except ValueError:
                # fq_key is not a directory key: assume a file key.
                if (fq_key == cls.KeyName.MATE1.value
                        or fq_key == cls.KeyName.MATE2.value):
                    # fq_key is for files of mate 1 or mate 2 reads.
                    if path1_strs is None:
                        # Because both keys for mate 1 and mate 2 must
                        # be read together to generate the FASTQ units,
                        # confirm that both keys were given.
                        # To prevent yielding each of these FASTQ units
                        # twice (once when fq_key is mate 1 and again
                        # when fq_key is mate 2), this block runs only
                        # if path1_strs is None (which is no longer true
                        # after the first time this block runs).
                        try:
                            path1_strs = fastq_args[cls.KeyName.MATE1.value]
                            path2_strs = fastq_args[cls.KeyName.MATE2.value]
                        except KeyError:
                            raise KeyError(
                                f"Must give both '{cls.KeyName.MATE1.value}' "
                                f"and '{cls.KeyName.MATE2.value}', or neither")
                        #
                        yield from cls._get_sample_pairs(phred_enc=phred_enc,
                                                         path1_strs=path1_strs,
                                                         path2_strs=path2_strs)
                else:
                    # fq_key is for files of single-end or interleaved
                    # paired-end reads. Yield every file from the key.
                    yield from cls._get_sample_files(phred_enc=phred_enc,
                                                     path_strs=path_strs,
                                                     key=fq_key)
            else:
                # The FASTQ key is a key for a directory.
                # For each path string, parse the directory path.
                for dir_path in map(dir_type.parse, path_strs):
                    # Yield each FastqUnit depending on the FASTQ key.
                    if fq_key == cls.KeySampleDirName.MATE12.value:
                        # The directory contains paired reads with mate
                        # 1 and mate 2 reads in separate FASTQ files.
                        yield from cls._get_demult_pairs(phred_enc=phred_enc,
                                                         dir_path=dir_path)
                    else:
                        # The directory contains FASTQ files that may be
                        # processed separately (single-end/interleaved).
                        yield from cls._get_demult_files(phred_enc=phred_enc,
                                                         dir_path=dir_path,
                                                         key=fq_key)

    @classmethod
    def from_strs(cls, /, *,
                  phred_enc: int,
                  no_dup_samples: bool = True,
                  **fastq_args: tuple[str]):
        """
        Yield a FastqUnit for each FASTQ file (or each pair of mate 1
        and mate 2 FASTQ files) whose paths are given as strings.

        Parameters
        ----------
        phred_enc: int
            ASCII offset for encoding Phred scores
        no_dup_samples: bool (default: True)
            If a sample name occurs more than once among the given FASTQ
            files, then log a warning and yield only the first FASTQ
            with that sample name.
        fastq_args: tuple[str]
            FASTQ files, given as tuples of file path string:
            - fastqs: FASTQ files of single-end reads
            - fastqi: FASTQ files of interleaved paired-end reads
            - fastq1: FASTQ files of mate 1 paired-end reads; must
                      correspond 1-for-1 (in order) with fastq2
            - fastq2: FASTQ files of mate 2 paired-end reads; must
                      correspond 1-for-1 (in order) with fastq1
            - fastqs_dir: Directory of FASTQ files of single-end reads
            - fastqi_dir: Directory of FASTQ files of interleaved
                          paired-end reads
            - fastq12_dir: Directory of FASTQ files of separate mate 1
                           and mate 2 paired-end reads; for every FASTQ
                           file of mate 1 reads, there must be a FASTQ
                           file of mate 2 reads with the same sample
                           name, and vice versa.

        Yield
        -----
        FastqUnit
            FastqUnit representing the FASTQ or pair of FASTQ files.
            The order is determined primarily by the order of keyword
            arguments; within each keyword argument, by the order of
            file or directory paths; and for directories, by the order
            in which ```pathlib.Path.iterdir``` returns file paths.
        """
        samples = set()
        for fq_unit in cls._from_strs(phred_enc=phred_enc, **fastq_args):
            if no_dup_samples:
                if (sample := fq_unit.sample) in samples:
                    logging.warning(f"Skipping duplicate sample: '{sample}'")
                    continue
                samples.add(sample)
            yield fq_unit

    def __str__(self):
        return " and ".join(self.str_paths)


class ReadsFileBase(object):
    _module: path.Module | None = None
    _step: str | None = None
    _ext: str | None = None

    def __init__(self, /, *,
                 input_path: (FastqUnit |
                              path.AbstractRefsetAlignmentFilePath |
                              path.AbstractOneRefAlignmentFilePath),
                 top_dir: str,
                 n_procs: int,
                 save_temp: bool,
                 resume: bool):
        self._input = input_path
        self._top_dir = path.TopDirPath.parse(top_dir)
        self._n_procs = n_procs
        self._save_temp = save_temp
        self._resume = resume
        self._output_files: list[path.BasePath] = list()
        self._index_files: list[path.BasePath] = list()

    @property
    def output_fields(self):
        return {**self._input.dict(),
                **dict(top=self._top_dir.top,
                       module=self._module,
                       step=self._step,
                       ext=self._ext)}

    @property
    def output(self):
        return self._input.edit(**self.output_fields)

    def _setup(self):
        self.output.parent.path.mkdir(parents=True, exist_ok=True)

    def _run(self, **kwargs):
        """ Run this step. Each class implements a unique function. """
        logging.error("The base class does not implement a run method")

    def _clean_index_files(self):
        while self._index_files:
            self._index_files.pop().path.unlink(missing_ok=True)

    def _needs_to_run(self):
        return not self._resume or not self.output.path.is_file()

    def run(self, **kwargs):
        """ Wrapper for internal run method that first checks if the
        step actually needs to be run because its output files might
        already exist. """
        if self._needs_to_run():
            # If the output of this step does not exist, or it does
            # (because, e.g. the alignment pipeline crashed after at
            # least one step had written a temporary output file), but
            # the pipeline is not allowed to resume from the leftover
            # temporary files, then run this step.
            self._setup()
            self._run(**kwargs)
            if not self._save_temp:
                self._clean_index_files()
        # Return the output of this step so that it can be fed into the
        # next step.
        return self.output

    def clean(self):
        if self._step is not None and not self._save_temp:
            while self._output_files:
                self._output_files.pop().path.unlink(missing_ok=True)


class FastqBase(ReadsFileBase):
    _module = path.Module.ALIGN

    def __init__(self, /, *, fq_unit: FastqUnit, **kwargs):
        super().__init__(**kwargs, input_path=fq_unit)

    @property
    def sample(self):
        return self._input.sample

    @staticmethod
    def _qc(fastq_paths: list[path.BasePath], extract: bool):
        cmd = [FASTQC_CMD]
        if extract:
            cmd.append("--extract")
        cmd.extend(fastq_paths)
        run_cmd(cmd)

    def qc(self, extract: bool):
        self._qc(list(self._input.inputs.values()), extract)


class FastqTrimmer(FastqBase):
    _step = path.Step.ALIGN_TRIM
    _cutadapt_output_flags = "-o", "-p"

    @cached_property
    def output(self) -> FastqUnit:
        return self._input.edit(**self.output_fields)

    @property
    def _cutadapt_output_args(self):
        """ Return arguments that specify Cutadapt output files. """
        return tuple(itertools.chain(*zip(self._cutadapt_output_flags,
                                          self.output.paths,
                                          strict=False)))

    @docdef.auto()
    def _cutadapt(self, /, *,
                  cut_q1: int,
                  cut_q2: int,
                  cut_g1: str,
                  cut_a1: str,
                  cut_g2: str,
                  cut_a2: str,
                  cut_o: int,
                  cut_e: float,
                  cut_indels: bool,
                  cut_nextseq: bool,
                  cut_discard_trimmed: bool,
                  cut_discard_untrimmed: bool,
                  cut_m: int):
        """ Trim adapters and low-quality bases with Cutadapt. """
        cmd = [CUTADAPT_CMD]
        cmd.extend(["--cores", self._n_procs])
        if cut_nextseq:
            if cut_q1 > 0:
                cmd.extend(["--nextseq-trim", cut_q1])
        else:
            if cut_q1 > 0:
                cmd.extend(["-q", cut_q1])
            if cut_q2 > 0:
                cmd.extend(["-Q", cut_q2])
        adapters = {"g": cut_g1, "a": cut_a1,
                    "G": cut_g2, "A": cut_a2}
        for arg, adapter in adapters.items():
            if adapter and (self._input.paired or arg.islower()):
                for adapt in adapter:
                    cmd.extend([f"-{arg}", adapt])
        cmd.extend(["-O", cut_o])
        cmd.extend(["-e", cut_e])
        cmd.extend(["-m", cut_m])
        if not cut_indels:
            cmd.append("--no-indels")
        if cut_discard_trimmed:
            cmd.append("--discard-trimmed")
        if cut_discard_untrimmed:
            cmd.append("--discard-untrimmed")
        if self._input.interleaved:
            cmd.append("--interleaved")
        cmd.extend(["--report", "minimal"])
        cmd.extend(self._cutadapt_output_args)
        cmd.extend(self._input.cutadapt_input_args)
        run_cmd(cmd)
        self._output_files.extend(self.output.inputs.values())

    def _needs_to_run(self):
        return not self._resume or not all(output.is_file()
                                           for output in self.output.paths)

    def _run(self, **kwargs):
        self._cutadapt(**kwargs)


class FastqAligner(FastqBase):
    _step = path.Step.ALIGN_ALIGN
    _ext = path.SAM_EXT

    def __init__(self, /, *,
                 fasta: (path.AbstractRefsetSeqFilePath |
                         path.OneRefSeqStepFilePath),
                 **kwargs):
        super().__init__(**kwargs)
        if not self._input.compatible_fasta(fasta):
            raise TypeError(f"FASTA {fasta} is not compatible with "
                            f"FASTQ file(s) {self._input}")
        self._fasta = fasta

    @cached_property
    def output(self):
        return path.create(**{**self._input.dict(),
                              **self._fasta.dict(),
                              **self.output_fields})

    @cached_property
    def _bowtie2_index_files(self):
        return [self._fasta.edit(ext=ext)
                for ext in path.BOWTIE2_INDEX_EXTS]

    @cached_property
    def _fasta_prefix(self):
        prefix = str(self._fasta.path.with_suffix(""))
        if bad_idxs := [index for index in map(str, self._bowtie2_index_files)
                        if not index.startswith(prefix)]:
            logging.critical(f"Bowtie2 index files {bad_idxs} do not start "
                             f"with FASTA prefix '{prefix}'")
        return prefix

    @property
    def _missing_bowtie2_index_files(self):
        return [index for index in self._bowtie2_index_files
                if not index.path.is_file()]

    def _bowtie2_build(self):
        """ Build an index of a reference using Bowtie 2. """
        cmd = [BOWTIE2_BUILD_CMD, self._fasta, self._fasta_prefix]
        run_cmd(cmd)
        self._index_files.extend(self._bowtie2_index_files)

    @docdef.auto()
    def _bowtie2(self, /, *,
                 bt2_local: bool,
                 bt2_discordant: bool,
                 bt2_mixed: bool,
                 bt2_dovetail: bool,
                 bt2_contain: bool,
                 bt2_unal: bool,
                 bt2_score_min: str,
                 bt2_i: int,
                 bt2_x: int,
                 bt2_gbar: int,
                 bt2_l: int,
                 bt2_s: str,
                 bt2_d: int,
                 bt2_r: int,
                 bt2_dpad: int,
                 bt2_orient: str):
        """ Align reads to the reference with Bowtie 2. """
        if missing := self._missing_bowtie2_index_files:
            logging.critical("Bowtie2 index files do not exist:\n\n"
                             + "\n".join(map(str, missing)))
            return
        cmd = [BOWTIE2_CMD]
        # Resources
        cmd.extend(["--threads", self._n_procs])
        # Alignment
        if bt2_local:
            cmd.append("--local")
        cmd.extend(["--gbar", bt2_gbar])
        cmd.extend(["--dpad", bt2_dpad])
        cmd.extend(["-L", bt2_l])
        cmd.extend(["-i", bt2_s])
        cmd.extend(["-D", bt2_d])
        cmd.extend(["-R", bt2_r])
        # Scoring
        cmd.append(self._input.phred_arg)
        cmd.append("--ignore-quals")
        cmd.extend(["--ma", MATCH_BONUS])
        cmd.extend(["--mp", MISMATCH_PENALTY])
        cmd.extend(["--np", N_PENALTY])
        cmd.extend(["--rfg", REF_GAP_PENALTY])
        cmd.extend(["--rdg", READ_GAP_PENALTY])
        # Filtering
        if not bt2_unal:
            cmd.append("--no-unal")
        cmd.extend(["--score-min", bt2_score_min])
        cmd.extend(["-I", bt2_i])
        cmd.extend(["-X", bt2_x])
        # Mate pair orientation
        orientations = tuple(op.value for op in cli.MateOrientationOption)
        if bt2_orient in orientations:
            cmd.append(f"--{bt2_orient}")
        else:
            cmd.append(f"--{orientations[0]}")
            logging.warning(f"Invalid mate orientation: '{bt2_orient}'; "
                            f"defaulting to '{orientations[0]}'")
        if not bt2_discordant:
            cmd.append("--no-discordant")
        if not bt2_contain:
            cmd.append("--no-contain")
        if bt2_dovetail:
            cmd.append("--dovetail")
        if not bt2_mixed:
            cmd.append("--no-mixed")
        # Formatting
        cmd.append("--xeq")
        # Input and output files
        cmd.extend(["-S", self.output.path])
        cmd.extend(["-x", self._fasta_prefix])
        cmd.extend(self._input.bowtie2_inputs)
        # Run alignment.
        run_cmd(cmd)
        self._output_files.append(self.output)

    def _run(self, **kwargs):
        if self._missing_bowtie2_index_files:
            self._bowtie2_build()
        self._bowtie2(**kwargs)


class XamBase(ReadsFileBase):
    def view(self, output: (path.AbstractRefsetAlignmentFilePath |
                            path.AbstractOneRefAlignmentFilePath |
                            path.AbstractRegionAlignmentFilePath |
                            path.BasePath)):
        cmd = [SAMTOOLS_CMD, "view"]
        if self._ext == path.BAM_EXT:
            cmd.append("-b")
        cmd.extend(("-o", output, self._input))
        if isinstance(output, path.AbstractRefsetAlignmentFilePath):
            # View all reads. No index is required.
            pass
        elif isinstance(output, path.AbstractOneRefAlignmentFilePath):
            # View only reads that aligned to the reference.
            cmd.append(output.ref)
            # This restriction requires an index.
            self._build_index(self._input)
        elif isinstance(output, (path.RegionAlignmentInFilePath |
                                 path.RegionAlignmentStepFilePath |
                                 path.RegionAlignmentOutFilePath)):
            # View only reads that aligned to a region of the reference.
            cmd.append(f"{output.ref}:{output.end5}-{output.end3}")
            # This restriction requires an index.
            self._build_index(self._input)
        else:
            raise TypeError(output)
        run_cmd(cmd)
        if output.path.is_file():
            self._output_files.append(output)
        else:
            logging.critical(f"Failed to view file: {self._input} -> {output}")

    @staticmethod
    def _get_index_path(bam: path.BasePath):
        return bam.edit(ext=path.BAI_EXT)

    def _build_index(self, bam: path.BasePath):
        bam_index = self._get_index_path(bam)
        if not self._resume or not bam_index.path.is_file():
            index_bam_file(bam)
            self._index_files.append(bam_index)


class SamRemoveEqualMappers(XamBase):
    _module = path.Module.ALIGN
    _step = path.Step.ALIGN_REMEQ
    _ext = path.SAM_EXT

    pattern_a = re.compile(SAM_ALIGN_SCORE + rb"(\d+)")
    pattern_x = re.compile(SAM_EXTRA_SCORE + rb"(\d+)")

    _MIN_SAM_FIELDS = 11
    _MAX_SAM_FLAG = 4095  # 2^12 - 1

    @staticmethod
    def _get_score(line: bytes, ptn: re.Pattern[bytes]):
        return (float(match.groups()[0])
                if (match := ptn.search(line)) else None)

    @classmethod
    def _is_best_alignment(cls, line: bytes):
        return ((score_x := cls._get_score(line, cls.pattern_x)) is None
                or score_x < cls._get_score(line, cls.pattern_a))

    @classmethod
    def _read_is_paired(cls, line: bytes):
        info = line.split()
        if len(info) < cls._MIN_SAM_FIELDS:
            raise ValueError(f"Invalid SAM line:\n{line.decode()}")
        flag = int(info[1])
        if 0 <= flag <= cls._MAX_SAM_FLAG:
            return bool(flag % 2)
        raise ValueError(f"Invalid SAM flag: {flag}")

    def _iter_paired(self, sam: BinaryIO, line: bytes):
        """ For each pair of reads, yield the pair of alignments for
        which both the forward alignment and the reverse alignment in
        the pair scored best among all alignments for the forward and
        reverse reads, respectively. Exclude pairs for which the forward
        and/or reverse read aligned equally well to multiple locations,
        or for which the best alignments for the forward and reverse
        reads individually are not part of the same alignment pair. """
        total_lines = 0
        yield_lines = 0
        while line:
            total_lines += 1
            line2 = sam.readline()
            total_lines += bool(line2)
            if self._is_best_alignment(line) and self._is_best_alignment(line2):
                yield_lines += 1
                yield line
                yield_lines += 1
                yield line2
            line = sam.readline()
        if total_lines % 2:
            logging.error(f"SAM file {self._input} was paired but had an odd "
                          f"number of lines ({total_lines})")

    def _iter_single(self, sam: BinaryIO, line: bytes):
        """ For each read, yield the best-scoring alignment, excluding
        reads that aligned equally well to multiple locations. """
        while line:
            if self._is_best_alignment(line):
                yield line
            line = sam.readline()

    def _remove_equal_mappers(self, buffer_length):
        with (open(self._input.path, "rb") as sami,
              open(self.output.path, "wb") as samo):
            # Copy the header from the input to the output SAM file.
            while (line := sami.readline()).startswith(SAM_HEADER):
                samo.write(line)
            if line:
                if self._read_is_paired(line):
                    lines = self._iter_paired(sami, line)
                else:
                    lines = self._iter_single(sami, line)
                while text := b"".join(itertools.islice(lines, buffer_length)):
                    samo.write(text)

    def _run(self, rem_buffer: int = cli.opt_rem_buffer.default):
        logging.info("\nRemoving Reads Mapping Equally to Multiple Locations"
                     f" in {self._input.path}\n")
        self._remove_equal_mappers(rem_buffer)
        self._output_files.append(self.output)


class XamSorter(XamBase):
    def _sort(self, name: bool):
        cmd = [SAMTOOLS_CMD, "sort"]
        if name:
            cmd.append("-n")
        cmd.extend(["-o", self.output.path, self._input.path])
        run_cmd(cmd)
        self._output_files.append(self.output)

    def _run(self, name: bool = False):
        logging.info(f"\nSorting {self._input.path} by Reference and Coordinate\n")
        self._sort(name)


class BamAlignSorter(XamSorter):
    _module = path.Module.ALIGN
    _step = path.Step.ALIGN_SORT
    _ext = path.BAM_EXT


class SamVectorSorter(XamSorter):
    _module = path.Module.VECTOR
    _step = path.Step.VECTOR_SORT
    _ext = path.SAM_EXT


class BamSplitter(XamBase):
    _module = path.Module.ALIGN
    _step = path.Step.ALIGN_SPLIT
    _ext = path.BAM_EXT

    def __init__(self, /, *,
                 fasta: path.RefsetSeqInFilePath | path.OneRefSeqStepFilePath,
                 **kwargs):
        super().__init__(**kwargs)
        self.fasta = fasta

    @cached_property
    def refs(self):
        return tuple(ref for ref, _ in FastaParser(self.fasta.path).parse())

    @cached_property
    def output(self):
        output = [self._input.edit(**{**self.output_fields,
                                      **dict(refset=None, ref=ref)})
                  for ref in self.refs]
        return output

    def _split(self):
        for output in self.output:
            self.view(output)

    def _needs_to_run(self):
        return not self._resume or not all(output.path.is_file()
                                           for output in self.output)

    def _setup(self):
        for output in self.output:
            output.parent.path.mkdir(parents=True, exist_ok=True)

    def _run(self):
        logging.info(f"\nSplitting {self._input} by reference\n")
        self._build_index(self._input)
        self._split()


class BamOutputter(XamBase):
    _module = path.Module.ALIGN
    _ext = path.BAM_EXT

    def _run(self):
        logging.info(f"\nOutputting {self._input} to {self.output}\n")
        self.view(self.output)


class BamVectorSelector(XamBase):
    _module = path.Module.VECTOR
    _step = path.Step.VECTOR_SELECT
    _ext = path.BAM_EXT

    def __init__(self, /, *, ref: str, end5: int, end3: int, **kwargs):
        super().__init__(**kwargs)
        self.ref = ref
        self.end5 = end5
        self.end3 = end3

    @staticmethod
    def ref_coords(ref: str, end5: int, end3: int):
        return f"{ref}:{end5}-{end3}"

    @property
    def fields(self):
        return dict(**super().output_fields,
                    end5=self.end5,
                    end3=self.end3)

    def _run(self):
        cmd = [SAMTOOLS_CMD, "view", "-h", "-o", self.output,
               self._input, self.ref_coords(self.ref, self.end5, self.end3)]
        run_cmd(cmd)


'''
primer1 = "CAGCACTCAGAGCTAATACGACTCACTATA"
primer1rc = "TATAGTGAGTCGTATTAGCTCTGAGTGCTG"
primer2 = "TGAAGAGCTGGAACGCTTCACTGA"
primer2rc = "TCAGTGAAGCGTTCCAGCTCTTCA"
adapters5 = (primer1, primer2rc)
adapters3 = (primer2, primer1rc)
'''
