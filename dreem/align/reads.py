from abc import ABC, abstractmethod
from enum import Enum
import itertools
import logging
import re
from functools import cached_property
from typing import BinaryIO, Iterable

from ..util import path
from ..util.cli import MateOrientationOption
from ..util.dflt import BUFFER_LENGTH
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
        SINGLE = path.SampleReadsInFilePath
        INTER = path.SampleReadsInFilePath
        MATE1 = path.SampleReads1InFilePath
        MATE2 = path.SampleReads2InFilePath

    class KeyDemultType(Enum):
        SINGLE = path.DemultReadsInFilePath
        INTER = path.DemultReadsInFilePath
        MATE1 = path.DemultReads1InFilePath
        MATE2 = path.DemultReads2InFilePath

    class KeySampleDirName(Enum):
        SINGLE = "fastqs_dir"  # single-end reads
        INTER = "fastqi_dir"  # interleaved paired-end reads
        MATE12 = "fastq12_dir"  # mate 1 and mate 2 paired-end reads

    class KeySampleDirType(Enum):
        SINGLE = path.SampleInDirPath
        INTER = path.SampleInDirPath
        MATE12 = path.SampleInDirPath

    def __init__(self, *,
                 fastqs: (path.SampleReadsInFilePath |
                          path.DemultReadsInFilePath |
                          None) = None,
                 fastqi: (path.SampleReadsInFilePath |
                          path.DemultReadsInFilePath |
                          None) = None,
                 fastq1: (path.SampleReads1InFilePath |
                          path.DemultReads1InFilePath |
                          None) = None,
                 fastq2: (path.SampleReads2InFilePath |
                          path.DemultReads2InFilePath |
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
        if self.demult:
            _ = self.ref

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

    @property
    def _inputs(self):
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
    def _test_by_sample(cls, keys: tuple[str, ...], types: tuple[type, ...]):
        """ Return whether the given path types match the path types expected
        for FASTQ files representing an entire sample. """
        return types == tuple(cls.KeySampleType[cls.KeyName(key).name].value
                              for key in keys)

    @classmethod
    def _test_by_ref(cls, keys: tuple[str, ...], types: tuple[type, ...]):
        """ Return whether the given path types match the path types expected
        for FASTQ files representing one reference from a sample. """
        return types == tuple(cls.KeyDemultType[cls.KeyName(key).name].value
                              for key in keys)

    @cached_property
    def _input_types(self):
        """ Validate that the types of paths given as input arguments match the
        expected types. """
        # Get the types of the input arguments.
        types = tuple(map(type, self._inputs.values()))
        # The types must match one of the expected set of path types given the
        # keyword arguments -- either the sample- or ref-based paths.
        if (self._test_by_sample(self._input_keys, types)
                == self._test_by_ref(self._input_keys, types)):
            raise TypeError(f"Got invalid FASTQ types: {types}")
        # Now the types are validated and can validate other properties.
        return types

    @cached_property
    def inputs(self):
        """ Return a dictionary of keyword arguments and path instances of the
        input FASTQ files. Validated. """
        inputs: dict[str, (path.SampleReadsInFilePath |
                           path.SampleReads1InFilePath |
                           path.SampleReads2InFilePath |
                           path.DemultReadsInFilePath |
                           path.DemultReads1InFilePath |
                           path.DemultReads2InFilePath)] = dict()
        for (key, path_inst), exp_type in zip(self._inputs.items(),
                                              self._input_types,
                                              strict=True):
            # Confirm that every item of self._inputs has the correct type.
            if (path_type := type(path_inst)) is not exp_type:
                raise TypeError(
                    f"Expected {exp_type} for key '{key}', but got {path_type}")
            inputs[key] = path_inst
        # Now the inputs are validated.
        return inputs

    @property
    def demult(self):
        return self._test_by_ref(self._input_keys, self._input_types)

    @property
    def paths(self):
        return tuple(inp.path for inp in self.inputs.values())

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
        if not self.demult:
            raise TypeError("Whole-sample FASTQ files do not define reference")
        refs: tuple[str] = tuple({fq.ref for fq in self.inputs.values()})
        if len(refs) != 1:
            raise ValueError(f"Required exactly one reference, but got {refs}")
        return refs[0]

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

    def trans(self, trans: type[path.PathTypeTranslator], **kwargs):
        """
        Return a new FastqUnit by translating all the FASTQ paths using
        the given PathTypeTranslator and optionally extra arguments.

        Parameters
        ----------
        trans: PathTypeTranslator
            Translation from one FASTQ type to another
        **kwargs: Any
            Keyword arguments passed to the ```trans_inst``` method of
            the path type translator

        Returns
        -------
        FastqUnit
            New FastqUnit with the same keys and translated paths
        """
        return FastqUnit(**{key: trans.trans_inst(inp, **kwargs)
                            for key, inp in self.inputs.items()},
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
        parse_type = cls.KeyDemultType[cls.KeyName(key).name].value
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
        mates1: dict[str, path.DemultReads1InFilePath] = dict()
        mates2: dict[str, path.DemultReads2InFilePath] = dict()
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
                fq = path.DemultReads1InFilePath.parse(fpath)
                if fq.ref in mates1:
                    raise ValueError(f"Got >1 FASTQ file for ref '{fq.ref}'")
                # Add the path to the dict of mate 1 files, keyed by reference.
                mates1[fq.ref] = fq
            elif is2:
                # The file name matched an extension for a mate 2 FASTQ file.
                fq = path.DemultReads2InFilePath.parse(fpath)
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
        # Determine if any reference names are present for either mate 1 or
        # mate 2 but not both mates.
        refs1 = set(mates1)
        refs2 = set(mates2)
        # Find union of all reference names from both mate 1 and mate 2 files.
        refs = refs1 | refs2
        # Find any reference names that exist for one mate but not both.
        if missing := refs - (refs1 & refs2):
            # If any reference names are missing for either mate, this is an
            # error because the FASTQ files should be generated in pairs, and
            # if not that indicates a problem, such as the accidental deletion
            # of one or more FASTQ files in the directory or the directory
            # actually containing single-end or interleaved FASTQ files.
            raise ValueError(f"Missing mate 1/2 for refs {', '.join(missing)}")
        for ref in refs:
            # Yield a FastqUnit for each pair of mated FASTQ files.
            yield cls(fastq1=mates1[ref], fastq2=mates2[ref],
                      phred_enc=phred_enc)

    @classmethod
    def _get_sample_files(cls, phred_enc: int, path_strs: tuple[str], key: str):
        file_type = cls.KeySampleType[cls.KeyName(key).name].value
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
                key.value: cls.KeySampleType[key.name].value.parse(path_str)
                for key, path_str in zip(keys, path_strs, strict=True)
            }
            yield FastqUnit(**paths, phred_enc=phred_enc)

    @classmethod
    def _from_strs(cls, *, phred_enc: int, **fastq_args: tuple[str]):
        """
        Yield a FastqUnit for each FASTQ file (or each pair of mate 1
        and mate 2 FASTQ files) whose paths are given as strings.
        """
        path1_strs = ()
        for fq_key, path_strs in fastq_args.items():
            try:
                # Check if fq_key is a key for a directory or file.
                dir_type = cls.KeySampleDirType[
                    cls.KeySampleDirName(fq_key).name].value
            except ValueError:
                # fq_key is not a directory key: assume a file key.
                if (fq_key == cls.KeyName.MATE1.value
                        or fq_key == cls.KeyName.MATE2.value):
                    if not path1_strs:
                        try:
                            path1_strs = fastq_args[cls.KeyName.MATE1.value]
                            path2_strs = fastq_args[cls.KeyName.MATE2.value]
                        except KeyError:
                            raise KeyError(
                                f"Must give both {cls.KeyName.MATE1.value} and "
                                f"{cls.KeyName.MATE2.value} or neither")
                        yield from cls._get_sample_pairs(phred_enc=phred_enc,
                                                         path1_strs=path1_strs,
                                                         path2_strs=path2_strs)
                else:
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
    def from_strs(cls, *, phred_enc: int, no_dup_samples: bool = True,
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
            FASTQ files, given as tuples of file path strings. At least
            one of the following keywords must be given:
            * fastqs: FASTQ files of single-end reads
            * fastqi: FASTQ files of interleaved paired-end reads
            * fastq1: FASTQ files of mate 1 paired-end reads; must
                      correspond 1-for-1 (in order) with fastq2
            * fastq2: FASTQ files of mate 2 paired-end reads; must
                      correspond 1-for-1 (in order) with fastq1
            * fastqs_dir: Directory of FASTQ files of single-end reads
            * fastqi_dir: Directory of FASTQ files of interleaved
                          paired-end reads
            * fastq12_dir: Directory of FASTQ files of separate mate 1
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


class ReadsFileBase(ABC):
    _module = path.Module.ALIGN
    _step = ""
    _ext = ""

    def __init__(self, *, top_dir: str, num_cpus: int,
                 save_temp: bool, resume: bool):
        self._top_dir = path.TopDirPath.parse(top_dir)
        self._num_cpus = num_cpus
        self._save_temp = save_temp
        self._resume = resume
        self._output_files: list[path.TopDirPath] = list()
        self._index_files: list[path.TopDirPath] = list()

    @classmethod
    def write_to_temp(cls):
        return bool(cls._step)

    @property
    @abstractmethod
    def _sample(self):
        raise NotImplementedError

    @property
    def _demult(self):
        raise NotImplementedError

    @property
    def fields(self):
        fields = {"top": self._top_dir.top,
                  "module": self._module,
                  "sample": self._sample}
        if self.write_to_temp():
            fields["step"] = self._step
        return fields

    @cached_property
    def output_dir(self):
        """ Return the directory in which files (temporary or final
        output) will be written. """
        # Use the convention that the class attribute "step" is the
        # empty string if the class outputs to the final output
        # directory, and otherwise is the name of the step.
        if "step" in self.fields:
            return path.SampleStepDirPath(**self.fields)
        else:
            return path.SampleOutDirPath(**self.fields)

    @cached_property
    @abstractmethod
    def output(self):
        raise NotImplementedError

    def _setup(self):
        self.output_dir.path.mkdir(parents=True, exist_ok=True)

    @abstractmethod
    def _run(self, **kwargs):
        """ Run this step. Each class implements a unique function. """
        raise NotImplementedError

    def _clean_index_files(self):
        while self._index_files:
            self._index_files.pop().path.unlink(missing_ok=True)

    def run(self, **kwargs):
        """ Wrapper for internal run method that first checks if the
        step actually needs to be run because its output files might
        already exist. """
        if not self._resume or not self.output.path.is_file():
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
        if self.write_to_temp() and not self._save_temp:
            while self._output_files:
                self._output_files.pop().path.unlink(missing_ok=True)


class FastqBase(ReadsFileBase, ABC):
    def __init__(self, *, fastq: FastqUnit, **kwargs):
        super().__init__(**kwargs)
        self._fastq = fastq

    @property
    def _sample(self):
        return self._fastq.sample

    @property
    def _demult(self):
        return self._fastq.demult

    @staticmethod
    def _qc(fastq_paths: list[path.TopDirPath], extract: bool):
        cmd = [FASTQC_CMD]
        if extract:
            cmd.append("--extract")
        cmd.extend(fastq_paths)
        run_cmd(cmd)

    def qc_input(self, extract: bool):
        self._qc(list(self._fastq.inputs.values()), extract)

    def qc_output(self, extract: bool):
        if self._save_temp or not self.write_to_temp():
            self._qc(list(self.output.inputs.values()), extract)


class FastqTrimmer(FastqBase):
    _step = path.Step.ALIGN_TRIM

    @cached_property
    def output(self):
        return self._fastq.trans(path.ReadsInToReadsStep,
                                 **self.fields,
                                 preserve_type=True)

    @property
    def _cutadapt_output_flags(self):
        return "-o", "-p"

    @property
    def _cutadapt_output_args(self):
        return tuple(itertools.chain(*zip(self._cutadapt_output_flags,
                                          self.output.paths,
                                          strict=False)))

    def _cutadapt(self,
                  cut_q1: int,
                  cut_q2: int,
                  cut_g1: tuple[str],
                  cut_a1: tuple[str],
                  cut_g2: tuple[str],
                  cut_a2: tuple[str],
                  cut_o: int,
                  cut_e: float,
                  cut_indels: bool,
                  cut_nextseq: bool,
                  cut_discard_trimmed: bool,
                  cut_discard_untrimmed: bool,
                  cut_m: int):
        cmd = [CUTADAPT_CMD]
        cmd.extend(["--cores", self._num_cpus])
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
            if adapter and (self._fastq.paired or arg.islower()):
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
        if self._fastq.interleaved:
            cmd.append("--interleaved")
        cmd.extend(["--report", "minimal"])
        cmd.extend(self._cutadapt_output_args)
        cmd.extend(self._fastq.cutadapt_input_args)
        run_cmd(cmd)
        self._output_files.extend(self.output.inputs.values())

    def _run(self, **kwargs):
        self._cutadapt(**kwargs)


class FastqAligner(FastqBase):
    _step = path.Step.ALIGN_ALIGN
    _ext = path.SAM_EXT

    def __init__(self, *,
                 fasta: path.RefsetSeqInFilePath | path.OneRefSeqStepFilePath,
                 **kwargs):
        super().__init__(**kwargs)
        if isinstance(fasta, path.RefsetSeqInFilePath):
            if self._demult:
                raise TypeError("Got a multi-FASTA but a demultiplexed FASTQ")
        elif isinstance(fasta, path.OneRefSeqStepFilePath):
            if not self._demult:
                raise TypeError("Got a single-FASTA but a multiplexed FASTQ")
        else:
            raise TypeError(fasta)
        self._fasta = fasta

    @cached_property
    def output(self):
        # First, fasta provides either the field 'refset' or 'ref'.
        # Then, the fields 'top', 'partition', 'module', 'step',
        # and 'sample' of this class provide the output directory.
        # Lastly, the 'ext' field of this class adds the file extension.
        fields = {**self._fasta.dict(),
                  **self.fields,
                  path.EXT_KEY: self._ext}
        outputs = list()
        for fq in self._fastq.inputs.values():
            temp_inst = path.ReadsInToAlignmentStep.trans_inst(fq, **fields)
            in_type = path.AlignmentInToAlignmentStep.inverse().trans_type(
                type(temp_inst))
            in_inst = in_type.parse(temp_inst.path)
            outputs.append(in_inst)
        if not outputs:
            raise ValueError("No output files")
        if any(output != outputs[0] for output in outputs[1:]):
            raise ValueError("Inconsistent output files")
        return outputs[0]

    @cached_property
    def _bowtie2_index_files(self):
        return [path.RefSeqToBowtie2Index.trans_inst(self._fasta, ext=ext)
                for ext in path.BOWTIE2_INDEX_EXTS]

    @cached_property
    def _fasta_prefix(self):
        prefix = str(self._fasta.path.with_suffix(""))
        if indexes := [index for index in map(str, self._bowtie2_index_files)
                       if not index.startswith(prefix)]:
            logging.critical(f"Bowtie2 index files {indexes} do not start "
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

    def _bowtie2(self,
                 bt2_local: bool,
                 bt2_discordant: bool,
                 bt2_mixed: bool,
                 bt2_dovetail: bool,
                 bt2_contain: bool,
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
        """
        Run alignment with Bowtie 2 (command ```bowtie2```).

        Parameters
        ----------
        bt2_local: bool
            Whether to align reads locally (True) or end-to-end (False)
        bt2_gbar: int
            Bar gaps within this many positions from the end of a read
        bt2_dpad: int
            Width of padding on each side of the matrix, to allow gaps
        bt2_l: int
            Length of each alignment seed during multiseed alignment
        bt2_s: str (function)
            Space between alignment seeds during multiseed alignment
        bt2_d: int
            Maximum number of failed seed extensions before moving on
        bt2_r: int
            Maximum number of attempts to re-seed repetitive seeds
        bt2_score_min: str (function)
            Minimum score to output an alignment
        bt2_i: int (paired-end reads only)
            Minimum length of the fragment containing both mates
        bt2_x: int (paired-end reads only)
            Maximum length of the fragment containing both mates
        bt2_orient: str (paired-end reads only)
            Upstream/downstream mate orientations for a valid alignment
        bt2_discordant: bool (paired-end reads only)
            Whether to output reads that align discordantly
        bt2_contain: bool (paired-end reads only)
            Whether to call a mate that contains the other concordant
        bt2_dovetail: bool (paired-end reads only)
            Whether to call dovetailed alignments concordant
        bt2_mixed: bool (paired-end reads only)
            Whether to align mates individually if a pair fails to align

        Returns
        -------
        RefsetAlignmentFilePath | OneRefAlignmentFilePath
            Path of the output SAM file
        """
        if missing := self._missing_bowtie2_index_files:
            logging.critical("Bowtie2 index files do not exist:\n\n"
                             + "\n".join(map(str, missing)))
            return
        cmd = [BOWTIE2_CMD]
        # Resources
        cmd.extend(["-p", self._num_cpus])
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
        cmd.append(self._fastq.phred_arg)
        cmd.append("--ignore-quals")
        cmd.extend(["--ma", MATCH_BONUS])
        cmd.extend(["--mp", MISMATCH_PENALTY])
        cmd.extend(["--np", N_PENALTY])
        cmd.extend(["--rfg", REF_GAP_PENALTY])
        cmd.extend(["--rdg", READ_GAP_PENALTY])
        # Filtering
        cmd.extend(["--score-min", bt2_score_min])
        cmd.extend(["-I", bt2_i])
        cmd.extend(["-X", bt2_x])
        cmd.append("--no-unal")
        # Mate pair orientation
        orientations = tuple(op.value for op in MateOrientationOption)
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
        cmd.extend(self._fastq.bowtie2_inputs)
        # Run alignment.
        run_cmd(cmd)
        self._output_files.append(self.output)

    def _run(self, **kwargs):
        if self._missing_bowtie2_index_files:
            self._bowtie2_build()
        self._bowtie2(**kwargs)


class XamBase(ReadsFileBase, ABC):
    def __init__(self, *,
                 xam: (path.RefsetAlignmentInFilePath |
                       path.OneRefAlignmentInFilePath),
                 **kwargs):
        super().__init__(**kwargs)
        self.input = xam

    @property
    def _sample(self):
        return self.input.sample

    @property
    def _demult(self):
        if isinstance(self.input, path.RefsetAlignmentInFilePath):
            return False
        if isinstance(self.input, path.OneRefAlignmentInFilePath):
            return True
        raise TypeError(self.input)

    @property
    def _output_trans(self):
        return (path.AlignmentInToAlignmentStep if self.write_to_temp()
                else path.AlignmentInToAlignmentOut)

    @cached_property
    def output(self):
        return self._output_trans.trans_inst(self.input,
                                             **self.fields,
                                             ext=self._ext,
                                             preserve_type=True)

    def view(self, output: (path.OneRefAlignmentInFilePath |
                            path.OneRefAlignmentStepFilePath |
                            path.OneRefAlignmentOutFilePath |
                            path.TopDirPath)):
        if self._demult and self.input.ext == self._ext:
            self.input.path.rename(output.path)
        else:
            cmd = [SAMTOOLS_CMD, "view"]
            if self._ext == path.BAM_EXT:
                cmd.append("-b")
            cmd.extend(("-o", output, self.input))
            if not self._demult:
                cmd.append(output.ref)
            run_cmd(cmd)
        if output.path.is_file():
            self._output_files.append(output)
        else:
            logging.critical(f"Failed to view file: {self.input} -> {output}")

    @staticmethod
    def _get_index_path(xam: path.TopDirPath):
        return path.AlignmentToAlignmentIndex.trans_inst(xam, ext=path.BAI_EXT)

    def _build_index(self, xam: path.TopDirPath):
        cmd = [SAMTOOLS_CMD, "index", xam]
        run_cmd(cmd)
        xam_index = self._get_index_path(xam)
        if not xam_index.path.is_file():
            logging.critical(f"Failed to index file: {xam}")
        self._index_files.append(xam_index)


class BamIndexer(XamBase):
    _ext = path.BAM_EXT

    def __init__(self,
                 xam: (path.RefsetAlignmentInFilePath |
                       path.OneRefAlignmentInFilePath),
                 num_cpus: int,
                 resume: bool):
        super().__init__(xam=xam,
                         top_dir=xam.top,
                         num_cpus=num_cpus,
                         resume=resume,
                         save_temp=True)

    def view(self, _):
        logging.error("BamIndexer does not implement view")

    @cached_property
    def output(self):
        return self._get_index_path(self.input)

    def _run(self, **kwargs):
        self._build_index(self.input)


class SamRemoveEqualMappers(XamBase):
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
            logging.error(f"SAM file {self.input} was paired but had an odd "
                          f"number of lines ({total_lines})")

    def _iter_single(self, sam: BinaryIO, line: bytes):
        """ For each read, yield the best-scoring alignment, excluding
        reads that aligned equally well to multiple locations. """
        while line:
            if self._is_best_alignment(line):
                yield line
            line = sam.readline()

    def _remove_equal_mappers(self, buffer_length=BUFFER_LENGTH):
        with (open(self.input.path, "rb") as sami,
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

    def _run(self):
        logging.info("\nRemoving Reads Mapping Equally to Multiple Locations"
                     f" in {self.input.path}\n")
        self._remove_equal_mappers()
        self._output_files.append(self.output)


class XamSorter(XamBase):
    def _sort(self, name: bool):
        cmd = [SAMTOOLS_CMD, "sort"]
        if name:
            cmd.append("-n")
        cmd.extend(["-o", self.output.path, self.input.path])
        run_cmd(cmd)
        self._output_files.append(self.output)

    def _run(self, name: bool = False):
        logging.info(f"\nSorting {self.input.path} by Reference and Coordinate\n")
        self._sort(name)


class BamAlignSorter(XamSorter):
    _step = path.Step.ALIGN_SORT
    _ext = path.BAM_EXT


class SamVectorSorter(XamSorter):
    _module = path.Module.VECTOR
    _step = path.Step.VECTOR_SORT
    _ext = path.SAM_EXT


class BamSplitter(XamBase):
    _step = path.Step.ALIGN_SPLIT
    _ext = path.BAM_EXT

    def __init__(self, *,
                 fasta: (path.RefsetSeqInFilePath |
                         path.OneRefSeqStepFilePath),
                 **kwargs):
        super().__init__(**kwargs)
        if isinstance(fasta, path.RefsetSeqInFilePath):
            if self._demult:
                raise TypeError("Got multi-FASTA but demultiplexed BAM")
        elif isinstance(fasta, path.OneRefSeqStepFilePath):
            if not self._demult:
                raise TypeError("Got single-FASTA but multiplexed BAM")
        else:
            raise TypeError(fasta)
        self.fasta = fasta

    @cached_property
    def refs(self):
        return tuple(ref for ref, _ in FastaParser(self.fasta.path).parse())

    @cached_property
    def output(self):
        assembly_type = path.AlignmentToOneRefAlignment.trans_type(
            self._output_trans.trans_type(type(self.input)))
        output = list()
        for ref in self.refs:
            assembly_path = assembly_type(**self.fields, ref=ref, ext=self._ext)
            output.append(path.OneRefAlignmentInFilePath.parse(assembly_path))
        if self._demult and len(output) != 1:
            logging.error(f"Got {len(output)} output files for a demultiplexed "
                          "input file (expected 1)")
        return output

    def _split(self):
        for output in self.output:
            self.view(output)

    def _run(self):
        logging.info(f"\nSplitting {self.input} by reference\n")
        self._build_index(self.input)
        self._split()


class BamOutputter(XamBase):
    _ext = path.BAM_EXT

    def _run(self):
        logging.info(f"\nOutputting {self.input} to {self.output}\n")
        self.view(self.output)


class BamVectorSelector(XamBase):
    _module = path.Module.VECTOR
    _step = path.Step.VECTOR_SELECT
    _ext = path.BAM_EXT

    def __init__(self, *, ref: str, end5: int, end3: int, **kwargs):
        super().__init__(**kwargs)
        self.ref = ref
        self.end5 = end5
        self.end3 = end3

    @staticmethod
    def ref_coords(ref: str, end5: int, end3: int):
        return f"{ref}:{end5}-{end3}"

    @property
    def fields(self):
        return {**super().fields, "ref": self.ref}

    @cached_property
    def output_dir(self):
        if "step" in self.fields:
            return path.RefStepDirPath(**self.fields)
        else:
            return path.RefOutDirPath(**self.fields)

    @cached_property
    def output(self):
        if self.write_to_temp():
            trans = path.AlignmentInToRegionAlignmentStep
        else:
            trans = path.AlignmentInToRegionAlignmentOut
        return trans.trans_inst(self.input,
                                **self.fields,
                                end5=self.end5,
                                end3=self.end3,
                                ext=self._ext,
                                preserve_type=True)

    def _run(self):
        cmd = [SAMTOOLS_CMD, "view", "-h", "-o", self.output,
               self.input, self.ref_coords(self.ref, self.end5, self.end3)]
        run_cmd(cmd)


'''
primer1 = "CAGCACTCAGAGCTAATACGACTCACTATA"
primer1rc = "TATAGTGAGTCGTATTAGCTCTGAGTGCTG"
primer2 = "TGAAGAGCTGGAACGCTTCACTGA"
primer2rc = "TCAGTGAAGCGTTCCAGCTCTTCA"
adapters5 = (primer1, primer2rc)
adapters3 = (primer2, primer1rc)
'''
