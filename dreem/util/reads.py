import pathlib
from abc import ABC, abstractmethod
from enum import Enum, StrEnum
import itertools
import logging
import os
import re
from functools import cached_property
from typing import BinaryIO

from dreem.util import path
from dreem.util.cli import DEFAULT_LOCAL, DEFAULT_UNALIGNED, DEFAULT_DISCORDANT, DEFAULT_MIXED, DEFAULT_DOVETAIL, \
    DEFAULT_CONTAIN, DEFAULT_FRAG_LEN_MIN, DEFAULT_FRAG_LEN_MAX, DEFAULT_SEED_INTERVAL, \
    DEFAULT_GAP_BAR, DEFAULT_SEED_SIZE, DEFAULT_EXTENSIONS, DEFAULT_RESEED, DEFAULT_PADDING, \
    MATCH_BONUS, MISMATCH_PENALTY, N_PENALTY, REF_GAP_PENALTY, READ_GAP_PENALTY, MateOrientation
from dreem.util.cli import DEFAULT_MIN_BASE_QUALITY, DEFAULT_ILLUMINA_ADAPTER, DEFAULT_MIN_OVERLAP, DEFAULT_MAX_ERROR, \
    DEFAULT_INDELS, DEFAULT_NEXTSEQ_TRIM, DEFAULT_DISCARD_TRIMMED, DEFAULT_DISCARD_UNTRIMMED, DEFAULT_MIN_LENGTH, \
    DEFAULT_SCORE_MIN
from dreem.util.dflt import BUFFER_LENGTH
from dreem.util.excmd import FASTQC_CMD, CUTADAPT_CMD, BOWTIE2_CMD, \
    BOWTIE2_BUILD_CMD, SAMTOOLS_CMD, run_cmd
from dreem.util.seq import FastaParser

# General parameters
DEFAULT_INTERLEAVED = False
SAM_HEADER = b"@"
SAM_ALIGN_SCORE = b"AS:i:"
SAM_EXTRA_SCORE = b"XS:i:"
FASTQ_REC_LENGTH = 4
DEFAULT_MIN_MAPQ = 30

# FastQC parameters
DEFAULT_EXTRACT = False


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

    class KeyName(StrEnum):
        SINGLE = "fastqs"  # single-end reads
        INTER = "fastqi"  # interleaved paired-end reads
        MATE1 = "fastq1"  # mate 1 paired-end reads
        MATE2 = "fastq2"  # mate 2 paired-end reads

    class Bowtie2Flag(StrEnum):
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

    class KeySampleDirName(StrEnum):
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
        return tuple(map(str, cls.KeyName))

    @classmethod
    def _test_single(cls, keys: tuple[str, ...]):
        """ Return whether the arguments match a single-end FASTQ.
        No validation is performed. """
        return keys == (cls.KeyName.SINGLE,)

    @classmethod
    def _test_interleaved(cls, keys: tuple[str, ...]):
        """ Return whether the arguments match an interleaved FASTQ.
        No validation is performed. """
        return keys == (cls.KeyName.INTER,)

    @classmethod
    def _test_separate_mates(cls, keys: tuple[str, ...]):
        """ Return whether the arguments match two separate paired-end FASTQs.
        No validation is performed. """
        return keys == (cls.KeyName.MATE1, cls.KeyName.MATE2)

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
        inputs: dict[str, path.SampleReadsInFilePath |
                          path.SampleReads1InFilePath |
                          path.SampleReads2InFilePath |
                          path.DemultReadsInFilePath |
                          path.DemultReads1InFilePath |
                          path.DemultReads2InFilePath] = dict()
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
        return tuple(str(self.Bowtie2Flag[self.KeyName(key).name])
                     for key in self._input_keys)

    @property
    def bowtie2_inputs(self):
        return tuple(itertools.chain(*map(list, zip(self._bowtie2_flags,
                                                    self.paths,
                                                    strict=True))))

    def trans(self, trans: type[path.PathTypeTranslator], **kwargs):
        """
        Return a new FastqUnit by translating all the FASTQ paths using the
        given PathTypeTranslator and optionally extra arguments.

        Parameters
        ----------
        trans: PathTypeTranslator
            Translation from one FASTQ type to another
        **kwargs: Any
            Keyword arguments passed to the ```trans_inst``` method of the
            path type translator

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
        Yield a FastqUnit for each demultiplexed FASTQ file in a directory.

        Parameters
        ----------
        phred_enc: int
            Phred score encoding
        dir_path: SampleInDirPath
            Directory containing the FASTQ files
        key: str
            Keyword argument for every FASTQ file ('fastqs' or 'fastqi')

        Return
        ------
        Iterable[FastqUnit]
            One for each FASTQ file in the directory
        """
        dp = dir_path.path
        parse_type = cls.KeyDemultType[cls.KeyName(key).name].value
        for fname in dp.iterdir():
            yield cls(**{key: parse_type.parse_path(dp.joinpath(fname))},
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
                fq = path.DemultReads1InFilePath.parse_path(fpath)
                if fq.ref in mates1:
                    raise ValueError(f"Got >1 FASTQ file for ref '{fq.ref}'")
                # Add the path to the dict of mate 1 files, keyed by reference.
                mates1[fq.ref] = fq
            elif is2:
                # The file name matched an extension for a mate 2 FASTQ file.
                fq = path.DemultReads2InFilePath.parse_path(fpath)
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
            yield FastqUnit(**{key: file_type.parse_path(path_str)},
                            phred_enc=phred_enc)

    @classmethod
    def _get_sample_pairs(cls, phred_enc: int,
                          path1_strs: tuple[str],
                          path2_strs: tuple[str]):
        keys = (cls.KeyName.MATE1, cls.KeyName.MATE2)
        for path_strs in zip(path1_strs, path2_strs, strict=True):
            paths = {key.value:
                         cls.KeySampleType[key.name].value.parse_path(path_str)
                     for key, path_str in zip(keys, path_strs, strict=True)}
            yield FastqUnit(**paths, phred_enc=phred_enc)

    @classmethod
    def from_strs(cls, *, phred_enc: int, **fastq_args: tuple[str]):
        path1_strs = ()
        for fq_key, path_strs in fastq_args.items():
            try:
                # Check if fq_key is a key for a directory or file.
                dir_type = cls.KeySampleDirType[
                    cls.KeySampleDirName(fq_key).name].value
            except ValueError:
                # fq_key is not a key for a directory: assume it is a file key.
                if fq_key == cls.KeyName.MATE1 or fq_key == cls.KeyName.MATE2:
                    if not path1_strs:
                        try:
                            path1_strs = fastq_args[cls.KeyName.MATE1]
                            path2_strs = fastq_args[cls.KeyName.MATE2]
                        except KeyError:
                            raise KeyError(
                                f"Must give both {cls.KeyName.MATE1} and "
                                f"{cls.KeyName.MATE2} or neither, not only one")
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
                for dir_path in map(dir_type.parse_path, path_strs):
                    # Yield each FastqUnit depending on the FASTQ key.
                    if fq_key == cls.KeySampleDirName.MATE12:
                        # The key indicates that the directory contains paired
                        # FASTQ files split into mate 1 and mate 2 reads.
                        yield from cls._get_demult_pairs(phred_enc=phred_enc,
                                                         dir_path=dir_path)
                    else:
                        # The key indicates that the directory contains FASTQ
                        # files that may be processed individually.
                        yield from cls._get_demult_files(phred_enc=phred_enc,
                                                         dir_path=dir_path,
                                                         key=fq_key)


class ReadsFileBase(ABC):
    partition = path.Partition.TEMP
    module = path.Module.ALIGN
    step = ""
    ext = ""

    def __init__(self, top_dir: path.TopDirPath, max_cpus: int):
        self.top_dir = top_dir
        self.max_cpus = max_cpus

    @property
    @abstractmethod
    def sample(self):
        raise NotImplementedError

    @property
    def demult(self):
        raise NotImplementedError

    @property
    def fields(self):
        fields = {"top": self.top_dir.top,
                  "partition": self.partition,
                  "module": self.module,
                  "sample": self.sample}
        if self.partition == path.Partition.TEMP:
            fields["step"] = self.step
        return fields

    @cached_property
    def output_dir(self):
        if self.partition == path.Partition.TEMP:
            return path.SampleTempDirPath(**self.fields)
        if self.partition == path.Partition.OUTPUT:
            return path.SampleOutDirPath(**self.fields)
        raise ValueError(self.partition)

    @cached_property
    @abstractmethod
    def output(self):
        raise NotImplementedError

    def setup(self):
        self.output_dir.path.mkdir(parents=True, exist_ok=True)

    @abstractmethod
    def run(self):
        raise NotImplementedError

    @abstractmethod
    def clean(self):
        raise NotImplementedError


class FastqBase(ReadsFileBase):
    def __init__(self,
                 top_dir: path.TopDirPath,
                 max_cpus: int,
                 fastq: FastqUnit):
        super().__init__(top_dir, max_cpus)
        self.fastq = fastq

    @property
    def sample(self):
        return self.fastq.sample

    @property
    def demult(self):
        return self.fastq.demult

    def qc(self, extract: bool = DEFAULT_EXTRACT):
        cmd = [FASTQC_CMD]
        if extract:
            cmd.append("--extract")
        cmd.extend(self.fastq.paths)
        run_cmd(cmd)


class FastqTrimmer(FastqBase):
    step = path.TempStep.ALIGN_TRIM

    @cached_property
    def output(self):
        return self.fastq.trans(path.ReadsInToReadsTemp,
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
                  qual1: int = DEFAULT_MIN_BASE_QUALITY,
                  qual2: int = 0,
                  adapters15: tuple[str] = (),
                  adapters13: tuple[str] = (DEFAULT_ILLUMINA_ADAPTER,),
                  adapters25: tuple[str] = (),
                  adapters23: tuple[str] = (DEFAULT_ILLUMINA_ADAPTER,),
                  min_overlap: int = DEFAULT_MIN_OVERLAP,
                  max_error: float = DEFAULT_MAX_ERROR,
                  indels: bool = DEFAULT_INDELS,
                  nextseq_trim: bool = DEFAULT_NEXTSEQ_TRIM,
                  discard_trimmed: bool = DEFAULT_DISCARD_TRIMMED,
                  discard_untrimmed: bool = DEFAULT_DISCARD_UNTRIMMED,
                  min_length: bool = DEFAULT_MIN_LENGTH):
        cmd = [CUTADAPT_CMD]
        cmd.extend(["--cores", self.max_cpus])
        if nextseq_trim:
            if qual1 > 0:
                cmd.extend(["--nextseq-trim", qual1])
        else:
            if qual1 > 0:
                cmd.extend(["-q", qual1])
            if qual2 > 0:
                cmd.extend(["-Q", qual2])
        adapters = {"g": adapters15, "a": adapters13,
                    "G": adapters25, "A": adapters23}
        for arg, adapter in adapters.items():
            if adapter and (self.fastq.paired or arg.islower()):
                for adapt in adapter:
                    cmd.extend([f"-{arg}", adapt])
        if min_overlap >= 0:
            cmd.extend(["-O", min_overlap])
        if max_error >= 0:
            cmd.extend(["-e", max_error])
        if not indels:
            cmd.append("--no-indels")
        if discard_trimmed:
            cmd.append("--discard-trimmed")
        if discard_untrimmed:
            cmd.append("--discard-untrimmed")
        if min_length:
            cmd.extend(["-m", min_length])
        cmd.extend(["--report", "minimal"])
        if self.fastq.interleaved:
            cmd.append("--interleaved")
        cmd.extend(self._cutadapt_output_args)
        cmd.extend(self.fastq.cutadapt_input_args)
        run_cmd(cmd)
        return self.output

    def run(self, **kwargs):
        self.setup()
        return self._cutadapt(**kwargs)

    def clean(self):
        for out_path in self.output.paths:
            out_path.unlink()


class FastqAligner(FastqBase):
    step = path.TempStep.ALIGN_ALIGN
    ext = path.SAM_EXT

    def __init__(self,
                 top_dir: path.TopDirPath,
                 max_cpus: int,
                 fastq: FastqUnit,
                 fasta: path.RefsetSeqInFilePath | path.OneRefSeqTempFilePath):
        super().__init__(top_dir, max_cpus, fastq)
        if isinstance(fasta, path.RefsetSeqInFilePath):
            if self.demult:
                raise TypeError("Got a multi-FASTA but a demultiplexed FASTQ")
        elif isinstance(fasta, path.OneRefSeqTempFilePath):
            if not self.demult:
                raise TypeError("Got a single-FASTA but a multiplexed FASTQ")
        else:
            raise TypeError(fasta)
        self.fasta = fasta

    @property
    def fasta_prefix(self):
        return self.fasta.path.with_suffix("")

    @cached_property
    def output(self):
        # fasta provides either 'refset' or 'ref'
        # output_dir provides 'top', 'partition', 'module', 'step', and 'sample'
        # ext provides the file extension
        fields = {**self.fasta.dict(),
                  **self.fields,
                  path.EXT_KEY: self.ext}
        outputs = list()
        for fq in self.fastq.inputs.values():
            temp_inst = path.ReadsInToAlignmentTemp.trans_inst(fq, **fields)
            in_type = path.AlignmentInToAlignmentTemp.inverse().trans_type(
                type(temp_inst))
            in_inst = in_type.parse_path(temp_inst.path)
            outputs.append(in_inst)
        if not outputs:
            raise ValueError("No output files")
        if any(output != outputs[0] for output in outputs[1:]):
            raise ValueError("Inconsistent output files")
        return outputs[0]

    def _bowtie2_build(self):
        """ Build an index of a reference genome using Bowtie 2. """
        cmd = [BOWTIE2_BUILD_CMD, "-q", self.fasta.path, self.fasta_prefix]
        run_cmd(cmd)

    def _bowtie2(self,
                 local=DEFAULT_LOCAL,
                 unaligned=DEFAULT_UNALIGNED,
                 discordant=DEFAULT_DISCORDANT,
                 mixed=DEFAULT_MIXED,
                 dovetail=DEFAULT_DOVETAIL,
                 contain=DEFAULT_CONTAIN,
                 score_min=DEFAULT_SCORE_MIN,
                 frag_len_min=DEFAULT_FRAG_LEN_MIN,
                 frag_len_max=DEFAULT_FRAG_LEN_MAX,
                 gap_bar=DEFAULT_GAP_BAR,
                 seed_size=DEFAULT_SEED_SIZE,
                 seed_interval=DEFAULT_SEED_INTERVAL,
                 extensions=DEFAULT_EXTENSIONS,
                 reseed=DEFAULT_RESEED,
                 padding=DEFAULT_PADDING,
                 orientation=MateOrientation.FR):
        """
        Run alignment with Bowtie 2 (command ```bowtie2```).

        Parameters
        ----------
        local: bool
            Whether to find local (True) or end-to-end (False) read alignments
        gap_bar: int
            Forbid gaps within this many positions from the end of a read
        padding: int
            Number of positions to pad on each side of the matrix to allow gaps
        seed_size: int
            Length of the alignment seeds during multiseed alignment
        seed_interval: str (function)
            Interval between alignment seeds during multiseed alignment

        extensions: int
            Maximum number of failed seed extensions before moving on
        reseed: int
            Maximum number of attempts to re-seed repetitive seeds

        score_min: str (function)
            Minimum score to consider an alignment "good enough" to output
        frag_len_min: int (paired-end reads only)
            Minimum length of the contiguous fragment containing both mates
        frag_len_max: int (paired-end reads only)
            Maximum length of the contiguous fragment containing both mates
        unaligned: bool
            Whether to output reads that do not align

        orientation: str (paired-end reads only)
            Upstream/downstream mate orientations for a valid alignment
        discordant: bool (paired-end reads only)
            Whether to output reads that align discordantly
        contain: bool (paired-end reads only)
            Whether to consider a mate that contains the other to be concordant
        dovetail: bool (paired-end reads only)
            Whether to consider dovetailed alignments to be concordant
        mixed: bool (paired-end reads only)
            Whether to output a mate that aligns if the other mate does not

        Returns
        -------
        RefsetAlignmentFilePath | OneRefAlignmentFilePath
            Path of the output SAM file
        """
        cmd = [BOWTIE2_CMD]
        # Resources
        cmd.extend(["-p", self.max_cpus])
        # Alignment
        if local:
            cmd.append("--local")
        if gap_bar:
            cmd.extend(["--gbar", gap_bar])
        if padding:
            cmd.extend(["--dpad", padding])
        if seed_size:
            cmd.extend(["-L", seed_size])
        if seed_interval:
            cmd.extend(["-i", seed_interval])
        if extensions:
            cmd.extend(["-D", extensions])
        if reseed:
            cmd.extend(["-R", reseed])
        # Scoring
        cmd.append(self.fastq.phred_arg)
        cmd.append("--ignore-quals")
        cmd.extend(["--ma", MATCH_BONUS])
        cmd.extend(["--mp", MISMATCH_PENALTY])
        cmd.extend(["--np", N_PENALTY])
        cmd.extend(["--rfg", REF_GAP_PENALTY])
        cmd.extend(["--rdg", READ_GAP_PENALTY])
        # Filtering
        if score_min:
            cmd.extend(["--score-min", score_min])
        if frag_len_min:
            cmd.extend(["-I", frag_len_min])
        if frag_len_max:
            cmd.extend(["-X", frag_len_max])
        if not unaligned:
            cmd.append("--no-unal")
        # Mate pair orientation
        if orientation in set(MateOrientation):
            cmd.append(f"--{orientation}")
        if not discordant:
            cmd.append("--no-discordant")
        if not contain:
            cmd.append("--no-contain")
        if dovetail:
            cmd.append("--dovetail")
        if not mixed:
            cmd.append("--no-mixed")
        # Formatting
        cmd.append("--xeq")
        # Inputs and outputs
        cmd.extend(["-S", self.output.path])
        cmd.extend(["-x", self.fasta_prefix])
        cmd.extend(self.fastq.bowtie2_inputs)
        # Run alignment
        run_cmd(cmd)
        return self.output

    def run(self, **kwargs):
        self.setup()
        self._bowtie2_build()
        return self._bowtie2(**kwargs)

    def clean(self):
        self.output.path.unlink()


class XamBase(ReadsFileBase):
    def __init__(self,
                 top_dir: path.TopDirPath,
                 max_cpus: int,
                 xam: path.RefsetAlignmentInFilePath |
                      path.OneRefAlignmentInFilePath):
        super().__init__(top_dir, max_cpus)
        self.xam = xam

    @property
    def sample(self):
        return self.xam.sample

    @property
    def demult(self):
        if isinstance(self.xam, path.RefsetAlignmentInFilePath):
            return False
        if isinstance(self.xam, path.OneRefAlignmentInFilePath):
            return True
        raise TypeError(self.xam)

    @cached_property
    def output(self):
        if self.partition == path.Partition.TEMP:
            trans = path.AlignmentInToAlignmentTemp
        elif self.partition == path.Partition.OUTPUT:
            trans = path.AlignmentInToAlignmentOut
        else:
            raise ValueError(self.partition)
        return trans.trans_inst(self.xam,
                                **self.fields,
                                ext=self.ext,
                                preserve_type=True)

    @staticmethod
    def _build_index(xam: path.TopDirPath):
        cmd = [SAMTOOLS_CMD, "index", xam.path]
        run_cmd(cmd)
        xam_index = xam.replace(ext=path.BAI_EXT)
        if not xam_index.path.is_file():
            raise FileNotFoundError(xam_index.path)
        return xam_index

    def index_input(self):
        return self._build_index(self.xam)

    def index_output(self):
        return self._build_index(self.output)

    def clean(self):
        self.output.path.unlink(missing_ok=True)


class SamRemoveEqualMappers(XamBase):
    step = path.TempStep.ALIGN_REMEQ
    ext = path.SAM_EXT

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
        for line2 in sam:
            if self._is_best_alignment(line) or self._is_best_alignment(line2):
                yield b"".join((line, line2))
            line = sam.readline()

    def _iter_single(self, sam: BinaryIO, line: bytes):
        while line:
            if self._is_best_alignment(line):
                yield line
            line = sam.readline()

    def _remove_equal_mappers(self, buffer_length=BUFFER_LENGTH):
        with (open(self.xam.path, "rb") as sami,
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
        return self.output

    def run(self):
        logging.info("\nRemoving Reads Mapping Equally to Multiple Locations"
                     f" in {self.xam.path}\n")
        self.setup()
        return self._remove_equal_mappers()


class XamSorter(XamBase):
    def _sort(self, name: bool):
        cmd = [SAMTOOLS_CMD, "sort"]
        if name:
            cmd.append("-n")
        cmd.extend(["-o", self.output.path, self.xam.path])
        run_cmd(cmd)
        return self.output

    def run(self, name: bool = False):
        logging.info(f"\nSorting {self.xam.path} by Reference and Coordinate\n")
        self.setup()
        return self._sort(name)


class BamAlignSorter(XamSorter):
    step = path.TempStep.ALIGN_SORT
    ext = path.BAM_EXT


class SamVectorSorter(XamSorter):
    module = path.Module.VECTOR
    step = path.TempStep.VECTOR_SORT
    ext = path.SAM_EXT


class BamSplitter(XamBase):
    partition = path.Partition.OUTPUT
    step = path.TempStep.ALIGN_SPLIT
    ext = path.BAM_EXT

    def __init__(self,
                 top_dir: path.TopDirPath,
                 max_cpus: int,
                 xam: path.RefsetAlignmentInFilePath |
                      path.OneRefAlignmentInFilePath,
                 fasta: path.RefsetSeqInFilePath |
                        path.OneRefSeqTempFilePath):
        super().__init__(top_dir, max_cpus, xam)
        if isinstance(fasta, path.RefsetSeqInFilePath):
            if self.demult:
                raise TypeError("Got multi-FASTA but demultiplexed BAM")
        elif isinstance(fasta, path.OneRefSeqTempFilePath):
            if not self.demult:
                raise TypeError("Got single-FASTA but multiplexed BAM")
        else:
            raise TypeError(fasta)
        self.fasta = fasta

    @cached_property
    def refs(self):
        return tuple(ref for ref, _ in FastaParser(self.fasta.path).parse())

    def _get_cmd(self, output: pathlib.Path, ref: str = ""):
        cmd = [SAMTOOLS_CMD, "view"]
        if self.ext == path.BAM_EXT:
            cmd.append("-b")
        cmd.extend(("-o", output, self.xam.path))
        if ref:
            cmd.append(ref)
        return cmd

    def _extract_one_ref(self, ref: str):
        output = path.OneRefAlignmentOutFilePath(**self.fields, ref=ref,
                                                 ext=self.ext)
        cmd = self._get_cmd(output.path, ref=ref)
        run_cmd(cmd)
        return output

    def _split_xam(self):
        self.index_input()
        return tuple(map(self._extract_one_ref, self.refs))

    def _move_xam(self):
        if self.ext == self.xam.ext:
            os.rename(self.xam.path, self.output.path)
        else:
            cmd = self._get_cmd(self.output.path)
            run_cmd(cmd)
        return self.output,

    def _build_outputs(self):
        if self.demult:
            return self._move_xam()
        else:
            return self._split_xam()

    def _build_indexes(self, bams: tuple[path.TopDirPath]):
        return tuple(map(self._build_index, bams))

    def run(self):
        logging.info(f"\nSplitting {self.xam.path} by reference\n")
        self.setup()
        bams = self._build_outputs()
        self._build_indexes(bams)
        return bams

    def clean(self):
        return


class BamVectorSelector(XamBase):
    module = path.Module.VECTOR
    step = path.TempStep.VECTOR_SELECT
    ext = path.BAM_EXT

    def __init__(self,
                 top_dir: path.TopDirPath,
                 max_cpus: int,
                 xam: path.OneRefAlignmentInFilePath,
                 ref: str,
                 end5: int,
                 end3: int):
        super().__init__(top_dir, max_cpus, xam)
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
        if self.partition == path.Partition.TEMP:
            return path.RefTempDirPath(**self.fields)
        if self.partition == path.Partition.OUTPUT:
            return path.RefOutDirPath(**self.fields)
        raise ValueError(self.partition)

    @cached_property
    def output(self):
        if self.partition == path.Partition.TEMP:
            trans = path.AlignmentInToRegionAlignmentTemp
        elif self.partition == path.Partition.OUTPUT:
            trans = path.AlignmentInToRegionAlignmentOut
        else:
            raise ValueError(self.partition)
        return trans.trans_inst(self.xam,
                                **self.fields,
                                end5=self.end5,
                                end3=self.end3,
                                ext=self.ext,
                                preserve_type=True)

    def _select(self):
        cmd = [SAMTOOLS_CMD, "view", "-h", "-o", self.output.path,
               self.xam.path, self.ref_coords(self.ref, self.end5, self.end3)]
        run_cmd(cmd)
        return self.output

    def run(self):
        self.setup()
        return self._select()


'''
primer1 = "CAGCACTCAGAGCTAATACGACTCACTATA"
primer1rc = "TATAGTGAGTCGTATTAGCTCTGAGTGCTG"
primer2 = "TGAAGAGCTGGAACGCTTCACTGA"
primer2rc = "TCAGTGAAGCGTTCCAGCTCTTCA"
adapters5 = (primer1, primer2rc)
adapters3 = (primer2, primer1rc)
'''
