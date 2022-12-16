from collections import defaultdict
from io import BufferedReader
import itertools
from multiprocessing import Pool
import os
import re
import shutil

from typing import List, Optional, Tuple

from dreem.util.util import FASTQC_CMD, CUTADAPT_CMD, BOWTIE2_CMD, PASTE_CMD, BASEN, \
    BOWTIE2_BUILD_CMD, TEMP_DIR, OUTPUT_DIR, DEFAULT_PROCESSES, run_cmd, \
    try_remove, try_rmdir, switch_directory, PHRED_ENCODING, SAMTOOLS_CMD, FastaParser


# General parameters
DEFAULT_INTERLEAVED = False
SAM_HEADER = b"@"
SAM_ALIGN_SCORE = b"AS:i:"
SAM_EXTRA_SCORE = b"XS:i:"
FASTQ_REC_LENGTH = 4
DEFAULT_MIN_MAPQ = 30

# FastQC parameters
DEFAULT_EXTRACT = False

# Cutadapt parameters
DEFAULT_MIN_BASE_QUALITY = 25
DEFAULT_ILLUMINA_ADAPTER = "AGATCGGAAGAGC"
DEFAULT_MIN_OVERLAP = 6
DEFAULT_MAX_ERROR = 0.1
DEFAULT_INDELS = True
DEFAULT_NEXTSEQ = True
DEFAULT_DISCARD_TRIMMED = False
DEFAULT_DISCARD_UNTRIMMED = False
DEFAULT_MIN_LENGTH = 50

# Bowtie 2 parameters
DEFAULT_LOCAL = True
DEFAULT_UNALIGNED = False
DEFAULT_DISCORDANT = False
DEFAULT_MIXED = False
DEFAULT_DOVETAIL = False
DEFAULT_CONTAIN = True
DEFAULT_FRAG_LEN_MIN = 0
DEFAULT_FRAG_LEN_MAX = 300  # maximum length of a 150 x 150 read
DEFAULT_N_CEILING = "L,0,0.05"
DEFAULT_SEED_INTERVAL = "L,1,0.1"
DEFAULT_GAP_BAR = 4
DEFAULT_SEED_SIZE = 20
DEFAULT_EXTENSIONS = 5
DEFAULT_RESEED = 1
DEFAULT_PADDING = 4
DEFAULT_ALIGN_THREADS = os.cpu_count()
MATCH_BONUS = "1"
MISMATCH_PENALTY = "1,1"
N_PENALTY = "0"
REF_GAP_PENALTY = "0,1"
READ_GAP_PENALTY = "0,1"
IGNORE_QUALS = True
SAM_EXT = ".sam"
BAM_EXT = ".bam"

# Text processing parameters
DEFAULT_BUFFER_LENGTH = 10000


def get_diffs(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences were different lengths: "
                         f"'{seq1}' and '{seq2}'")
    diffs = [i for i, (x1, x2) in enumerate(zip(seq1, seq2)) if x1 != x2]
    return diffs


def get_fastq_name(fastq: str, fastq2: str = ""):
    exts = (".fq", ".fastq")
    reads = ("mate", "r")
    base = os.path.basename(fastq)
    counts = {ext: count for ext in exts if (count := base.count(ext))}
    if not counts:
        raise ValueError(f"{fastq} had no FASTQ extension")
    if sum(counts.values()) > 1:
        raise ValueError(f"{fastq} had multiple FASTQ extensions")
    name = base[:base.index(list(counts)[0])]
    if fastq2:
        name2 = get_fastq_name(fastq2)
        diffs = get_diffs(name, name2)
        if len(diffs) != 1:
            raise ValueError("FASTQ names must differ by exactly 1 character.")
        idx = diffs[0]
        if name[idx] != "1" or name2[idx] != "2":
            raise ValueError("FASTQs must be named '1' and '2', respectively.")
        name = name[:idx]
        for read in reads:
            if name.rstrip("_").lower().endswith(read):
                name = name[:name.lower().rindex(read)]
                break
        name = name.rstrip("_")
    return name


def get_fastq_dir(fastq: str, fastq2: Optional[str] = None):
    fq_dir = os.path.dirname(fastq)
    if fastq2 and os.path.dirname(fastq2) != fq_dir:
        raise ValueError("FASTQs are not in the same directory.")
    return os.path.basename(fq_dir)


def get_fastq_pairs(fq_dir: str):
    lengths = defaultdict(set)
    for file in os.listdir(fq_dir):
        lengths[len(file)].add(file)
    pairs = dict()
    for fq_files in lengths.values():
        while fq_files:
            file1 = fq_files.pop()
            diffs = defaultdict(set)
            for file2 in fq_files:
                diffs[len(get_diffs(file1, file2))].add(file2)
            try:
                file2 = diffs[1].pop()
            except KeyError:
                raise ValueError(f"Found no mate for {file1} in {fq_dir}")
            if diffs[1]:
                raise ValueError(f"Found >1 mate for {file1} in {fq_dir}")
            fq_files.pop(file2)
            pair = (os.path.join(fq_dir, file1), os.path.join(fq_dir, file2))
            pairs[get_fastq_name(file1, file2)] = pair
    return pairs


class SeqFileBase(object):
    _operation_dir = ""

    def __init__(self,
                 root_dir: str,
                 ref_file: str,
                 sample: str) -> None:
        self._root_dir = root_dir
        self._ref_file = ref_file
        self._sample = sample
    
    @property
    def operation_dir(self):
        if not self._operation_dir:
            raise NotImplementedError
        return self._operation_dir
    
    @property
    def root_dir(self):
        return self._root_dir
    
    @property
    def ref_file(self):
        return self._ref_file
    
    @property
    def ref_prefix(self):
        return os.path.splitext(self.ref_file)[0]
    
    @property
    def ref_filename(self):
        return os.path.basename(self.ref_prefix)
    
    @property
    def sample(self):
        return self._sample
    
    def _get_dir(self, dest: str):
        return os.path.join(self.root_dir, dest,
                            self.operation_dir, self.sample)
    
    @property
    def output_dir(self):
        return self._get_dir(TEMP_DIR)
    
    def _make_output_dir(self):
        os.makedirs(self.output_dir, exist_ok=True)
    
    def run(self, *args, **kwargs):
        raise NotImplementedError
    
    def clean(self):
        shutil.rmtree(self.output_dir, ignore_errors=True)


class FastqBase(SeqFileBase):
    def __init__(self,
                 root_dir: str,
                 ref_file: str,
                 sample: str,
                 paired: bool,
                 fastq: str,
                 fastq2: str = "",
                 interleave: bool = True,
                 encoding: int = PHRED_ENCODING) -> None:
        if fastq2 and not paired:
            raise ValueError("fastq2 can only be given if reads are paired")
        super().__init__(root_dir, ref_file, sample)
        self._fq_in = fastq
        self._fq2_in = fastq2
        self._paired = paired
        self._interleave_out = interleave
        self._encoding = encoding
        self._name = get_fastq_name(fastq, fastq2)
    
    @property
    def fq_in(self):
        return self._fq_in
    
    @property
    def fq1_in(self):
        return self._fq_in
        
    @property
    def fq2_in(self):
        return self._fq2_in

    @property
    def paired(self):
        return self._paired

    @property
    def interleaved_in(self):
        return self.paired and not self.fq2_in
    
    @property
    def interleave_out(self):
        return self._interleave_out and self.paired

    @property
    def encoding(self):
        return self._encoding
    
    @property
    def encoding_arg(self):
        return f"--phred{self.encoding}"
    
    @property
    def name(self):
        return self._name
    
    @property
    def inputs(self):
        fqs = [self.fq_in]
        if self.fq2_in:
            fqs.append(self.fq2_in)
        return fqs
    
    @property
    def outputs(self):
        if self.paired and not self.interleave_out:
            fqs = ["{self.name}_R{r}" for r in ("1", "2")]
        else:
            fqs = [self.name]
        return [os.path.join(self.output_dir, f"{fq}.fq") for fq in fqs]
    
    @staticmethod
    def _qc(files: List[str], extract: bool):
        cmd = [FASTQC_CMD]
        if extract:
            cmd.append("--extract")
        cmd.extend(files)
        run_cmd(cmd)
    
    def qc_inputs(self, extract: bool = DEFAULT_EXTRACT):
        return self._qc(self.inputs, extract)

    def qc_outputs(self, extract: bool = DEFAULT_EXTRACT):
        return self._qc(self.outputs, extract)


class FastqTrimmer(FastqBase):
    _operation_dir = "alignment/1_trim"

    def _cutadapt(self,
                  qual1=DEFAULT_MIN_BASE_QUALITY,
                  qual2=None,
                  adapters15: Tuple[str] = (),
                  adapters13: Tuple[str] = (DEFAULT_ILLUMINA_ADAPTER,),
                  adapters25: Tuple[str] = (),
                  adapters23: Tuple[str] = (DEFAULT_ILLUMINA_ADAPTER,),
                  min_overlap=DEFAULT_MIN_OVERLAP,
                  max_error=DEFAULT_MAX_ERROR,
                  indels=DEFAULT_INDELS,
                  nextseq=DEFAULT_NEXTSEQ,
                  discard_trimmed=DEFAULT_DISCARD_TRIMMED,
                  discard_untrimmed=DEFAULT_DISCARD_UNTRIMMED,
                  min_length=DEFAULT_MIN_LENGTH,
                  cores=DEFAULT_PROCESSES):
        cmd = [CUTADAPT_CMD]
        if cores >= 0:
            cmd.extend(["--cores", str(cores)])
        if nextseq:
            nextseq_qual = qual1 if qual1 else DEFAULT_MIN_BASE_QUALITY
            cmd.extend(["--nextseq-trim", str(nextseq_qual)])
        else:
            if qual1 is not None:
                cmd.extend(["-q", str(qual1)])
            if qual2 is not None:
                self.fq2_in
                cmd.extend(["-Q", str(qual2)])
        adapters = {"g": adapters15, "a": adapters13,
                    "G": adapters25, "A": adapters23}
        for arg, adapter in adapters.items():
            if adapter and (arg.islower() or self.paired):
                if isinstance(adapter, str):
                    adapter = (adapter,)
                if not isinstance(adapter, tuple):
                    raise ValueError("adapters must be str or tuple")
                for adapt in adapter:
                    cmd.extend([f"-{arg}", adapt])
        if min_overlap >= 0:
            cmd.extend(["-O", str(min_overlap)])
        if max_error >= 0:
            cmd.extend(["-e", str(max_error)])
        if not indels:
            cmd.append("--no-indels")
        if discard_trimmed:
            cmd.append("--discard-trimmed")
        if discard_untrimmed:
            cmd.append("--discard-untrimmed")
        if min_length:
            cmd.extend(["-m", str(min_length)])
        cmd.extend(["--report", "minimal"])
        if self.interleaved_in or self.interleave_out:
            cmd.append("--interleaved")
        for flag, output in zip(["-o", "-p"], self.outputs):
            cmd.extend([flag, output])
        cmd.extend(self.inputs)
        run_cmd(cmd)
        return self.outputs
    
    def run(self, **kwargs):
        print(f"\nTrimming Adapters from {self.fq_in}\n")
        self._make_output_dir()
        return self._cutadapt(**kwargs)


class FastqAligner(FastqBase):
    _operation_dir = "alignment/2_align"

    @property
    def sam_out(self):
        return os.path.join(self.output_dir, f"{self.ref_filename}.sam")
    
    @property
    def outputs(self):
        return [self.sam_out]

    def _bowtie2_build(self):
        """
        Build an index of a reference genome using Bowtie 2.
        :param ref: (str) path to the reference genome FASTA file
        :return: None
        """
        cmd = [BOWTIE2_BUILD_CMD, "-q", self.ref_file, self.ref_prefix]
        run_cmd(cmd)
    
    def _bowtie2(self,
                 local=DEFAULT_LOCAL,
                 unaligned=DEFAULT_UNALIGNED,
                 discordant=DEFAULT_DISCORDANT,
                 mixed=DEFAULT_MIXED,
                 dovetail=DEFAULT_DOVETAIL,
                 contain=DEFAULT_CONTAIN,
                 frag_len_min=DEFAULT_FRAG_LEN_MIN,
                 frag_len_max=DEFAULT_FRAG_LEN_MAX,
                 n_ceil=DEFAULT_N_CEILING,
                 gap_bar=DEFAULT_GAP_BAR,
                 seed_size=DEFAULT_SEED_SIZE,
                 seed_interval=DEFAULT_SEED_INTERVAL,
                 extensions=DEFAULT_EXTENSIONS,
                 reseed=DEFAULT_RESEED,
                 padding=DEFAULT_PADDING,
                 threads=DEFAULT_ALIGN_THREADS):
        cmd = [BOWTIE2_CMD]
        if self.interleaved_in:
            cmd.extend(["--interleaved", self.fq_in])
        elif self.paired:
            cmd.extend(["-1", self.fq1_in, "-2", self.fq2_in])
        else:
            cmd.extend(["-U", self.fq_in])
        cmd.extend(["-x", self.ref_prefix])
        cmd.extend(["-S", self.sam_out])
        cmd.append(self.encoding_arg)
        cmd.append("--xeq")
        cmd.extend(["--ma", MATCH_BONUS])
        cmd.extend(["--mp", MISMATCH_PENALTY])
        cmd.extend(["--np", N_PENALTY])
        cmd.extend(["--rfg", REF_GAP_PENALTY])
        cmd.extend(["--rdg", READ_GAP_PENALTY])
        if local:
            cmd.append("--local")
        if not unaligned:
            cmd.append("--no-unal")
        if not discordant:
            cmd.append("--no-discordant")
        if not mixed:
            cmd.append("--no-mixed")
        if dovetail:
            cmd.append("--dovetail")
        if not contain:
            cmd.append("--no-contain")
        if frag_len_min:
            cmd.extend(["-I", str(frag_len_min)])
        if frag_len_max:
            cmd.extend(["-X", str(frag_len_max)])
        if n_ceil:
            cmd.extend(["--n-ceil", n_ceil])
        if gap_bar:
            cmd.extend(["--gbar", str(gap_bar)])
        if seed_size:
            cmd.extend(["-L", str(seed_size)])
        if seed_interval:
            cmd.extend(["-i", str(seed_interval)])
        if extensions:
            cmd.extend(["-D", str(extensions)])
        if reseed:
            cmd.extend(["-R", str(reseed)])
        if padding:
            cmd.extend(["--dpad", str(padding)])
        if threads:
            cmd.extend(["-p", str(threads)])
        if IGNORE_QUALS:
            cmd.append("--ignore-quals")
        run_cmd(cmd)
        return self.sam_out
    
    def run(self, **kwargs):
        print(f"\nAligning Reads {self.fq_in} to Reference {self.ref_file}\n")
        self._make_output_dir()
        self._bowtie2_build()
        return self._bowtie2(**kwargs)


class SamBase(SeqFileBase):
    def __init__(self,
                 root_dir: str,
                 ref_file: str,
                 sample: str,
                 xam_file: str) -> None:
        super().__init__(root_dir, ref_file, sample)
        self._ref_names = {ref for ref, _ in FastaParser(ref_file).parse()}
        self._xam_in = xam_file
    
    @property
    def refs(self):
        return self._ref_names
    
    @property
    def xam_in(self):
        return self._xam_in
    
    @property
    def name(self):
        return os.path.splitext(os.path.basename(self.xam_in))[0]
    
    @property
    def sam_name(self):
        return f"{self.name}.sam"

    @property
    def bam_name(self):
        return f"{self.name}{BAM_EXT}"
    
    @property
    def sam_out(self):
        return os.path.join(self.output_dir, self.sam_name)
    
    @property
    def bam_out(self):
        return os.path.join(self.output_dir, self.bam_name)
    
    def _get_bam_split(self, ref: bytes):
        return os.path.join(self.output_dir, f"{ref}{BAM_EXT}")
    
    @staticmethod
    def _index_bam(bam_file: str):
        cmd = [SAMTOOLS_CMD, "index", bam_file]
        run_cmd(cmd)

    def index_bam_in(self):
        self._index_bam(self.xam_in)
    
    def index_bam_out(self):
        self._index_bam(self.bam_out)


class SamRemoveEqualMappers(SamBase):
    _operation_dir = "alignment/3_rem"

    pattern_a = re.compile(SAM_ALIGN_SCORE + rb"(\d+)")
    pattern_x = re.compile(SAM_EXTRA_SCORE + rb"(\d+)")

    def __init__(self, root_dir: str, ref_file: str, sample: str,
                 xam_file: str, paired: bool) -> None:
        super().__init__(root_dir, ref_file, sample, xam_file)
        self._paired = paired
    
    @property
    def paired(self):
        return self._paired

    @staticmethod
    def get_score(line: bytes, ptn: re.Pattern[bytes]):
        return (float(match.groups()[0])
                if (match := ptn.search(line)) else None)
    
    @classmethod
    def is_best_alignment(cls, line: bytes):
        return ((score_x := cls.get_score(line, cls.pattern_x)) is None
                or score_x < cls.get_score(line, cls.pattern_a))
    
    def _iter_paired(self, sam: BufferedReader, line: bytes):
        for line2 in sam:
            if self.is_best_alignment(line) or self.is_best_alignment(line2):
                yield b"".join((line, line2))
            line = sam.readline()
    
    def _iter_single(self, sam: BufferedReader, line: bytes):
        while line:
            if self.is_best_alignment(line):
                yield(line)
            line = sam.readline()

    def _remove_equal_mappers(self, buffer_length=DEFAULT_BUFFER_LENGTH):
        with open(self.xam_in, "rb") as sami, open(self.sam_out, "wb") as samo:
            # Copy the header from the input to the output SAM file.
            while (line := sami.readline()).startswith(SAM_HEADER):
                samo.write(line)
            iter_sam = self._iter_paired if self.paired else self._iter_single
            lines = iter_sam(sami, line)
            while text := b"".join(itertools.islice(lines, buffer_length)):
                samo.write(text)
        return self.sam_out
    
    def run(self):
        print("\nRemoving Reads Mapping Equally to Multiple Locations in "
              f"{self.xam_in}\n")
        self._make_output_dir()
        return self._remove_equal_mappers()


class SamSorter(SamBase):
    _operation_dir = "alignment/4_sort"

    def _sort(self, name: bool = False):
        cmd = [SAMTOOLS_CMD, "sort"]
        if name:
            cmd.append("-n")
        cmd.extend(["-o", self.bam_out, self.xam_in])
        run_cmd(cmd)
        return self.bam_out
    
    def run(self, name: bool = False):
        print(f"\nSorting {self.xam_in} by Reference and Coordinate\n")
        self._make_output_dir()
        return self._sort(name)


class SamSplitter(SamBase):
    _operation_dir = "alignment/5_split"

    def _get_bam_ref(self, ref: bytes):
        return os.path.join(self.output_dir, f"{ref.decode()}{BAM_EXT}")
    
    @property
    def bams_out(self):
        return list(map(self._get_bam_ref, self.refs))
    
    def _output_bam_ref(self, ref: bytes):
        output = self._get_bam_ref(ref)
        cmd = [SAMTOOLS_CMD, "view", "-b", "-o", output,
               self.xam_in, ref.decode()]
        run_cmd(cmd)
        return output
    
    def _split_bam(self):
        self.index_bam_in()
        return list(map(self._output_bam_ref, self.refs))
    
    def run(self):
        print(f"\nSplitting {self.xam_in} into Individual References\n")
        self._make_output_dir()
        return self._split_bam()


class SamOutputter(SamBase):
    _operation_dir = "alignment"

    @property
    def output_dir(self):
        return self._get_dir(OUTPUT_DIR)
    
    def run(self):
        print(f"\nOutputting Cleaned BAM files to {self.output_dir}\n")
        self._make_output_dir()
        output = switch_directory(self.xam_in, self.output_dir)
        os.rename(self.xam_in, output)
        return output


'''
primer1 = "CAGCACTCAGAGCTAATACGACTCACTATA"
primer1rc = "TATAGTGAGTCGTATTAGCTCTGAGTGCTG"
primer2 = "TGAAGAGCTGGAACGCTTCACTGA"
primer2rc = "TCAGTGAAGCGTTCCAGCTCTTCA"
adapters5 = (primer1, primer2rc)
adapters3 = (primer2, primer1rc)
'''
