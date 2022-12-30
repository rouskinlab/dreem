from collections import defaultdict
from multiprocessing import Pool
import os

from typing import List, Optional, Tuple

from dreem.util.cmd import FASTQC_CMD, CUTADAPT_CMD, BOWTIE2_CMD, \
    BOWTIE2_BUILD_CMD, run_cmd
from dreem.util.dflt import NUM_PROCESSES, PHRED_ENCODING
from dreem.util.ngs import NgsFileBase
from dreem.util.seq import get_diffs
from dreem.util.cli import DEFAULT_MIN_BASE_QUALITY, DEFAULT_ILLUMINA_ADAPTER, DEFAULT_MIN_OVERLAP, DEFAULT_MAX_ERROR, DEFAULT_INDELS, DEFAULT_NEXTSEQ_TRIM, DEFAULT_DISCARD_TRIMMED, DEFAULT_DISCARD_UNTRIMMED, DEFAULT_MIN_LENGTH, DEFAULT_SCORE_MIN
from dreem.util.cli import DEFAULT_LOCAL, DEFAULT_UNALIGNED, DEFAULT_DISCORDANT, DEFAULT_MIXED, DEFAULT_DOVETAIL, DEFAULT_CONTAIN, DEFAULT_FRAG_LEN_MIN, DEFAULT_FRAG_LEN_MAX, DEFAULT_N_CEILING, DEFAULT_SEED_INTERVAL, DEFAULT_GAP_BAR, DEFAULT_SEED_SIZE, DEFAULT_EXTENSIONS, DEFAULT_RESEED, DEFAULT_PADDING, DEFAULT_ALIGN_THREADS, MATCH_BONUS, MISMATCH_PENALTY, N_PENALTY, REF_GAP_PENALTY, READ_GAP_PENALTY, IGNORE_QUALS
from dreem.util.cli import DEFAULT_INTERLEAVE_OUTPUT, FASTQ2


# FastQC parameters
DEFAULT_EXTRACT = False


def get_fastq_name(fastq: str, fastq2: str=""):
    exts = (".fq", ".fastq")
    reads = ("mate", "r")
    base = os.path.basename(fastq)
    counts = {ext: count for ext in exts if (count := base.count(ext))}
    if not counts:
        raise ValueError(f"File '{fastq}' had no FASTQ extension.")
    if sum(counts.values()) > 1:
        raise ValueError(f"File '{fastq}' had multiple FASTQ extensions.")
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


def get_fastq_dir(fastq: str, fastq2: Optional[str]=None):
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
                raise ValueError(f"Found ≥2 mates for {file1} in {fq_dir}")
            fq_files.pop(file2)
            pair = (os.path.join(fq_dir, file1), os.path.join(fq_dir, file2))
            pairs[get_fastq_name(file1, file2)] = pair
    return pairs


class FastqBase(NgsFileBase):
    def __init__(self,
                 root_dir: str,
                 ref_file: str,
                 sample: str,
                 paired: bool,
                 fastq: str,
                 fastq2: str=FASTQ2,
                 interleave_out: bool=DEFAULT_INTERLEAVE_OUTPUT,
                 encoding: int = PHRED_ENCODING) -> None:
        if fastq2 and not paired:
            raise ValueError("fastq2 can only be given if reads are paired")
        super().__init__(root_dir, ref_file, sample)
        self._fq_in = fastq
        self._fq2_in = fastq2
        self._paired = paired
        self._interleave_out = interleave_out
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
            fqs = [f"{self.name}_R{r}" for r in ("1", "2")]
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
    
    def qc_inputs(self, extract: bool=DEFAULT_EXTRACT):
        return self._qc(self.inputs, extract)

    def qc_outputs(self, extract: bool=DEFAULT_EXTRACT):
        return self._qc(self.outputs, extract)


class FastqTrimmer(FastqBase):
    _operation_dir = "alignment/1_trim"

    def _cutadapt(self,
                  qual1: int=DEFAULT_MIN_BASE_QUALITY,
                  qual2: int=0,
                  adapters15: Tuple[str]=(),
                  adapters13: Tuple[str]=(DEFAULT_ILLUMINA_ADAPTER,),
                  adapters25: Tuple[str]=(),
                  adapters23: Tuple[str]=(DEFAULT_ILLUMINA_ADAPTER,),
                  min_overlap: int=DEFAULT_MIN_OVERLAP,
                  max_error: float=DEFAULT_MAX_ERROR,
                  indels: bool=DEFAULT_INDELS,
                  nextseq_trim: bool=DEFAULT_NEXTSEQ_TRIM,
                  discard_trimmed: bool=DEFAULT_DISCARD_TRIMMED,
                  discard_untrimmed: bool=DEFAULT_DISCARD_UNTRIMMED,
                  min_length: bool=DEFAULT_MIN_LENGTH,
                  cores: int=NUM_PROCESSES):
        cmd = [CUTADAPT_CMD]
        if cores >= 0:
            cmd.extend(["--cores", str(cores)])
        if nextseq_trim:
            if qual1 > 0:
                cmd.extend(["--nextseq-trim", str(qual1)])
        else:
            if qual1 > 0:
                cmd.extend(["-q", str(qual1)])
            if qual2 > 0:
                cmd.extend(["-Q", str(qual2)])
        adapters = {"g": adapters15, "a": adapters13,
                    "G": adapters25, "A": adapters23}
        for arg, adapter in adapters.items():
            if adapter and (self.paired or arg.islower()):
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
    
    def run(self, **kwargs):
        print(f"\nTrimming Adapters from {self.fq_in}\n")
        self._make_output_dir()
        self._cutadapt(**kwargs)


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
                 score_min=DEFAULT_SCORE_MIN,
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
        if score_min:
            cmd.extend(["--score-min", score_min])
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
    
    def run(self, **kwargs):
        print(f"\nAligning Reads {self.fq_in} to Reference {self.ref_file}\n")
        self._make_output_dir()
        self._bowtie2_build()
        self._bowtie2(**kwargs)


'''
primer1 = "CAGCACTCAGAGCTAATACGACTCACTATA"
primer1rc = "TATAGTGAGTCGTATTAGCTCTGAGTGCTG"
primer2 = "TGAAGAGCTGGAACGCTTCACTGA"
primer2rc = "TCAGTGAAGCGTTCCAGCTCTTCA"
adapters5 = (primer1, primer2rc)
adapters3 = (primer2, primer1rc)
'''
