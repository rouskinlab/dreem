import itertools
from multiprocessing import Pool
import os
import re
from typing import List, Optional, Tuple

from dreem.util.util import FASTQC_CMD, CUTADAPT_CMD, BOWTIE2_CMD, PASTE_CMD, BASEN, \
    BOWTIE2_BUILD_CMD, TEMP_DIR, OUTPUT_DIR, DEFAULT_PROCESSES, run_cmd, \
    try_remove, switch_directory, PHRED_ENCODING, SAMTOOLS_CMD



# General parameters
DEFAULT_INTERLEAVED = False
SAM_HEADER = b"@"
SAM_ALIGN_SCORE = b"AS:i:"
SAM_EXTRA_SCORE = b"XS:i:"

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
DEFAULT_TRIM_CORES = cpus if (cpus := os.cpu_count()) else 1
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
DEFAULT_MAP_QUAL_MIN = 30
DEFAULT_N_CEIL = "L,0,0.05"
DEFAULT_SEED = 12
DEFAULT_EXTENSIONS = 6
DEFAULT_RESEED = 1
DEFAULT_PADDING = 4
DEFAULT_ALIGN_THREADS = os.cpu_count()
DEFAULT_METRICS = 60
MATCH_BONUS = "1"
MISMATCH_PENALTY = "-1,-1"
N_PENALTY = "0"
REF_GAP_PENALTY = "0,-1"
READ_GAP_PENALTY = "0,-1"
IGNORE_QUALS = True
SAM_EXT = ".sam"
BAM_EXT = ".bam"


def get_fastq_name(fastq: str, fastq2: Optional[str] = None):
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
        if len(name) != len(name2):
            raise ValueError("FASTQ names were different lengths: "
                             f"'{name}' and '{name2}'")
        diffs = [i for i, (n1, n2) in enumerate(zip(name, name2)) if n1 != n2]
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


class SeqFileBase(object):
    _operation_dir = ""

    def __init__(self,
                 root_dir: str,
                 ref_file: str,
                 paired: bool,
                 encoding: int = PHRED_ENCODING) -> None:
        self._root_dir = root_dir
        self._ref_file = ref_file
        self._paired = paired
        self._encoding = encoding
    
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
    def paired(self):
        return self._paired
    
    @property
    def encoding(self):
        return self._encoding
    
    @property
    def encoding_arg(self):
        return f"--phred{self.encoding}"
    
    def run(output_dir: str, **kwargs):
        raise NotImplementedError


class FastqBase(SeqFileBase):
    def __init__(self,
                 root_dir: str,
                 ref_file: str,
                 samples: List[str],
                 fastqs: List[str],
                 paired: bool,
                 encoding: int = PHRED_ENCODING) -> None:
        super().__init__(root_dir, ref_file, paired, encoding=encoding)
        if len(fastqs) != len(samples):
            raise ValueError("fastqs and samples had different lengths.")
        self._samples = samples
        self._fastqs = fastqs

    @property
    def samples(self):
        return self._samples
    
    @property
    def fastqs(self):
        return self._fastqs
        
    def _get_temp_dir(self, sample: str):
        return os.path.join(self.root_dir, TEMP_DIR, self.operation_dir, sample)
    
    @property
    def temp_dirs(self):
        return list(map(self._get_temp_dir, self.samples))
    
    @property
    def temp_outputs(self):
        return list(map(switch_directory, self.fastqs, self.temp_dirs))
    
    def qc(self, extract: bool = DEFAULT_EXTRACT):
        cmd = [FASTQC_CMD]
        if extract:
            cmd.append("--extract")
        cmd.extend(self.fastqs)
        run_cmd(cmd)


class FastqInterleaver(FastqBase):
    _operation_dir = "align/interleave"

    def __init__(self,
                 root_dir: str,
                 ref_file: str,
                 samples: List[str],
                 fastq1s: List[str],
                 fastq2s: List[str],
                 encoding: int = PHRED_ENCODING) -> None:
        if len(fastq1s) != len(fastq2s):
            raise ValueError("fastq1s and fastq2s had different lengths.")
        super().__init__(root_dir, ref_file, samples, fastq1s, paired=True,
                         encoding=encoding)
        self._fastq2s = fastq2s
    
    @property
    def fastq1s(self):
        return self._fastqs
    
    @property
    def fastq2s(self):
        return self._fastq2s
    
    @property
    def fastqs(self):
        return [f"{sample}.fq" for sample in self.samples]

    def run(self):
        outputs = self.temp_outputs
        for fq1, fq2, output in zip(self.fastq1s, self.fastq2s, outputs):
            cmd = [PASTE_CMD, fq1, fq2, ">", output]
            run_cmd(cmd)
        return outputs


class FastqMasker(FastqBase):
    _operation_dir = "align/mask"

    def _mask(self, fq_in: str, fq_out: str, min_qual: int):
        min_code = min_qual + self.encoding
        NL = b"\n"[0]
        BN = BASEN[0]
        with open(fq_in, "rb") as fqi, open(fq_out, "wb") as fqo:
            for seq_header in fqi:
                masked = bytearray(seq_header)
                seq = fqi.readline()
                qual_header = fqi.readline()
                quals = fqi.readline()
                if len(seq) != len(quals):
                    raise ValueError("seq and qual have different lengths")
                masked.extend(base if qual >= min_code or base == NL else BN
                              for base, qual in zip(seq, quals))
                masked.extend(qual_header)
                masked.extend(quals)
                fqo.write(masked)

    def run(self, min_qual: int = DEFAULT_MIN_BASE_QUALITY,
            parallel: bool = False):
        args = zip(self.fastqs, outputs := self.temp_outputs,
                   [min_qual] * len(self.fastqs))
        if parallel:
            with Pool(DEFAULT_PROCESSES) as pool:
                pool.starmap(self._mask, args)
        else:
            for arg in args:
                self._mask(*arg)
        return outputs


class FastqTrimmer(FastqBase):
    _operation_dir = "align/trim"

    def _cutadapt(self,
                  fq_in: str,
                  fq_out: str,
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
                  cores=DEFAULT_TRIM_CORES):
        cmd = [CUTADAPT_CMD]
        if cores >= 0:
            cmd.append(f"--cores {cores}")
        if nextseq:
            nextseq_qual = qual1 if qual1 else DEFAULT_MIN_BASE_QUALITY
            cmd.append(f"--nextseq-trim {nextseq_qual}")
        else:
            if qual1 is not None:
                cmd.append(f"-q {qual1}")
            if qual2 is not None:
                self.fastq2
                cmd.append(f"-Q {qual2}")
        adapters = {"g": adapters15, "a": adapters13,
                    "G": adapters25, "A": adapters23}
        for arg, adapter in adapters.items():
            if adapter and (arg.islower() or self.paired):
                if isinstance(adapter, str):
                    adapter = (adapter,)
                if not isinstance(adapter, tuple):
                    raise ValueError("adapters must be str or tuple")
                for adapt in adapter:
                    cmd.append(f"-{arg} {adapt}")
        if min_overlap >= 0:
            cmd.append(f"-O {min_overlap}")
        if max_error >= 0:
            cmd.append(f"-e {max_error}")
        if not indels:
            cmd.append("--no-indels")
        if discard_trimmed:
            cmd.append("--discard-trimmed")
        if discard_untrimmed:
            cmd.append("--discard-untrimmed")
        if min_length:
            cmd.append(f"-m {min_length}")
        cmd.append("--interleaved")
        cmd.append(f"-o {fq_out}")
        cmd.append(fq_in)
        run_cmd(cmd)
    
    def run(self, parallel: bool = False, **kwargs):
        args = zip(self.fastqs, outputs := self.temp_outputs)
        if parallel:
            raise NotImplementedError()
        else:
            for fq_in, fq_out in args:
                self._cutadapt(fq_in, fq_out, **kwargs)
        return outputs


class SamBase(SeqFileBase):
    def __init__(self,
                 root_dir: str,
                 ref_file: str,
                 
                 paired: bool,
                 encoding: int = PHRED_ENCODING) -> None:
        super().__init__(root_dir, ref_file, paired, encoding)
        


class FastqAligner(FastqBase):
    operation_dir = "align/bowtie2"

    def _bowtie2_build(self):
        """
        Build an index of a reference genome using Bowtie 2.
        :param ref: (str) path to the reference genome FASTA file
        :return: None
        """
        cmd = [BOWTIE2_BUILD_CMD, self.ref_file, self.ref_prefix]
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
                 map_qual_min=DEFAULT_MAP_QUAL_MIN,
                 n_ceil=DEFAULT_N_CEIL,
                 seed=DEFAULT_SEED,
                 extensions=DEFAULT_EXTENSIONS,
                 reseed=DEFAULT_RESEED,
                 padding=DEFAULT_PADDING,
                 threads=DEFAULT_ALIGN_THREADS,
                 metrics=DEFAULT_METRICS):
        cmd = [BOWTIE2_CMD]
        fastq_flag = "--interleaved" if self.paired else "-U"
        cmd.extend([fastq_flag, ",".join(self.fastqs)])
        cmd.append(f"-x {self.ref_prefix}")
        cmd.append(f"-S {self.sam_temp}")
        cmd.append(f"{self.encoding_arg}")
        cmd.append("--xeq")
        cmd.append(f"--ma {MATCH_BONUS}")
        cmd.append(f"--mp {MISMATCH_PENALTY}")
        cmd.append(f"--np {N_PENALTY}")
        cmd.append(f"--rfg {REF_GAP_PENALTY}")
        cmd.append(f"--rdg {READ_GAP_PENALTY}")
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
            cmd.append(f"-I {frag_len_min}")
        if frag_len_max:
            cmd.append(f"-X {frag_len_max}")
        if map_qual_min:
            cmd.append(f"--score-min C,{map_qual_min}")
        if n_ceil:
            cmd.append(f"--n-ceil {n_ceil}")
        if seed:
            cmd.append(f"-L {seed}")
        if extensions:
            cmd.append(f"-D {extensions}")
        if reseed:
            cmd.append(f"-R {reseed}")
        if padding:
            cmd.append(f"--dpad {padding}")
        if threads:
            cmd.append(f"-p {threads}")
        if metrics:
            cmd.append(f"--met-stderr --met {metrics}")
        if IGNORE_QUALS:
            cmd.append("--ignore-quals")
        run_cmd(cmd)
    
    def run(self, **kwargs):
        self._bowtie2_build()
        self._bowtie2(**kwargs)


class AlignmentCleaner(SeqFileBase):
    operation_dir = "align/cleaned"

    def remove_equal_mappers(self):
        print("\n\nRemoving reads mapping equally to multiple locations")
        pattern_a = re.compile(SAM_ALIGN_SCORE + rb"(\d+)")
        pattern_x = re.compile(SAM_EXTRA_SCORE + rb"(\d+)")

        def get_score(line, ptn):
            return (float(match.groups()[0])
                    if (match := ptn.search(line)) else None)

        def is_best_alignment(line):
            return ((score_x := get_score(line, pattern_x)) is None
                    or score_x < get_score(line, pattern_a))

        kept = 0
        removed = 0
        with open(self.sam_input, "rb") as sami, open(self.sam_temp, "wb") as samo:
            # Copy the header from the input to the output SAM file.
            while (line := sami.readline()).startswith(SAM_HEADER):
                samo.write(line)
            if self.paired:
                for line2 in sami:
                    if is_best_alignment(line) or is_best_alignment(line2):
                        samo.write(line)
                        samo.write(line2)
                        kept += 1
                    else:
                        removed += 1
                    line = sami.readline()
            else:
                while line:
                    if is_best_alignment(line):
                        samo.write(line)
                        kept += 1
                    else:
                        removed += 1
                    line = sami.readline()
        total = kept + removed
        items = "Pairs" if self.paired else "Reads"
        print(f"\n{items} processed: {total}")
        print(f"{items} kept: {kept} ({round(100 * kept / total, 2)}%)")
        print(f"{items} lost: {removed} ({round(100 * removed / total, 2)}%)")
        return kept, removed


class AlignmentFinisher(SeqFileBase):
    operation_dir = "align"

    def sort(self):
        cmd = [SAMTOOLS_CMD, "sort", "-o", self.bam_out, self.sam_input]
        run_cmd(cmd)
    
    def index(self):
        cmd = [SAMTOOLS_CMD, "index", "-b", self.bam_out]
        run_cmd(cmd)


def run(output_dir, ref, fastq1, fastq2=None, **kwargs):
    # Check the file extensions.
    primer1 = "CAGCACTCAGAGCTAATACGACTCACTATA"
    primer1rc = "TATAGTGAGTCGTATTAGCTCTGAGTGCTG"
    primer2 = "TGAAGAGCTGGAACGCTTCACTGA"
    primer2rc = "TCAGTGAAGCGTTCCAGCTCTTCA"
    adapters5 = (primer1, primer2rc)
    adapters3 = (primer2, primer1rc)
    primer_trimmed_reads = os.path.join(output_dir,
        f"__{sample}_primer_trimmed.fastq")
    trim(adapter_trimmed_reads, primer_trimmed_reads, interleaved=paired,
         adapters15=adapters5, adapters13=adapters3,
         adapters25=adapters5, adapters23=adapters3,
         min_overlap=6, min_length=10, discard_trimmed=True)
    try_remove(adapter_trimmed_reads)
    # Mask low-quality bases remaining in the trimmed file.
    masked_trimmed_reads = os.path.join(output_dir,
        f"__{sample}_masked_trimmed.fastq")
    mask_low_qual(primer_trimmed_reads, masked_trimmed_reads)
    try_remove(primer_trimmed_reads)
    # Quality check the trimmed FASTQ file.
    #fastqc(output_dir, masked_trimmed_reads)
    # Index the reference.
    prefix = index(ref)
    # Align the reads to the reference.
    sam_align = os.path.join(output_dir, f"__{sample}_align.sam")
    bowtie(sam_align, prefix, masked_trimmed_reads, interleaved=True)
    try_remove(masked_trimmed_reads)
    # Delete reads aligning equally well to multiple locations.
    sam_pruned = os.path.join(output_dir, f"__{sample}_pruned.sam")
    remove_equal_mappers(sam_align, sam_pruned)
    # Delete the non-pruned SAM file.
    try_remove(sam_align)
    # Sort the SAM file while converting to BAM format
    bam_file = os.path.join(output_dir, f"{sample}.bam")
    sort_sam(sam_pruned, bam_file)
    # Delete the SAM file.
    try_remove(sam_pruned)
    # Index the BAM file.
    index_bam(bam_file)
