import json
import os
import random
import shlex
import subprocess
import sys
from tempfile import NamedTemporaryFile
from typing import List, Set, Dict
import yaml

import numpy as np
import pandas as pd
import pyarrow.orc  # This prevents: AttributeError: module 'pyarrow' has no attribute 'orc'

from dreem.aggregate import poisson


# CONSTANTS
BASES = b"ACGT"
COMPS = b"TGCA"
RBASE = b"ACGU"
RCOMP = b"UGCA"
BASEN = b"N"
BASES_SET = set(BASES)
BASEN_SET = set(BASES + BASEN)

# BYTE ENCODINGS FOR MUTATION VECTORS
BLANK = b"\x00"  # 00000000 (000): no coverage at this position
MATCH = b"\x01"  # 00000001 (001): match with reference
DELET = b"\x02"  # 00000010 (002): deletion from reference
INS_5 = b"\x04"  # 00000100 (004): insertion 5' of base in reference
INS_3 = b"\x08"  # 00001000 (008): insertion 3' of base in reference
SUB_A = b"\x10"  # 00010000 (016): substitution to A
SUB_C = b"\x20"  # 00100000 (032): substitution to C
SUB_G = b"\x40"  # 01000000 (064): substitution to G
SUB_T = b"\x80"  # 10000000 (128): substitution to T
AMBIG = b"\xff"  # 11111111 (255): could be anything
SUB_N = (SUB_A[0] | SUB_C[0] | SUB_G[0] | SUB_T[0]).to_bytes()
ANY_N = (MATCH[0] | SUB_N[0]).to_bytes()
MADEL = (MATCH[0] | DELET[0]).to_bytes()
NOSUB = (AMBIG[0] ^ SUB_N[0]).to_bytes()
UNAMB_SET = set(b"".join([BLANK, MATCH, DELET, INS_5, INS_3,
                          SUB_A, SUB_C, SUB_G, SUB_T]))

# COMMANDS
BOWTIE2_CMD = "bowtie2"
BOWTIE2_BUILD_CMD = "bowtie2-build"
CUTADAPT_CMD = "cutadapt"
FASTQC_CMD = "fastqc"
PASTE_CMD = "paste"
SAMTOOLS_CMD = "samtools"

# GLOBAL SETTINGS
DEFAULT_PROCESSES = cpus if (cpus := os.cpu_count()) else 1
BASE_COLORS = {"A": "#D3822A", "C": "#5AB4E5", "G": "#EAE958", "T": "#357766"}
PHRED_ENCODING = 33
OUTPUT_DIR = "output"
TEMP_DIR = "temp"


def make_folder(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)
    return folder

def make_cmd(args, module):
    cmd = 'dreem-' + module + ' '
    for key, value in args.items():
        if type(value) in (list, tuple):
            for v in value:
                cmd += '--' + key + ' ' + str(v) + ' '
        elif value is not None:
            cmd += '--' + key + ' ' + str(value) + ' '
    return cmd

def run_cmd(args: List[str], shell: bool = False):
    subprocess.run(args, check=True, shell=shell)
    cmd = shlex.join(args)
    return cmd

def name_temp_file(dirname: str, prefix: str, suffix: str):
    file = NamedTemporaryFile(dir=dirname, prefix=prefix, suffix=suffix,
                              mode="w", delete=False)
    file.close()
    return file.name

def get_filename(file_path):
    return os.path.splitext(os.path.basename(file_path))[0]

def switch_directory(old_path, new_dir):
    return os.path.join(new_dir, os.path.basename(old_path))

def try_remove(file):
    try:
        os.remove(file)
    except OSError:
        pass

def get_files(folder, ext):
    paths = []
    for construct in os.listdir(folder):
        if os.path.isdir(os.path.join(folder, construct)):
            for section in os.listdir(os.path.join(folder, construct)):
                if not section.endswith(ext):
                    continue
                paths.append(os.path.join(folder, construct, section))
    else:
        if construct.endswith(ext):
            paths.append(os.path.join(folder, construct))
    return paths


class Seq(bytes):
    __slots__ = []

    alph = b""
    comp = b""
    alphaset = set(alph)
    trans = alph.maketrans(alph, comp)

    def __init__(self, seq: bytes):
        self.validate_seq(seq)
        super().__init__()
    
    @classmethod
    def validate_seq(cls, seq):
        if not seq:
            raise ValueError("seq is empty")
        if set(seq) - cls.alphaset:
            raise ValueError(f"Invalid characters in seq: '{seq.decode()}'")

    @property
    def rc(self):
        return self.__class__(self[::-1].translate(self.trans))

    def __getitem__(self, item):
        return self.__class__(super().__getitem__(item))

    def __str__(self):
        return self.decode()


class DNA(Seq):
    alph = BASES
    comp = COMPS
    alphaset = set(alph)
    trans = alph.maketrans(alph, comp)

    def tr(self):
        """
        Transcribe DNA into RNA.
        """
        return RNA(self.replace(b"T", b"U"))


class RNA(Seq):
    alph = RBASE
    comp = RCOMP
    alphaset = set(alph)
    trans = alph.maketrans(alph, comp)

    def rt(self):
        """
        Reverse transcribe RNA into DNA.
        """
        return DNA(self.replace(b"U", b"T"))


class AmbigDNA(DNA):
    alph = BASES + BASEN
    comp = COMPS + BASEN
    alphaset = set(alph)
    trans = alph.maketrans(alph, comp)


class FastaIO(object):
    __slots__ = ["_path", "_refs"]

    defsymbol = b">"
    deftrunc = len(defsymbol)

    def __init__(self, path: str):
        self._path = path



class FastaParser(FastaIO):
    def __init__(self, path: str):
        super().__init__(path)
        self._refs: Set[bytes] = set()

    @classmethod
    def _parse_fasta_record(cls, fasta, line: bytes):
        if not line.startswith(cls.defsymbol):
            raise ValueError("FASTA definition line does not start with "
                             f"'{cls.defsymbol.decode()}'")
        name = line.split()[0][cls.deftrunc:]
        seq = bytearray()
        while (line := fasta.readline()) and not line.startswith(cls.defsymbol):
            seq.extend(line.rstrip())
        seq = DNA(bytes(seq))
        return line, name, seq

    def parse(self):
        with open(self._path, "rb") as f:
            line = f.readline()
            while line:
                line, name, seq = self._parse_fasta_record(f, line)
                if name in self._refs:
                    raise ValueError(
                        f"Duplicate entry in {self._path}: '{name.decode()}'")
                self._refs.add(name)
                yield name, seq


class FastaWriter(FastaIO):
    def __init__(self, path: str, refs: Dict[bytes, DNA]):
        super().__init__(path)
        self._refs = refs
    
    def write(self):
        with open(self._path, "wb") as f:
            for ref, seq in self._refs.items():
                f.write(b"".join(self.defsymbol, ref, b"\n", seq, b"\n"))


def fastq_to_df(fastq_file):    
    df, data = pd.DataFrame(), pd.read_csv(fastq_file, sep='\t', header=None)
    df['construct'] = data.iloc[np.arange(0,len(data),4)].reset_index(drop=True)
    df['sequence'] = data.iloc[np.arange(1,len(data),4)].reset_index(drop=True)
    df['quality'] = data.iloc[np.arange(3,len(data),4)].reset_index(drop=True)
    return df


def sam_to_df(path, skiprows=None):
    with open(path) as f:
        lines = f.readlines()
    lines = [l for l in lines if not l.startswith('@')]
    lines = [l.split('\t') for l in lines]
    df = pd.read_csv(path, sep='\t', skiprows=skiprows, header=None)
    df = df[df.columns[:12]]
    df.columns = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT','TLEN', 'SEQ', 'QUAL','OPT']
    return df


def sort_fastq_pairs(fq1s, fq2s):
    assert len(fq1s) >= len(fq2s), 'More fq2s than fq1s'
    for f2 in fq2s:
        if f2.replace('_R2', '_R1') not in fq1s:
            raise ValueError(f'No matching pair for {f2}')
    fq1s = [f2.replace('_R2', '_R1') for f2 in fq2s] + [f for f in fq1s if f not in fq2s]
    fq2s += [None for f in fq1s if f not in fq2s]
    samples = [f.split('/')[-1].split('.')[0].split('_')[:-1] for f in fq1s]
    samples = ['_'.join(s) for s in samples]
    return fq1s, fq2s, samples
