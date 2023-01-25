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

BLANK_INT = 0
MATCH_INT = 1
DELET_INT = 2
INS_5_INT = 4
INS_3_INT = 8
SUB_A_INT = 16
SUB_C_INT = 32
SUB_G_INT = 64
SUB_T_INT = 128
SUB_N_INT = SUB_A_INT + SUB_C_INT + SUB_G_INT + SUB_T_INT
AMBIG_INT = 255

B_MATCH = 0
B_DELET = 1
B_INS_5 = 2
B_INS_3 = 3
B_SUB_A = 4
B_SUB_C = 5
B_SUB_G = 6
B_SUB_T = 7

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
DEFAULT_PHRED_ENCODING = 33
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
    os.system(' '.join(args))#, check=True, shell=shell)
    cmd = shlex.join(args)
    return cmd

def name_temp_file(dirname: str, prefix: str, suffix: str):
    file = NamedTemporaryFile(dir=dirname, prefix=prefix, suffix=suffix,
                              mode="w", delete=False)
    file.close()
    return file.name

def get_filename(file_path: str):
    return os.path.splitext(os.path.basename(file_path))[0]

def switch_directory(old_path: str, new_dir: str):
    return os.path.join(new_dir, os.path.basename(old_path))

def try_remove(file: str):
    try:
        os.remove(file)
    except OSError:
        pass

def try_rmdir(dir: str):
    try:
        os.rmdir(dir)
    except OSError:
        pass

def get_files(folder: str, ext: str):
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


def sam_to_df(path):
    skiprows = 0
    with open(path) as f:
        while f.readline().startswith('@'):
            skiprows += 1
    df = pd.read_csv(path, sep='\t', skiprows=skiprows, header=None)
    df = df[df.columns[:11]]
    df.columns = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT','TLEN', 'SEQ', 'QUAL']
    return df


def sort_fastq_pairs(fq1s, fq2s):
    if type(fq1s) is str:
        fq1s = [fq1s]
    if type(fq2s) is str:
        fq2s = [fq2s]
    assert len(fq1s) >= len(fq2s), 'More fq2s than fq1s'
    for f2 in fq2s:
        if f2.replace('_R2', '_R1') not in fq1s:
            raise ValueError(f'No matching pair for {f2}')
    fq1s = [f2.replace('_R2', '_R1') for f2 in fq2s] + [f for f in fq1s if f.replace('_R1','_R2') not in fq2s]
    fq2s += [None for f in fq1s if f.replace('_R1','_R2') not in fq2s]
    samples = [f.split('/')[-1].split('.')[0].split('_')[:-1] for f in fq1s]
    samples = ['_'.join(s) for s in samples]
    return fq1s, fq2s, samples


def query_muts(muts: np.ndarray, bits: int, sum_up = True, axis=0, set_type = 'superset'):
    """
    Count the number of times a query mutation occurs in each column
    or one column of a set of mutation vectors.
    The counting operation comprises three steps:
    1. bitwise AND to confirm at least one "1" bit is shared, e.g.
       bits: 11110000 & muts: 00100000 -> 00100000 (True)
       bits: 11110000 & muts: 00100010 -> 00100000 (True)
       bits: 11110000 & muts: 00000000 -> 00000000 (False)
    2. bitwise OR to confirm no "1" bit in muts is not in bits, e.g.
       bits: 11110000 | muts: 00100000 -> 11110000 =? 11110000 (True)
       bits: 11110000 | muts: 00100010 -> 11110010 =? 11110000 (False)
       bits: 11110000 | muts: 00000000 -> 11110000 =? 11110000 (True)
    3. logical AND to confirm that both tests pass, e.g.
       bits: 11110000, muts: 00100000 -> True  AND True  (True)
       bits: 11110000, muts: 00100010 -> True  AND False (False)
       bits: 11110000, muts: 00000000 -> False AND True  (False)

    Arguments
    muts: NDArray of a set of mutation vectors (2-dimensional)
          or one column in a set of mutation vectors (1-dimensional).
          Data type must be uint8.
    bits: One-byte int in the range [0, 256) representing the mutation
          to be queried. The bits in the int encode the mutation as
          defined above, e.g.
          - 00000010 (int 2) is a deletion
          - 11010001 (int 209) is either substitution to A, G, or T
                               or a match to C
    
    Returns
    if sum_up: 
        int of the number of times the query mutation occurs in muts
        count: If muts is 1-dimensional, int of the number of times the
            query mutation occurs in muts.
            If muts is 2-dimensional, NDArray with one int for each
            column in muts.
    if not sum_up:
        bool NDArray with one bool for each column in muts.
        True if the query mutation occurs in the column, False otherwise.
    """
    if not muts.dtype == np.uint8:
        raise TypeError('muts must be of type uint8 and not {}'.format(muts.dtype))
    #print(np.logical_and(muts & bits, (muts | bits) == bits).sum(axis=0))
    assert isinstance(bits, int) and 0 <= bits < 256

        
    if set_type == 'subset':
        logic_fun = lambda m, b: np.logical_and(m & b, (m | b) == b)
    if set_type == 'superset':
        logic_fun = lambda m, b: np.array(m & b, dtype=bool)
    
    if sum_up:
        return logic_fun(muts, bits).sum(axis=axis)
    else:
        return logic_fun(muts, bits)
     
