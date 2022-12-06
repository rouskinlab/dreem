import yaml, sys, os, random
import pandas as pd
import numpy as np
import subprocess
import json
import pyarrow.orc  # This prevents: AttributeError: module 'pyarrow' has no attribute 'orc'
from dreem.aggregate import poisson

def make_folder(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)
    return folder

def clear_folder(folder):
    os.system('rm -fr ' + folder)
    os.makedirs(folder)

def run_cmd(cmd):
    return os.system(cmd), cmd

def make_cmd(args, module):
    cmd = 'dreem-' + module + ' '
    for key, value in args.items():
        if type(value) in (list, tuple):
            for v in value:
                cmd += '--' + key + ' ' + str(v) + ' '
        elif value is not None:
            cmd += '--' + key + ' ' + str(value) + ' '
    return cmd

class Seq(bytes):
    __slots__ = []

    alph = b""
    comp = b""
    low_qual = b"N"

    def __init__(self, seq: bytes):
        if not seq:
            raise ValueError("seq is empty")
        if any(base not in self.alph for base in seq):
            raise ValueError(f"Invalid characters in seq: '{seq.decode()}'")
        super().__init__()

    @property
    def rc(self):
        return self.__class__(
            self[::-1].translate(self.maketrans(self.alph, self.comp)))

    def __getitem__(self, item):
        return self.__class__(super().__getitem__(item))

    def __str__(self):
        return self.decode()


class DNA(Seq):
    alph = b"ACGT"
    comp = b"TGCA"

    @property
    def tr(self):
        return RNA(self.replace(b"T", b"U"))


class RNA(Seq):
    alph = b"ACGU"
    comp = b"UGCA"

    @property
    def rt(self):
        return DNA(self.replace(b"U", b"T"))


class Primer(str):
    @property
    def as_dna(self):
        return DNA(self.encode())




def fastq_to_df(fastq_file):    
    df, data = pd.DataFrame(), pd.read_csv(fastq_file, sep='\t', header=None)
    df['construct'] = data.iloc[np.arange(0,len(data),4)].reset_index(drop=True)
    df['sequence'] = data.iloc[np.arange(1,len(data),4)].reset_index(drop=True)
    df['quality'] = data.iloc[np.arange(3,len(data),4)].reset_index(drop=True)
    return df


def sam_to_df(path):
    with open(path) as f:
        lines = f.readlines()
    lines = [l for l in lines if not l.startswith('@')]
    lines = [l.split('\t') for l in lines]
    df = pd.DataFrame(lines)
    df.columns = ['QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT', 'TLEN', 'SEQ']
    return df


