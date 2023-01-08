import shlex
import subprocess,os

from typing import List


# Commands for external applications

BOWTIE2_CMD = "bowtie2"
BOWTIE2_BUILD_CMD = "bowtie2-build"
CUTADAPT_CMD = "cutadapt"
FASTQC_CMD = "fastqc"
SAMTOOLS_CMD = "samtools"


# Command utility functions

def run_cmd(args: List[str], shell: bool = False):
    os.system(' '.join(args))#, check=True, shell=shell)
    cmd = shlex.join(args)
    return cmd
