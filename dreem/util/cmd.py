import shlex
import subprocess

from typing import Any, List


# Commands for external applications

BOWTIE2_CMD = "bowtie2"
BOWTIE2_BUILD_CMD = "bowtie2-build"
CUTADAPT_CMD = "cutadapt"
FASTQC_CMD = "fastqc"
SAMTOOLS_CMD = "samtools"


# Command utility functions

def run_cmd(args: List[Any], shell: bool = False):
    args_str = tuple(map(str, args))
    subprocess.run(args_str, check=True, shell=shell)
    cmd = shlex.join(args_str)
    return cmd
