import shlex
import subprocess

from typing import List


# Commands for external applications

BOWTIE2_CMD = "bowtie2"
BOWTIE2_BUILD_CMD = "bowtie2-build"
CUTADAPT_CMD = "cutadapt"
FASTQC_CMD = "fastqc"
SAMTOOLS_CMD = "samtools"


# Command utility functions

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
