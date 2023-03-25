import logging
import shlex
import subprocess
from typing import Any, List

logger = logging.getLogger(__name__)


# Commands for external applications

BOWTIE2_CMD = "bowtie2"
BOWTIE2_BUILD_CMD = "bowtie2-build"
CUTADAPT_CMD = "cutadapt"
FASTQC_CMD = "fastqc"
SAMTOOLS_CMD = "samtools"


# Command utility functions

def run_cmd(args: List[Any], **kwargs):
    """ Run a command via subprocess.run(), with logging. """
    # Use shlex to place quotes around arguments containing whitespace.
    cmd = shlex.join(map(str, args))
    logger.debug(f"Shell $ {cmd}")
    return subprocess.run(cmd, check=True, shell=True, **kwargs)
