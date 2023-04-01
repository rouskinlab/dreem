import logging
import shlex
import subprocess
from typing import Any, List

from ..util.logs import log_process


logger = logging.getLogger(__name__)


# Commands for external applications
BOWTIE2_CMD = "bowtie2"
BOWTIE2_BUILD_CMD = "bowtie2-build"
CUTADAPT_CMD = "cutadapt"
FASTQC_CMD = "fastqc"
SAMTOOLS_CMD = "samtools"


# Command utility functions

def run_cmd(args: List[Any]):
    """ Run a command via subprocess.run(), with logging. """
    # Use shlex to place quotes around arguments containing whitespace.
    cmd = shlex.join(map(str, args))
    # Log the command with which the process was run.
    logger.debug(f"Shell $ {cmd}")
    # Run the process and capture the output.
    process = subprocess.run(cmd, check=True, shell=True, capture_output=True)
    # Log the output of the process.
    log_process(logger, process)
    return process
