from itertools import filterfalse
import logging
from pathlib import Path
import shlex
import subprocess
from typing import Any


logger = logging.getLogger(__name__)


# Commands for external applications
BOWTIE2_CMD = "bowtie2"
BOWTIE2_BUILD_CMD = "bowtie2-build"
CUTADAPT_CMD = "cutadapt"
FASTQC_CMD = "fastqc"
SAMTOOLS_CMD = "samtools"


# Command utility functions

def run_cmd(args: list[Any], verify_outputs: list[Path] | None = None):
    """ Run a command via subprocess.run(), with logging. """
    if verify_outputs is None:
        verify_outputs = list()
    verify_str = ", ".join(map(str, verify_outputs))
    # Use shlex to place quotes around arguments containing whitespace.
    cmd = shlex.join(map(str, args))
    # Check if any expected output file already exists.
    if exists := ", ".join(map(str, filter(Path.exists, verify_outputs))):
        raise FileExistsError("Expected output files already exist for "
                              f"{cmd}: {exists}")
    if verify_outputs:
        logger.debug(f"Verified outputs of {cmd} do not exist: {verify_str}")
    # Log the command with which the process was run.
    logger.debug(f"Shell $ {cmd}")
    # Run the process and capture the output.
    process = subprocess.run(cmd, check=True, shell=True, capture_output=True)
    # Log the output of the process.
    log_process(process)
    # Check if any expected output file is missing.
    if missing := ", ".join(map(str, filterfalse(Path.exists, verify_outputs))):
        raise FileNotFoundError("Failed to create expected output files of "
                                f"{cmd}: {missing}")
    if verify_outputs:
        logger.debug(f"Verified outputs of {cmd} exist: {verify_str}")
    return process


def log_process(process: subprocess.CompletedProcess):
    """ Log the output and error messages of a process. """
    if process.stdout:
        logger.info(f"STDOUT of {process.args}:\n{process.stdout.decode()}")
    if process.stderr:
        logger.info(f"STDERR of {process.args}:\n{process.stderr.decode()}")
