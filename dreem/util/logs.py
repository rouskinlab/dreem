"""
Logging configuration module of DREEM

Purpose
-------
Central manager of logging.
"""

import logging
import os
from subprocess import CompletedProcess
import sys

#import pkg_resources

#pkg_version = pkg_resources.get_distribution("dreem").version

#WELCOME = f"""
#Welcome to DREEM version {pkg_version}
#running on {sys.platform}
#with {os.cpu_count()} processors.
#"""



MAX_VERBOSE = 2
MAX_QUIET = 2
FILE_MSG_FORMAT = "LOGMSG>\t%(asctime)s\t%(name)s\t%(levelname)s\n%(message)s\n"
STREAM_MSG_FORMAT = "%(levelname)s>\t%(message)s"


LEVEL_NAME = {
    logging.DEBUG: "debug",
    logging.INFO: "info",
    logging.WARNING: "warning",
    logging.ERROR: "error",
    logging.CRITICAL: "critical",
}


def get_dreem_logger():
    """ Return the main DREEM logger. """
    if __name__ != (expect_name := "dreem.util.logs"):
        raise ValueError(
            f"{__file__} is named '{__name__}' (expected '{expect_name}')")
    dreem_logger_name = __name__.split(".")[0]
    return logging.getLogger(dreem_logger_name)


def get_verbosity(verbose: int = 0, quiet: int = 0):
    """ Get the logging level based on the verbose and quiet arguments.

    Parameters
    ----------
    verbose: int [0, 2]
        0 (): Log only warnings and errors
        1 (-v): Also log status updates
        2 (-vv): Also log detailed information (useful for debugging)
    quiet: int [0, 2]
        0 (): Suppress only status updates and detailed information
        1 (-q): Also suppress warnings
        2 (-qq): Also suppress non-critical error messages (discouraged)

    Giving both ```verbose``` and ```quiet``` flags causes the verbosity
    to default to ```verbose=0```, ```quiet=0```.
    """

    # Limit verbose and quiet to 2.
    if verbose > MAX_VERBOSE:
        logging.warning(f"Setting 'verbose' to {MAX_VERBOSE} (got {verbose})")
        verbose = MAX_VERBOSE
    if quiet > MAX_QUIET:
        logging.warning(f"Setting 'quiet' to {MAX_QUIET} (got {quiet})")
        quiet = MAX_QUIET

    # Set logging level based on verbose and quiet.
    if (verbose, quiet) == (2, 0):
        return logging.DEBUG
    if (verbose, quiet) == (1, 0):
        return logging.INFO
    if (verbose, quiet) == (0, 0):
        return logging.WARNING
    if (verbose, quiet) == (0, 1):
        return logging.ERROR
    if (verbose, quiet) == (0, 2):
        return logging.CRITICAL

    get_dreem_logger().warning(f"Invalid options: verbose={verbose}, "
                               f"quiet={quiet}. Setting both to 0")
    return get_verbosity(0, 0)


def config(verbose: int, quiet: int, log_file: str | None = None):
    """ Configure the main DREEM logger with handlers and verbosity. """
    # Set up logger.
    logger = get_dreem_logger()
    logger.setLevel(get_verbosity(verbose=MAX_VERBOSE))
    # Add stream handler.
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(get_verbosity(verbose, quiet))
    stream_handler.setFormatter(logging.Formatter(STREAM_MSG_FORMAT))
    logger.addHandler(stream_handler)
    # Add file handler.
    if log_file is not None:
        file_handler = logging.FileHandler(log_file, "a")
        file_handler.setLevel(get_verbosity(verbose=MAX_VERBOSE))
        file_handler.setFormatter(logging.Formatter(FILE_MSG_FORMAT))
        logger.addHandler(file_handler)
    #logger.info(WELCOME)


def log_process(logger: logging.Logger, process: CompletedProcess):
    """ Log the output and error messages of a process. """
    if process.stdout:
        logger.debug(f"STDOUT of {process.args}:\n{process.stdout.decode()}")
    if process.stderr:
        logger.debug(f"STDERR of {process.args}:\n{process.stderr.decode()}")
