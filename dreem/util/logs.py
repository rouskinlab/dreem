"""
Logging configuration module of DREEM

Purpose
-------
Central manager of logging.
"""

import logging
from subprocess import CompletedProcess


MAX_VERBOSE = 2
MAX_QUIET = 2
MIN_FILE_VERBOSITY = 2
MIN_STREAM_VERBOSITY = -2
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


def get_verbosity(verbose: int, quiet: int, minimum: int | None = None):
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
    minimum: int | None
        Ensure that the verbosity level is at least minimum, even if the
        verbose and quiet arguments specify less than this minimum. Must
        be an integer in [-2, 2], or None to set to minimum.

    Giving both ```verbose``` and ```quiet``` flags causes the verbosity
    to default to ```verbose=0```, ```quiet=0```.
    """

    # Limit verbose and quiet to 2.
    if minimum is None:
        minimum = -MAX_QUIET
    if verbose > MAX_VERBOSE:
        logging.warning(f"Setting 'verbose' to {MAX_VERBOSE} (got {verbose})")
        verbose = MAX_VERBOSE
    if quiet > MAX_QUIET:
        logging.warning(f"Setting 'quiet' to {MAX_QUIET} (got {quiet})")
        quiet = MAX_QUIET

    # Set logging level based on verbose and quiet.
    if minimum == 2 or (verbose, quiet) == (2, 0):
        return logging.DEBUG
    if minimum == 1 or (verbose, quiet) == (1, 0):
        return logging.INFO
    if minimum == 0 or (verbose, quiet) == (0, 0):
        return logging.WARNING
    if minimum == -1 or (verbose, quiet) == (0, 1):
        return logging.ERROR
    if minimum == -2 or (verbose, quiet) == (0, 2):
        return logging.CRITICAL

    get_dreem_logger().warning(f"Invalid options: verbose={verbose}, "
                               f"quiet={quiet}, minimum={minimum}.")
    return get_verbosity(0, 0)


def config(verbose: int, quiet: int, log_file: str | None = None):
    """ Configure the main DREEM logger with handlers and verbosity. """
    # Set up logger.
    logger = get_dreem_logger()
    logger.setLevel(get_verbosity(verbose, quiet, MAX_VERBOSE))
    # Add stream handler.
    shdlr = logging.StreamHandler()
    shdlr.setLevel(get_verbosity(verbose, quiet, MIN_STREAM_VERBOSITY))
    shdlr.setFormatter(logging.Formatter(STREAM_MSG_FORMAT))
    logger.addHandler(shdlr)
    # Add file handler.
    if log_file is not None:
        fhdlr = logging.FileHandler(log_file, "a")
        fhdlr.setLevel(get_verbosity(verbose, quiet, MIN_FILE_VERBOSITY))
        fhdlr.setFormatter(logging.Formatter(FILE_MSG_FORMAT))
        logger.addHandler(fhdlr)


def log_process(logger: logging.Logger, process: CompletedProcess):
    """ Log the output and error messages of a process. """
    if process.stdout:
        logger.debug(f"STDOUT:\n{process.stdout.decode()}")
    if process.stderr:
        logger.debug(f"STDERR:\n{process.stderr.decode()}")
