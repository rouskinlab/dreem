"""
Logging module of DREEM

Purpose
-------
Central manager of logging.
"""

import logging


MAX_VERBOSE = 2
MAX_QUIET = 2


def set_verbosity(verbose: int, quiet: int):
    """ Set the logging level based on the verbose and quiet arguments.

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
    if verbose == 0 and quiet == 2:
        level = logging.CRITICAL
    elif verbose == 0 and quiet == 1:
        level = logging.ERROR
    elif verbose == 0 and quiet == 0:
        level = logging.WARNING
    elif verbose == 1 and quiet == 0:
        level = logging.INFO
    elif verbose == 2 and quiet == 0:
        level = logging.DEBUG
    else:
        logging.warning(f"Invalid options: verbose={verbose}, quiet={quiet}. "
                        "Defaulting to verbose=0, quiet=0.")
        level = logging.WARNING
    logging.basicConfig(level=level)
    return level
