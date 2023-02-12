"""
Logging module of DREEM

Purpose
-------
Central manager of logging.
"""

import logging


def set_verbosity(verbose: int, quiet: int):
    """ Set the logging level based on the verbose and quiet arguments.

    Parameters
    ----------
    verbose: int (≥ 0, ≤ 2)
        Verbosity level
    quiet: int (≥ 0, ≤ 2)
        Quiet level
    """
    if 0 <= verbose <= 2 and 0 <= quiet <= 2 and (verbose == 0 or quiet == 0):
        if quiet == 2:  # verbose == 0
            level = logging.CRITICAL
        elif quiet == 1:  # verbose == 0
            level = logging.ERROR
        elif verbose == 1:  # quiet == 0
            level = logging.INFO
        elif verbose == 2:  # quiet == 0
            level = logging.DEBUG
        else:  # verbose == quiet == 0
            level = logging.WARNING
    else:
        logging.warning(f"Invalid options: verbose={verbose}, quiet={quiet}")
        level = logging.WARNING
    logging.basicConfig(level=level)
    return level
