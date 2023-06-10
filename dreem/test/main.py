from os.path import dirname
import unittest as ut

from click import command

from ..core import docdef, path
from ..core.cli import opt_temp_dir, opt_save_temp, opt_verbose
from ..core.parallel import lock_temp_dir


# Parameters for command line interface
params = [opt_temp_dir, opt_save_temp, opt_verbose]


@command(path.MOD_TEST, params=params)
def cli(**kwargs):
    """ Run all unit tests to ensure DREEM works properly. """
    return run(**kwargs)


@lock_temp_dir
@docdef.auto()
def run(temp_dir: str, save_temp: bool, verbose: int):
    """
    Run all unit tests for DREEM.
    """
    # Discover all unit test modules in DREEM.
    dreem_main_dir = dirname(dirname(__file__))
    # The line top_level_dir=dirname(dreem_main_dir) is required to make
    # Python consider dreem as a package, so that relative imports work.
    # Omitting this line raises an ImportError during every test.
    suite = ut.TestLoader().discover(dreem_main_dir,
                                     top_level_dir=dirname(dreem_main_dir))
    # Run all unit tests.
    runner = ut.TextTestRunner(verbosity=verbose)
    runner.run(suite)
