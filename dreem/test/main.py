from click import command, pass_obj

from ..util.cli import DreemCommandName, dreem_command
from .vector_tests import test_vectorize_read


@command(DreemCommandName.TEST.value)
# Pass context object
@pass_obj
# Turn into DREEM command
@dreem_command()
def cli():
    return run()


def run():
    test_vectorize_read.test_run()
