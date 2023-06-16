from click import group

#from . import canon as canon_mod
from . import pos as pos_mod, read as read_mod
from ..core import path


# Group for all graph commands
@group(path.MOD_GRAPH)
def cli():
    """ Graphing command line interface """


# Add graph commands to the CLI.
cli.add_command(pos_mod.cli)
cli.add_command(read_mod.cli)
