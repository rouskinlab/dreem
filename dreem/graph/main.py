from logging import getLogger
import os

from click import Context, group, pass_context

from . import seq as seq_mod
from ..core import docdef, path
from ..core.cli import merge_params

logger = getLogger(__name__)


params = merge_params(seq_mod.params)


# Group for all graph commands
@group(path.MOD_GRAPH, params=params,
       invoke_without_command=True,
       context_settings={"show_default": True})
@pass_context
def cli(ctx: Context, **kwargs):
    """ Graphing command line interface """
    # If no subcommand was given, then run the default graphs.
    if ctx.invoked_subcommand is None:
        run(**kwargs)


cli.add_command(seq_mod.cli)


@docdef.auto()
def run():
    return
