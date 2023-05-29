from logging import getLogger
import os

from click import Context, group, pass_context

from . import profile
from ..core import docdef, path
from ..core.cli import merge_params

logger = getLogger(__name__)


all_params = merge_params(profile.params)


# Group for all graph commands
@group(path.MOD_GRAPH, params=[],
       invoke_without_command=True,
       context_settings={"show_default": True})
@pass_context
def cli(ctx: Context, **kwargs):
    """ Graphing command line interface """
    # If no subcommand was given, then run the default graphs.
    if ctx.invoked_subcommand is None:
        run(**kwargs)


@docdef.auto()
def run():
    return
