from click import Context, group, pass_context

from . import align, cluster, demultiplex, test, vector
from .util.cli import opt_quiet, opt_verbose
from .util.logio import set_verbosity


# Group for all DREEM commands
@group(params=[opt_verbose, opt_quiet],
       chain=True,
       context_settings={"show_default": True})
@pass_context
def cli(ctx: Context, verbose: int, quiet: int):
    """ DREEM command line interface """
    # Set verbosity level for logging.
    set_verbosity(verbose, quiet)
    # Ensure context object exists and is a dict.
    ctx.ensure_object(dict)


# Add all commands to the DREEM CLI command group.
cli.add_command(test.cli)
cli.add_command(demultiplex.cli)
cli.add_command(align.cli)
cli.add_command(vector.cli)
cli.add_command(cluster.cli)


def run(**kwargs):
    """ Run DREEM pipeline """
    pass

if __name__ == "__main__":
    cli()
