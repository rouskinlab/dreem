import click

from dreem.util.cli import *
from dreem.stats.main import run


@click.command()
@opt_out_dir
@opt_stats_count
@opt_stats_frac
@opt_report
def cli(*args, **kwargs):
    run(*args, **kwargs)


if __name__ == "__main__":
    cli()
