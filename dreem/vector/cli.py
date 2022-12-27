import click
from main import run
from click_option_group import optgroup
from dreem.util.cli import *


@click.command()
@opti_library
@opti_coords
@opti_primers
@opti_fill
@opti_parallel
@opto_out_dir
@argi_fasta
@argi_bams
def cli(*args, **opts):
    run(*args, **opts)


if __name__ == '__main__':
    cli()
