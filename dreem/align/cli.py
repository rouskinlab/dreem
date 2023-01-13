import click
from main import run
from click_option_group import optgroup
from dreem.util.cli import *

@click.command()

@argi_fasta
@opti_fastqu
@opti_fastqi
@opti_fastq1
@opti_fastq2
@opto_out_dir

def cli(*args, **kwargs):
    run(*args, **kwargs)

if __name__ == '__main__':
    cli()
