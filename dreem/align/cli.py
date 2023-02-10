import click
from dreem.align.main import run
from dreem.util.cli import *

@click.command()
@arg_fasta
@opt_top_dir
@opt_fastqs
@opt_fastqi
@opt_fastq1
@opt_fastq2
@opt_fastqs_dir
@opt_fastqi_dir
@opt_fastq12_dir
@opt_phred_enc
@opt_parallel
@opt_cpus
def cli(*args, **kwargs):
    run(*args, **kwargs)


if __name__ == '__main__':
    cli()
