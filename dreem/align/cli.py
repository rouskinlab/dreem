import click
from dreem.align.main import run
from dreem.util.cli import argi_fasta, opti_fastqs, opti_fastqi, opti_fastq1, opti_fastq2, opti_fastqs_dir, opti_fastqi_dir, opti_fastq12_dir, opto_top_dir


@click.command()
@argi_fasta
@opti_fastqs
@opti_fastqi
@opti_fastq1
@opti_fastq2
@opti_fastqs_dir
@opti_fastqi_dir
@opti_fastq12_dir
@opto_top_dir
def cli(*args, **kwargs):
    run(*args, **kwargs)


if __name__ == '__main__':
    cli()
