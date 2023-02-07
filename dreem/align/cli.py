import click
from dreem.align.main import run
from dreem.util.cli import argi_fasta, opti_fastqs, opti_fastqi, opti_fastq1, opti_fastq2, opti_fastqs_dir, opti_fastqi_dir, opti_fastq12_dir, opto_top_dir, opti_phred_enc


@click.command()
@argi_fasta
@opto_top_dir
@opti_fastqs
@opti_fastqi
@opti_fastq1
@opti_fastq2
@opti_fastqs_dir
@opti_fastqi_dir
@opti_fastq12_dir
@opti_phred_enc
def cli(*args, **kwargs):
    run(*args, **kwargs)


if __name__ == '__main__':
    cli()
