import click
from main import run
from click_option_group import optgroup
from dreem.util.cli import *

@click.command()

@optgroup.group('I/O')
@fasta
@fastq1
@fastq2
@opto_out_dir

@optgroup.group('Alignment')
@demultiplexing

@optgroup.group('Selection')
@coords
@primers
@fill
@parallel

@optgroup.group('Miscellaneous')
@verbose

def cli(**args):
    run(**args)

if __name__ == '__main__':
    run()
