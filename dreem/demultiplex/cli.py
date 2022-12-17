import click
from main import run
from click_option_group import optgroup
from dreem.util.cli_args import *

@click.command()

@optgroup.group('I/O')
@library
@fastq1
@fastq2
@out_dir

@optgroup.group('Demultiplexing')
@barcode_start
@barcode_length
@max_barcode_mismatches

@optgroup.group('Selection')
@coords
@primers
@fill
@interleaved

@optgroup.group('Miscellaneous')
@verbose

def cli(**args):
    run(**args)

if __name__ == '__main__':
    cli()