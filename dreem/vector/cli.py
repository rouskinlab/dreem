import click
from main import run
from click_option_group import optgroup
from dreem.util.cli import *

@click.command()

@optgroup.group('I/O')
@fasta
@input_dir
@opto_out_dir
@opti_library

@optgroup.group('Selection')
@coords
@primers
@fill

@optgroup.group('Miscellaneous')
@verbose

def cli(**args):
    run(**args)
    
if __name__ == '__main__':
    run()