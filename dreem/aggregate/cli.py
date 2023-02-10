import click
from main import run
from click_option_group import optgroup
from dreem.util.cli import *

@click.command()

@optgroup.group('I/O')
@input_dir
@opt_library
@opt_top_dir
@samples
@sample
@clustering_file

@optgroup.group('Selection')
@coords
@primers
@fill

@optgroup.group('RNAstructure')
@rnastructure_path
@rnastructure_temperature
@rnastructure_fold_args
@rnastructure_dms
@rnastructure_dms_min_unpaired_value
@rnastructure_dms_max_paired_value
@rnastructure_partition
@rnastructure_probability

@optgroup.group('Miscellaneous')
@verbose


def cli(**args):
    run(**args)

if __name__ == '__main__':
    run()
