
from click import command, pass_obj


@click.command()

@optgroup.group('I/O')
@opt_library
@fastq1
@fastq2
@opt_out_dir

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
