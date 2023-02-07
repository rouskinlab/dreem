from dreem.vector.main import run
from dreem.util.cli import *


@click.command()
@opto_top_dir
@opti_library
@opti_coords
@opti_primers
@opti_fill
@opti_parallel
@opti_min_phred
@opti_phred_enc
@opto_rerun
@argi_fasta
@argi_bams
def cli(*args, **opts):
    run(*args, **opts)


if __name__ == '__main__':
    cli()
