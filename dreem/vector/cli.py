from ..vector.main import run
from ..util.cli import *


@click.command()
# Input files
@opt_fasta
@opt_bamf
# Output directories
@opt_out_dir
@opt_temp_dir
# File generation
@opt_rerun
@opt_resume
@opt_save_temp
# Read quality
@opt_phred_enc
@opt_min_phred
# Regions
@opt_library
@opt_cfill
@opt_coords
@opt_primers
@opt_primer_gap
# Logging
@opt_verbose
@opt_quiet
@opt_logfile
# Parallelization
@opt_parallel
@opt_max_procs
def cli(*args, **opts):
    run(*args, **opts)


if __name__ == '__main__':
    cli()
