from pipeline import run
from dreem.util.cli import *


@click.command()
@optgroup.group('I/O')
@opto_top_dir
@argi_fasta
@opti_fastqs
@opti_fastqi
@opti_fastq1
@opti_fastq2
@opti_library
@optgroup.group('Selection')
@opti_coords
@opti_primers
@opti_fill
@opti_parallel
@optgroup.group('Demultiplexing')
@demultiplexing
@barcode_start
@barcode_length
@max_barcode_mismatches
@optgroup.group('Clustering')
@clustering
@max_clusters
@min_iter
@signal_thresh
@info_thresh
@include_g_u
@include_del
@min_reads
@convergence_cutoff
@num_runs
@n_cpus
@optgroup.group('Aggregation')
@rnastructure_path
@rnastructure_temperature
@rnastructure_fold_args
@rnastructure_dms
@rnastructure_dms_min_unpaired_value
@rnastructure_dms_max_paired_value
@rnastructure_partition
@rnastructure_probability
@poisson
@optgroup.group('Misc')
@verbose
def cli(*args, **kwargs):
    run(*args, **kwargs)


if __name__ == '__main__':
    cli()
