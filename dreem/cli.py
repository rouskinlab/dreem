from pipeline import run
from dreem.util.cli import *


@click.command()
@optgroup.group('I/O')
@opt_top_dir
@arg_fasta
@opt_fastqs
@opt_fastqi
@opt_fastq1
@opt_fastq2
@opt_library
@opt_phred_enc
@opt_min_phred
@opt_rerun
@optgroup.group('Selection')
@opt_coords
@opt_primers
@opt_fill
@opt_parallel
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
@opt_cpus
@optgroup.group('Aggregation')
@rnastructure_path
@rnastructure_temperature
@rnastructure_fold_args
@rnastructure_dms
@rnastructure_dms_min_unpaired_value
@rnastructure_dms_max_paired_value
@rnastructure_partition
@rnastructure_probability


@optgroup.group('Misc')
@verbose
def cli(*args, **kwargs):
    run(*args, **kwargs)


if __name__ == '__main__':
    cli()
