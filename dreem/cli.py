import click
from pipeline import run
from click_option_group import optgroup
from dreem.util.cli import *

@click.command()

@optgroup.group('I/O')
@opto_out_dir
@fasta
@fastq1
@fastq2
@input_dir
@samples
@sample
@clustering_file
@opti_library    

@optgroup.group('Selection')
@coords
@primers
@fill
@parallel

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

def cli(**args):
    run(**args)

if __name__ == '__main__':
    run()

