from pipeline import run
from dreem.util.cli import *


@click.command()
@opt_out_dir
@opt_temp_dir
@arg_fasta
@opt_fastqs
@opt_fastqi
@opt_fastq1
@opt_fastq2
@opt_library
@opt_phred_enc
@opt_min_phred
@opt_rerun
@opt_resume
@opt_trim
@opt_trim_adapt13
@opt_trim_adapt15
@opt_trim_adapt23
@opt_trim_adapt25
@opt_trim_minover
@opt_trim_maxerr
@opt_trim_minq1
@opt_trim_minq2
@opt_trim_minlen
@opt_trim_indels
@opt_trim_xtrim
@opt_trim_xuntrim
@opt_trim_nextseq
@opt_align_local
@opt_align_unal
@opt_align_disc
@opt_align_mixed
@opt_align_dove
@opt_align_cont
@opt_align_minl
@opt_align_maxl
@opt_align_score
@opt_align_sint
@opt_align_slen
@opt_align_gbar
@opt_align_exten
@opt_align_reseed
@opt_align_pad
@opt_align_orient
@opt_coords
@opt_primers
@opt_spanall
@opt_parallel
@opt_samples
@opt_demultiplex
@opt_max_barcode_mismatches
@opt_cluster
@opt_max_clusters
@opt_min_iter
@opt_signal_thresh
@opt_info_thresh
@opt_include_gu
@opt_include_del
@opt_min_reads
@opt_convergence_cutoff
@opt_num_runs
@opt_max_procs
@rnastructure_path
@rnastructure_temperature
@rnastructure_fold_args
@rnastructure_dms
@rnastructure_dms_min_unpaired_value
@rnastructure_dms_max_paired_value
@rnastructure_partition
@rnastructure_probability
@opt_verbose
@opt_quiet
@opt_logfile
def cli(*args, **kwargs):
    run(*args, **kwargs)


if __name__ == '__main__':
    cli()
