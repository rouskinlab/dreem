from dreem.align.main import run
from dreem.util.cli import *


@click.command()
# Input files
@opt_fasta
@opt_fastqs
@opt_fastqi
@opt_fastq1
@opt_fastq2
# Output directories
@opt_out_dir
@opt_temp_dir
# File generation
@opt_rerun
@opt_resume
@opt_save_temp
# Quality control
@opt_phred_enc
@opt_fastqc
@opt_fastqc_extract
# Trimming
@opt_cutadapt
@opt_cut_a1
@opt_cut_g1
@opt_cut_a2
@opt_cut_g2
@opt_cut_o
@opt_cut_e
@opt_cut_q1
@opt_cut_q2
@opt_cut_m
@opt_cut_indels
@opt_cut_discard_trimmed
@opt_cut_discard_untrimmed
@opt_cut_nextseq
# Alignment
@opt_bt2_local
@opt_bt2_discordant
@opt_bt2_mixed
@opt_bt2_dovetail
@opt_bt2_contain
@opt_bt2_i
@opt_bt2_x
@opt_bt2_score_min
@opt_bt2_s
@opt_bt2_l
@opt_bt2_gbar
@opt_bt2_d
@opt_bt2_r
@opt_bt2_dpad
@opt_bt2_orient
# Logging
@opt_verbose
@opt_quiet
@opt_logfile
# Parallelization
@opt_parallel
@opt_max_procs
def cli(*args, **kwargs):
    run(*args, **kwargs)


if __name__ == '__main__':
    cli()
