import click
from dreem.align.main import run
from dreem.util.cli import *


@click.command()
@opt_out_dir
@opt_temp_dir
@arg_fasta
@opt_fastqs
@opt_fastqi
@opt_fastq1
@opt_fastq2
@opt_phred_enc
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
@opt_align_iseed
@opt_align_lseed
@opt_align_gbar
@opt_align_exten
@opt_align_reseed
@opt_align_pad
@opt_align_orient
@opt_parallel
@opt_max_procs
def cli(*args, **kwargs):
    run(*args, **kwargs)


if __name__ == '__main__':
    cli()
