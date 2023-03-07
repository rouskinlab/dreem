import logging
import os
from click import Context, group, pass_context

from . import align, cluster, demultiplex, test, vector, aggregate
from .util import docdef
from .util.cli import (merge_params, opt_demultiplex, opt_cluster, opt_quiet,
                       opt_verbose)
from .util.logio import set_verbosity
from .util import path

verbose_params = [
    opt_verbose,
    opt_quiet,
]

all_params = merge_params(verbose_params,
                          [opt_demultiplex],
                          # demultiplex.params,
                          align.params,
                          vector.params,
                          [opt_cluster],
                          # cluster.params,
                          # aggregate.params,
                          # draw.params
                          )


# Group for all DREEM commands
@group(params=all_params, 
       chain=True,
       invoke_without_command=True,
       context_settings={"show_default": True})
@pass_context
def cli(ctx: Context, verbose: int, quiet: int, **kwargs):
    """ DREEM command line interface """
    # Set verbosity level for logging.
    set_verbosity(verbose, quiet)
    # Ensure context object exists and is a dict.
    ctx.ensure_object(dict)
    # If no subcommand was given, then run the entire pipeline.
    if ctx.invoked_subcommand is None:
        run(**kwargs)
 

# Add all commands to the DREEM CLI command group.
#cli.add_command(test.cli)
cli.add_command(demultiplex.cli)
cli.add_command(align.cli)
cli.add_command(vector.cli)
cli.add_command(cluster.cli)


@docdef.auto()
def run(*,
        # General options
        out_dir: str,
        temp_dir: str,
        save_temp: bool,
        rerun: bool,
        resume: bool,
        max_procs: int,
        parallel: bool,
        library: str,
        fasta: str,
        fastqs: tuple[str],
        fastqi: tuple[str],
        fastq1: tuple[str],
        fastq2: tuple[str],
        phred_enc: int,
        # Demultiplexing options
        demult_on: bool,
        # FIXME: add parameters for demultiplexing
        # Alignment options
        fastqs_dir: tuple[str],
        fastqi_dir: tuple[str],
        fastq12_dir: tuple[str],
        fastqc: bool,
        fastqc_extract: bool,
        cut: bool,
        cut_q1: int,
        cut_q2: int,
        cut_g1: str,
        cut_a1: str,
        cut_g2: str,
        cut_a2: str,
        cut_o: int,
        cut_e: float,
        cut_indels: bool,
        cut_nextseq: bool,
        cut_discard_trimmed: bool,
        cut_discard_untrimmed: bool,
        cut_m: int,
        bt2_local: bool,
        bt2_discordant: bool,
        bt2_mixed: bool,
        bt2_dovetail: bool,
        bt2_contain: bool,
        bt2_unal: bool,
        bt2_score_min: str,
        bt2_i: int,
        bt2_x: int,
        bt2_gbar: int,
        bt2_l: int,
        bt2_s: str,
        bt2_d: int,
        bt2_r: int,
        bt2_dpad: int,
        bt2_orient: str,
        rem_buffer: int,
        # Vectoring
        bamf: tuple[str],
        bamd: tuple[str],
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, str, str], ...],
        primer_gap: int,
        cfill: bool,
        min_phred: int,
        ambid: bool,
        strict_pairs: bool,
        batch_size: float,
        # Clustering
        cluster_on: bool,
        # Aggregate
        samples: str,
        rnastructure_path: str,
        ):
    """ Run entire DREEM pipeline. """
    # Demultiplexing
    if demult_on:
        fastqs_dir_dm, fastqi_dir_dm, fastq12_dir_dm = demultiplex.run(
            # FIXME: add arguments
        )
        fastqs = ()
        fastqi = ()
        fastq1 = ()
        fastq2 = ()
        fastqs_dir = fastqs_dir + fastqs_dir_dm
        fastqi_dir = fastqi_dir + fastqi_dir_dm
        fastq12_dir = fastq12_dir + fastq12_dir_dm
    # Alignment
    bams_aln = align.run(
        out_dir=out_dir,
        temp_dir=temp_dir,
        save_temp=save_temp,
        rerun=rerun,
        resume=resume,
        max_procs=max_procs,
        parallel=parallel,
        fasta=fasta,
        fastqs=fastqs,
        fastqi=fastqi,
        fastq1=fastq1,
        fastq2=fastq2,
        fastqs_dir=fastqs_dir,
        fastqi_dir=fastqi_dir,
        fastq12_dir=fastq12_dir,
        phred_enc=phred_enc,
        fastqc=fastqc,
        fastqc_extract=fastqc_extract,
        cut=cut,
        cut_q1=cut_q1,
        cut_q2=cut_q2,
        cut_g1=cut_g1,
        cut_a1=cut_a1,
        cut_g2=cut_g2,
        cut_a2=cut_a2,
        cut_o=cut_o,
        cut_e=cut_e,
        cut_indels=cut_indels,
        cut_nextseq=cut_nextseq,
        cut_discard_trimmed=cut_discard_trimmed,
        cut_discard_untrimmed=cut_discard_untrimmed,
        cut_m=cut_m,
        bt2_local=bt2_local,
        bt2_discordant=bt2_discordant,
        bt2_mixed=bt2_mixed,
        bt2_dovetail=bt2_dovetail,
        bt2_contain=bt2_contain,
        bt2_unal=bt2_unal,
        bt2_score_min=bt2_score_min,
        bt2_i=bt2_i,
        bt2_x=bt2_x,
        bt2_gbar=bt2_gbar,
        bt2_l=bt2_l,
        bt2_s=bt2_s,
        bt2_d=bt2_d,
        bt2_r=bt2_r,
        bt2_dpad=bt2_dpad,
        bt2_orient=bt2_orient,
        rem_buffer=rem_buffer,
    )
    bamf = bamf + bams_aln
    # Vectoring
    profiles_vec = vector.run(
        out_dir=out_dir,
        temp_dir=temp_dir,
        save_temp=save_temp,
        rerun=rerun,
        resume=resume,
        max_procs=max_procs,
        parallel=parallel,
        fasta=fasta,
        bamf=bamf,
        bamd=bamd,
        library=library,
        cfill=cfill,
        coords=coords,
        primers=primers,
        primer_gap=primer_gap,
        phred_enc=phred_enc,
        min_phred=min_phred,
        ambid=ambid,
        strict_pairs=strict_pairs,
        batch_size=batch_size,
    )
    if cluster_on:
        cluster_results = cluster.run(
            # FIXME: add arguments
        )
    # Aggregate
    # TODO: make sample better
    sample = [f.split("/")[-1].split("_R1.f")[0] for f in fastq1]
    for s in sample:
        aggregate.main.run(
            out_dir=out_dir,
            fasta=fasta,
            library=library,
            samples=samples,
            clustering_file = cluster_results if cluster_on else None,
            bv_files = profiles_vec,
            sample = s,
            rnastructure_path = rnastructure_path,
        )

if __name__ == "__main__":
    cli()
