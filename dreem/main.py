import cProfile
import os

from click import Context, group, pass_context

from . import align, cluster, demultiplex, test, vector, aggregate
from .util import docdef
from .util.cli import (merge_params, opt_demultiplex, opt_cluster, opt_quiet,
                       opt_verbose, opt_profile)
from .util.logio import set_verbosity

logging_params = [
    opt_verbose,
    opt_quiet,
    opt_profile,
]

all_params = merge_params(logging_params,
                          [opt_demultiplex],
                          # demultiplex.params,
                          align.params,
                          vector.params,
                          [opt_cluster],
                          # cluster.params,
                          # aggregate.params,
                          aggregate.params,
                          # draw.params
                          )


# Group for all DREEM commands
@group(params=all_params,
       chain=True,
       invoke_without_command=True,
       context_settings={"show_default": True})
@pass_context
def cli(ctx: Context, verbose: int, quiet: int, profile: str, **kwargs):
    """ DREEM command line interface """
    # Set verbosity level for logging.
    set_verbosity(verbose, quiet)
    # Ensure context object exists and is a dict.
    ctx.ensure_object(dict)
    # If no subcommand was given, then run the entire pipeline.
    if ctx.invoked_subcommand is None:
        if profile:
            profile_path = os.path.abspath(profile)
            # Profile the program as it runs and write results to the
            # file given in the parameter profile.
            os.makedirs(os.path.dirname(profile_path), exist_ok=True)
            cProfile.runctx("run(**kwargs)",
                            globals=globals(),
                            locals=locals(),
                            filename=profile_path,
                            sort="time")
        else:
            # Run without profiling.
            run(**kwargs)


# Add all commands to the DREEM CLI command group.
cli.add_command(test.cli)
#cli.add_command(demultiplex.cli)
cli.add_command(align.cli)
cli.add_command(vector.cli)
cli.add_command(cluster.cli)
cli.add_command(aggregate.cli)


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
        parallel_demultiplexing:bool,
        clipped:int,
        mismatch_tolerence:int,
        index_tolerance:int,
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
        # Aggregation
        samples: str,
        rnastructure_path: str,
        bv_files: tuple[str],
        rnastructure_use_temp: bool,
        rnastructure_fold_args: str,
        rnastructure_use_dms: str,
        rnastructure_dms_min_unpaired_value: float,
        rnastructure_dms_max_paired_value: float,
        rnastructure_deltag_ensemble: bool,
        rnastructure_probability: bool,
        ):
    """ Run entire DREEM pipeline. """

    # Demultiplexing
    if demult_on:
        fastqs_dir_dm, fastqi_dir_dm, fastq12_dir_dm = demultiplex.run(
            library_csv=library,
            demulti_workspace=temp_dir,
            mixed_fastq1=fastq1,
            mixed_fastq2=fastq2,
            clipped=clipped,
            index_tolerance=index_tolerance,
            mismatch_tolerence=mismatch_tolerence,
            parallel=parallel_demultiplexing

        )
        fastqs = ()
        fastqi = ()
        fastq1 = ()
        fastq2 = ()
        
        fastqs_dir = fastqs_dir + fastqs_dir_dm
        fastqi_dir = fastqi_dir + fastqi_dir_dm
        fastq12_dir = fastq12_dir + fastq12_dir_dm
    # Alignment
    bamf += align.run(
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
        rem_buffer=rem_buffer)
    # Vectoring
    bv_files += vector.run(
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
    else:
        cluster_results = ""
    # Aggregate
    aggregate.main.run(
        out_dir=out_dir,
        temp_dir=temp_dir,
        save_temp=save_temp,
        library=library,
        samples=samples,
        fasta=fasta,
        bv_files=bv_files,
        clustering_file=cluster_results,
        rnastructure_path=rnastructure_path,
        rnastructure_use_temp=rnastructure_use_temp,
        rnastructure_fold_args=rnastructure_fold_args,
        rnastructure_use_dms=rnastructure_use_dms,
        rnastructure_dms_min_unpaired_value=rnastructure_dms_min_unpaired_value,
        rnastructure_dms_max_paired_value=rnastructure_dms_max_paired_value,
        rnastructure_deltag_ensemble=rnastructure_deltag_ensemble,
        rnastructure_probability=rnastructure_probability,
    )


if __name__ == "__main__":
    cli()
