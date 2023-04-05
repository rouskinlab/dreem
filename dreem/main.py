import cProfile
import os
import json

from click import Context, group, pass_context

from . import align, cluster, demultiplex, vector, aggregate, draw
from .util.dependencies import *
from .util import docdef, logs
from .util.cli import (merge_params, opt_demultiplex, opt_cluster,
                       opt_verbose, opt_quiet, opt_log, opt_profile,
                       opt_version)


logging_params = [
    opt_verbose,
    opt_quiet,
    opt_log,
    opt_profile,
]

misc_params = [
    opt_version,
]

all_params = merge_params(logging_params,
                          [opt_demultiplex],
                          demultiplex.params,
                          align.params,
                          vector.params,
                          [opt_cluster],
                          cluster.params,
                          aggregate.params,
                          aggregate.params,
                          draw.params,
                          misc_params,
                          )


# Group for all DREEM commands
@group(params=all_params,
       invoke_without_command=True,
       context_settings={"show_default": True})
@pass_context
def cli(ctx: Context, verbose: int, quiet: int, log: str, profile: str,
        **kwargs):
    """ DREEM command line interface """
    # Configure logging.
    logs.config(verbose, quiet, log_file=log)
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
# cli.add_command(test.cli)
cli.add_command(demultiplex.cli)
cli.add_command(align.cli)
cli.add_command(vector.cli)
cli.add_command(cluster.cli)
cli.add_command(aggregate.cli)
cli.add_command(draw.cli)


@docdef.auto()
def run(*,
        # General options
        out_dir: str,
        temp_dir: str,
        save_temp: bool,
        rerun: bool,
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
        demulti_overwrite: bool,
        demult_on: bool,
        parallel_demultiplexing: bool,
        clipped: int,
        mismatch_tolerence: int,
        index_tolerance: int,
        barcode_start: int,
        barcode_length: int,
        # Alignment options
        fastqs_dir: tuple[str],
        fastqi_dir: tuple[str],
        fastq12_dir: tuple[str],
        fastqc: bool,
        qc_extract: bool,
        cut: bool,
        cut_q1: int,
        cut_q2: int,
        cut_g1: tuple[str],
        cut_a1: tuple[str],
        cut_g2: tuple[str],
        cut_a2: tuple[str],
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
        # Vectoring
        bamf: tuple[str],
        bamd: tuple[str],
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, str, str], ...],
        primer_gap: int,
        autosect: bool,
        min_phred: int,
        ambid: bool,
        strict_pairs: bool,
        batch_size: float,
        # Clustering
        clust: bool,
        mp_report: tuple[str],
        max_clusters: int,
        num_runs: int,
        signal_thresh: float,
        include_gu: bool,
        include_del: bool,
        polya_max: int,
        min_iter: int,
        max_iter: int,
        convergence_cutoff: float,
        min_reads: int,
        # Aggregation
        samples: str,
        rnastructure_path: str,
        bv_files: tuple[str],
        clustering_file: str,
        rnastructure_use_temp: bool,
        rnastructure_fold_args: str,
        rnastructure_use_dms: str,
        rnastructure_dms_min_unpaired_value: float,
        rnastructure_dms_max_paired_value: float,
        rnastructure_deltag_ensemble: bool,
        rnastructure_probability: bool,
        # Drawing
        inpt: tuple[str],
        flat: bool,
        section: str,
        mutation_fraction: bool,
        mutation_fraction_identity: bool,
        base_coverage: bool,
        mutation_per_read_per_reference: bool,
        mutations_in_barcodes: bool,
        # Misc
        version: bool
        ):

    if version:
        #print(f"DREEM version {__version__}")
        return 0

    check_bowtie2_exists()
    check_cutadapt_exists()
    check_fastqc_exists()
    check_samtools_exists()
    check_rnastructure_exists(rnastructure_path)
    
    """ Run entire DREEM pipeline. """
    # Demultiplexing
    if demult_on:
        fastqs_dir_dm, fastqi_dir_dm, fastq12_dir_dm = demultiplex.run(
            fasta=fasta,
            library=library,
            out_dir=out_dir,
            temp_dir=temp_dir,
            demulti_overwrite=demulti_overwrite,
            fastq1=fastq1,
            fastq2=fastq2,
            clipped=clipped,
            index_tolerance=index_tolerance,
            mismatch_tolerence=mismatch_tolerence,
            parallel_demultiplexing=parallel_demultiplexing,
            barcode_start=barcode_start,
            barcode_length=barcode_length,

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
        qc_extract=qc_extract,
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
        bt2_orient=bt2_orient)
    # Vectoring
    bv_files += vector.run(
        out_dir=out_dir,
        temp_dir=temp_dir,
        save_temp=save_temp,
        rerun=rerun,
        max_procs=max_procs,
        parallel=parallel,
        fasta=fasta,
        bamf=bamf,
        bamd=bamd,
        library=library,
        autosect=autosect,
        coords=coords,
        primers=primers,
        primer_gap=primer_gap,
        phred_enc=phred_enc,
        min_phred=min_phred,
        ambid=ambid,
        strict_pairs=strict_pairs,
        batch_size=batch_size,
    )
    if clust:
        cluster_results = cluster.run(
            mp_report=bv_files,
            max_clusters=max_clusters,
            min_iter=min_iter,
            max_iter=max_iter,
            signal_thresh=signal_thresh,
            include_gu=include_gu,
            include_del=include_del,
            polya_max=polya_max,
            min_reads=min_reads,
            convergence_cutoff=convergence_cutoff,
            num_runs=num_runs,
            max_procs=max_procs
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

    draw.run(
        inpt=list(inpt) + [json.load(open(os.path.join(out_dir, f), 'r')) for f in os.listdir(out_dir) if
                           f.endswith(".json")],
        out_dir=out_dir,
        flat=flat,
        mutation_fraction=mutation_fraction,
        mutation_fraction_identity=mutation_fraction_identity,
        base_coverage=base_coverage,
        mutation_per_read_per_reference=mutation_per_read_per_reference,
        mutations_in_barcodes=mutations_in_barcodes,
        section=section
    )


if __name__ == "__main__":
    cli()
