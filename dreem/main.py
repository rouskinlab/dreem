import cProfile

from click import Context, group, pass_context

from . import (demultiplex as demultiplex_mod,
               align as align_mod,
               relate as relate_mod,
               mask as call_mod,
               cluster as cluster_mod,
               table as table_mod,
               fold as fold_mod,
               aggregate as agg_mod,
               draw as draw_mod)
from .core.dependencies import *
from .core import docdef, logs
from .core.cli import (merge_params, opt_demultiplex,
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
                          demultiplex_mod.params,
                          align_mod.params,
                          relate_mod.params,
                          call_mod.params,
                          cluster_mod.params,
                          table_mod.params,
                          fold_mod.params,
                          agg_mod.params,
                          draw_mod.params,
                          misc_params)


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
cli.add_command(demultiplex_mod.cli)
cli.add_command(align_mod.cli)
cli.add_command(relate_mod.cli)
cli.add_command(call_mod.cli)
cli.add_command(cluster_mod.cli)
cli.add_command(table_mod.cli)
cli.add_command(fold_mod.cli)
cli.add_command(agg_mod.cli)
cli.add_command(draw_mod.cli)


@docdef.auto()
def run(*,
        # General options
        out_dir: str,
        temp_dir: str,
        save_temp: bool,
        rerun: bool,
        max_procs: int,
        parallel: bool,
        # Input files
        fasta: str,
        fastqs: tuple[str, ...],
        fastqi: tuple[str, ...],
        fastqm: tuple[str, ...],
        phred_enc: int,
        bam: tuple[str, ...],
        rel: tuple[str, ...],
        call: tuple[str, ...],
        clust: tuple[str, ...],
        table: tuple[str, ...],
        # Demultiplexing
        demulti_overwrite: bool,
        demult_on: bool,
        parallel_demultiplexing: bool,
        clipped: int,
        mismatch_tolerence: int,
        index_tolerance: int,
        barcode_start: int,
        barcode_length: int,
        # Alignment
        dmfastqs: tuple[str, ...],
        dmfastqi: tuple[str, ...],
        dmfastqm: tuple[str, ...],
        fastqc: bool,
        qc_extract: bool,
        cut: bool,
        cut_q1: int,
        cut_q2: int,
        cut_g1: tuple[str, ...],
        cut_a1: tuple[str, ...],
        cut_g2: tuple[str, ...],
        cut_a2: tuple[str, ...],
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
        bt2_score_min_e2e: str,
        bt2_score_min_loc: str,
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
        min_phred: int,
        ambid: bool,
        batch_size: float,
        # Filtering
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, str, str], ...],
        primer_gap: int,
        library: str,
        count_del: bool,
        count_ins: bool,
        discount_mut: tuple[str, ...],
        exclude_polya: int,
        exclude_gu: bool,
        exclude_pos: tuple[tuple[str, int], ...],
        min_finfo_read: float,
        max_fmut_read: int,
        min_mut_gap: int,
        min_ninfo_pos: int,
        max_fmut_pos: float,
        # Clustering
        max_clusters: int,
        em_runs: int,
        min_em_iter: int,
        max_em_iter: int,
        em_thresh: float,
        # Folding
        dms_quantile: float,
        # Aggregation
        samples: str,
        rnastructure_path: str,
        rnastructure_use_temp: bool,
        rnastructure_fold_args: str,
        rnastructure_use_dms: str,
        rnastructure_dms_min_unpaired_value: float,
        rnastructure_dms_max_paired_value: float,
        rnastructure_deltag_ensemble: bool,
        rnastructure_probability: bool,
        # Drawing
        inpt: tuple[str, ...],
        flat: bool,
        section: str,
        mutation_fraction: bool,
        mutation_fraction_identity: bool,
        base_coverage: bool,
        mutation_per_read_per_reference: bool,
        mutations_in_barcodes: bool,
        # Misc
        version: bool):
    if version:
        # print(f"DREEM version {__version__}")
        return 0

    check_bowtie2_exists()
    check_cutadapt_exists()
    check_fastqc_exists()
    check_samtools_exists()
    check_rnastructure_exists(rnastructure_path)

    """ Run entire DREEM pipeline. """
    # Demultiplexing
    if demult_on:
        for dms, dmi, dmm in demultiplex_mod.run(
                fasta=fasta,
                library=library,
                out_dir=out_dir,
                temp_dir=temp_dir,
                demulti_overwrite=demulti_overwrite,
                fastqm=fastqm,
                clipped=clipped,
                index_tolerance=index_tolerance,
                mismatch_tolerence=mismatch_tolerence,
                parallel_demultiplexing=parallel_demultiplexing,
                barcode_start=barcode_start,
                barcode_length=barcode_length,
                phred_enc=phred_enc):
            dmfastqs = dmfastqs + dms
            dmfastqi = dmfastqi + dmi
            dmfastqm = dmfastqm + dmm
    # Alignment
    bam += tuple(map(str, align_mod.run(
        out_dir=out_dir,
        temp_dir=temp_dir,
        save_temp=save_temp,
        rerun=rerun,
        max_procs=max_procs,
        parallel=parallel,
        fasta=fasta,
        fastqs=fastqs,
        fastqi=fastqi,
        fastqm=fastqm,
        dmfastqs=dmfastqs,
        dmfastqi=dmfastqi,
        dmfastqm=dmfastqm,
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
        bt2_score_min_e2e=bt2_score_min_e2e,
        bt2_score_min_loc=bt2_score_min_loc,
        bt2_i=bt2_i,
        bt2_x=bt2_x,
        bt2_gbar=bt2_gbar,
        bt2_l=bt2_l,
        bt2_s=bt2_s,
        bt2_d=bt2_d,
        bt2_r=bt2_r,
        bt2_dpad=bt2_dpad,
        bt2_orient=bt2_orient
    )))
    # Mutation Vectoring
    rel += tuple(map(str, relate_mod.run(
        fasta=fasta,
        bam=bam,
        out_dir=out_dir,
        temp_dir=temp_dir,
        phred_enc=phred_enc,
        min_phred=min_phred,
        ambid=ambid,
        batch_size=batch_size,
        max_procs=max_procs,
        parallel=parallel,
        rerun=rerun,
        save_temp=save_temp,
    )))
    # Filtering
    call += tuple(map(str, call_mod.run(
        rel=rel,
        coords=coords,
        primers=primers,
        primer_gap=primer_gap,
        library=library,
        count_del=count_del,
        count_ins=count_ins,
        discount_mut=discount_mut,
        exclude_polya=exclude_polya,
        exclude_gu=exclude_gu,
        exclude_pos=exclude_pos,
        min_finfo_read=min_finfo_read,
        max_fmut_read=max_fmut_read,
        min_mut_gap=min_mut_gap,
        min_ninfo_pos=min_ninfo_pos,
        max_fmut_pos=max_fmut_pos,
        max_procs=max_procs,
        parallel=parallel,
        rerun=rerun,
    )))
    # Clustering
    clust += tuple(map(str, cluster_mod.run(
        call=call,
        max_clusters=max_clusters,
        em_runs=em_runs,
        min_em_iter=min_em_iter,
        max_em_iter=max_em_iter,
        em_thresh=em_thresh,
        max_procs=max_procs,
        parallel=parallel,
        rerun=rerun,
    )))
    # Table
    table += tuple(map(str, table_mod.run(
        rel=rel,
        call=call,
        clust=clust,
        max_procs=max_procs,
        parallel=parallel,
        rerun=rerun,
    )))
    fold = tuple(map(str, fold_mod.run(
        table=table,
        fasta=fasta,
        library=library,
        coords=coords,
        primers=primers,
        primer_gap=primer_gap,
        dms_quantile=dms_quantile,
        temp_dir=temp_dir,
        save_temp=save_temp,
        max_procs=max_procs,
        parallel=parallel,
        rerun=rerun,
    )))
    '''
    # Aggregate
    aggregate_mod.run(
        mv_file=mvec,
        clust_file=clust,
        out_dir=out_dir,
        temp_dir=temp_dir,
        save_temp=save_temp,
        coords=coords,
        primers=primers,
        primer_gap=primer_gap,
        library=library,
        samples=samples,
        rnastructure_path=rnastructure_path,
        rnastructure_use_temp=rnastructure_use_temp,
        rnastructure_fold_args=rnastructure_fold_args,
        rnastructure_use_dms=rnastructure_use_dms,
        rnastructure_dms_min_unpaired_value=rnastructure_dms_min_unpaired_value,
        rnastructure_dms_max_paired_value=rnastructure_dms_max_paired_value,
        rnastructure_deltag_ensemble=rnastructure_deltag_ensemble,
        rnastructure_probability=rnastructure_probability,
    )
    # Drawing
    draw_mod.run(
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
    '''


if __name__ == "__main__":
    cli()
