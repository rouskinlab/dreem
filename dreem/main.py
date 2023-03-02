from click import Context, group, pass_context

from .align import cli as align_cli
from .cluster import cli as cluster_cli
from .demultiplex import cli as demultiplex_cli
from .util.logio import set_verbosity
from .util.cli import opt_quiet, opt_verbose
from .vector import cli as vector_cli


# Group for all DREEM commands
@group(params=[opt_verbose, opt_quiet],
       chain=True,
       context_settings={"show_default": True})
@pass_context
def cli(ctx: Context, verbose: int, quiet: int):
    """ Main function of the DREEM command line interface (CLI) """
    # Set verbosity level for logging.
    set_verbosity(verbose, quiet)
    # Ensure context object exists and is a dict.
    ctx.ensure_object(dict)


cli.add_command(demultiplex_cli)
cli.add_command(align_cli)
cli.add_command(vector_cli)
cli.add_command(cluster_cli)

'''

@opt_fastqs
@opt_fastqi
@opt_fastq1
@opt_fastq2
@opt_library
@opt_phred_enc
@opt_min_phred
@opt_fastqc
@opt_fastqc_extract
@opt_rerun
@opt_resume
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
@opt_coords
@opt_primers
@opt_primer_gap
@opt_cfill
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
def run(out_dir: str,
        temp_dir: str,
        save_temp: bool,
        rerun: bool,
        resume: bool,
        demultiplex: bool,
        cluster: bool,
        fasta: str,
        fastqs: tuple[str],
        fastqi: tuple[str],
        fastq1: tuple[str],
        fastq2: tuple[str],
        library: str,
        samples: str,
        n_procs: int,
        phred_enc: int,
        min_phred: int,
        fastqc: bool,
        fastqc_extract: bool,
        trim: bool,
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
        bt2_unal: bool,
        bt2_discordant: bool,
        bt2_mixed: bool,
        bt2_dovetail: bool,
        bt2_contain: bool,
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
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, str, str], ...],
        primer_gap: int,
        cfill: bool,
        parallel: bool,
        max_barcode_mismatches: int,
        max_clusters: int,
        signal_thresh: float,
        info_thresh: float,
        include_gu: bool,
        include_del: bool,
        min_reads: int,
        convergence_cutoff: float,
        min_iter: int,
        num_runs: int,
        rnastructure_path: str,
        rnastructure_temperature: int,
        rnastructure_fold_args: str,
        rnastructure_dms: bool,
        rnastructure_dms_min_unpaired_value: int,
        rnastructure_dms_max_paired_value: int,
        rnastructure_partition: bool,
        rnastructure_probability: bool,
        verbose: int,
        quiet: int,
        log: str):
    """ Run the entire DREEM pipeline from command-line arguments.

    Parameters
    ----------
    out_dir: str (default: "./output")
        Path to the directory where all output files will be written
    temp_dir: str (default: "./temp")
        Path to the directory where all temporary files will be written
    rerun: bool (default: False)
        Whether to rerun the creation of output files that already exist
    resume: bool (default: False)
        Whether to resume processes that terminated before creating all
        output files by resuming at the most recent temporary file
    parallel: bool
        Whether to allow multiple tasks to run in parallel
    n_procs: int (≥ 1; default: os.cpu_count())
        Maximum number of CPUs. If 1, parallelization is turned off.
    demultiplex: bool (default: False)
        Whether to run demultiplexing
    cluster: bool (default: False)
        Whether to run clustering
    fasta: str
        Path to the FASTA file of reference sequences
    fastqs: tuple[str]
        Paths to FASTQ files of single-end reads
    fastqi: tuple[str]
        Paths to FASTQ files of interleaved paired-end reads
    fastq1: tuple[str]
        Paths to FASTQ files of mate 1 paired-end reads
    fastq2: tuple[str]
        Paths to FASTQ files of mate 2 paired-end reads
    library: str
        Path to the library file
    samples: str
        Path to the samples file
    coords: tuple[tuple[str, int, int], ...]
        Define each region by its reference name and 5' and 3' coord-
        inates (1-indexed, inclusive both first and last positions):
        ```--coords ref-1 end5-1 end3-1 --coords ref-2 end5-2 end3-2```
    primers: tuple[tuple[str, str, str], ...] (for amplicons only)
        Define each amplicon-based region by its reference name and
        sequences of the forward and reverse primers of the amplicon:
        ```--primers ref-1 fwd-1 rev-1 --primers ref-2 fwd-2 rev-2```
        Note: give actual reverse primers, not the reverse complements.
    cfill: bool (default: False)
        Whether to create a region that spans the entire sequence of
        each reference for which no coordinates or primers were given
    barcode_start: int
        Start position of the barcode for demultiplexing
    barcode_length: int
        Length of the barcode for demultiplexing
    max_barcode_mismatches: int
        Maximum number of mismatches in the barcode for demultiplexing
    max_clusters: int
        Maximum number of clusters
    signal_thresh: float
        Signal threshold for clustering
    info_thresh: float
        Information threshold for clustering
    include_gu: bool
        Include DMS reactivities of Gs and Us for clustering
    inlcude_del: bool
        Include deletions when counting mutations
    min_reads: int
        Minimum number of reads for clustering
    convergence_cutoff: float
        Convergence cutoff for clustering
    min_iterations: int
        Minimum number of iterations for clustering
    num_runs: int
        Number of runs for clustering
    rnastructure_path: str
        Path to RNAstructure, to predict structure and free energy.
    rnastructure_temperature: bool
        Use samples.csv temperature values for RNAstructure.
    rnastructure_fold_args: str
        Options to pass to RNAstructure fold.
    rnastructure_dms: bool
        Use the DMS signal to make predictions with RNAstructure.
    rnastructure_dms_min_unpaired_value: int
        Minimum unpaired value for using the dms signal as an input for RNAstructure.
    rnastructure_dms_max_paired_value: int
        Maximum paired value for using the dms signal as an input for RNAstructure.
    rnastructure_partition: bool
        Use RNAstructure partition function to predict free energy.
    rnastructure_probability: bool
        Use RNAstructure partition function to predict per-base mutation probability.
    


    # -----------------------------------------------------------------------------------------------------------------------

    ## Aggregate
    # -----------------------------------------------------------------------------------------------------------------------
    verbose_print('\naggregation \n------------------')
    for sample in samples_names:
        clustering_file = os.path.join(top_dir, 'output', 'clustering', sample + '.json') if cluster else None
        vect_out_dir = os.path.join(top_dir, 'output', 'vectoring', sample)
        constructs = [c for c in os.listdir(vect_out_dir) if os.path.isdir(os.path.join(vect_out_dir, c))]
        sections = [
            [s for s in os.listdir(os.path.join(vect_out_dir, c)) if os.path.isdir(os.path.join(vect_out_dir, c, s))]
            for c in constructs]
        cli.aggregation.run(
            bv_files=[os.path.join(top_dir, 'output', 'vectoring', sample, c, s) for c in constructs for s in
                      sections[constructs.index(c)]],
            fasta=fasta,
            library=library,
            sample=sample,
            samples=samples,
            clustering_file=clustering_file,
            out_dir=os.path.join(top_dir, 'output', 'aggregation'),
            rnastructure_path=rnastructure_path,
            rnastructure_temperature=rnastructure_temperature,
            rnastructure_fold_args=rnastructure_fold_args,
            rnastructure_dms=rnastructure_dms,
            rnastructure_dms_max_paired_value=rnastructure_dms_max_paired_value,
            rnastructure_dms_min_unpaired_value=rnastructure_dms_min_unpaired_value,
            rnastructure_partition=rnastructure_partition,
            rnastructure_probability=rnastructure_probability,
            coords=coords,
            primers=primers,
            fill=cfill)

    verbose_print('Done!')
    # -----------------------------------------------------------------------------------------------------------------------



'''


if __name__ == "__main__":
    cli()
