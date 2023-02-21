import logging
from collections import defaultdict

import dreem
from dreem.util.logio import set_verbosity


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
        max_procs: int,
        phred_enc: int,
        min_phred: int,
        fastqc: bool,
        fastqc_extract: bool,
        trim: bool,
        trim_minq1: int,
        trim_minq2: int,
        trim_adapt15: tuple[str],
        trim_adapt13: tuple[str],
        trim_adapt25: tuple[str],
        trim_adapt23: tuple[str],
        trim_minover: int,
        trim_maxerr: float,
        trim_indels: bool,
        trim_nextseq: bool,
        trim_discard_trimmed: bool,
        trim_discard_untrimmed: bool,
        trim_minlen: int,
        align_local: bool,
        align_unal: bool,
        align_disc: bool,
        align_mixed: bool,
        align_dove: bool,
        align_cont: bool,
        align_score: str,
        align_minl: int,
        align_maxl: int,
        align_gbar: int,
        align_lseed: int,
        align_iseed: str,
        align_exten: int,
        align_reseed: int,
        align_pad: int,
        align_orient: str,
        coords: tuple[tuple[str, int, int], ...],
        primers: tuple[tuple[str, str, str], ...],
        primer_gap: int,
        spanall: bool,
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
    max_procs: int (≥ 1; default: os.cpu_count())
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
    spanall: bool (default: False)
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
    
    verbose: int [0, 2]
        0 (): Log warnings and errors
        1 (-v): Also log program status updates
        2 (-vv): Also log detailed information about program operation
    quiet: int [0, 2]
        0 (): Log warnings and errors
        1 (-q): Suppress warnings
        2 (-qq): Suppress non-critical error messages (NOT recommended)
    """

    set_verbosity(verbose, quiet)

    if max_procs < 1:
        logging.warning("Max CPUs must be ≥ 1: setting to 1")
        max_procs = 1


    def verbose_print(*args):
        print(*args) if verbose else None

    # Run DREEM
    verbose_print("""

    ========================================

                RUNNING   DREEM

    ========================================

    """)
    # -----------------------------------------------------------------------------------------------------------------------

    fastq_args = dict(fastqs=fastqs, fastqi=fastqi,
                      fastq1=fastq1, fastq2=fastq2)

    ## Demultiplexing (optional)
    # -----------------------------------------------------------------------------------------------------------------------

    if demultiplex:
        verbose_print('\ndemultiplexing \n------------------')

        demult_fqs = dreem.demultiplexing.run(
            fastqs=fastqs, fastqi=fastqi, fastq1=fastq1, fastq2=fastq2,
            phred_enc=phred_enc, fasta=fasta, top_dir=out_dir, library=library,
            max_barcode_mismatches=max_barcode_mismatches)
        verbose_print('demultiplexing done')

        # Rearrange the outputs of demultiplexing
        # from {sample1: {ref1: FastqUnit1, ref2: FastqUnit2, ...},
        #       sample2: {ref1: FastqUnit3, ref2: FastqUnit4, ...},
        #       ...}
        # to {fastqs_dir:  (SampleDir1, SampleDir2, ...),
        #     fastqi_dir:  (SampleDir3, SampleDir4, ...),
        #     fastq12_dir: (SampleDir5, SampleDir6, ...)}
        # to match the format of align_fqs without demultiplexing.
        align_fqs = defaultdict(set)
        for sample, constructs in demult_fqs.items():
            for ref, fq_unit in constructs.items():
                for fq_key, fq_path in fq_unit.inputs.items():
                    # Add the parent directory (the sample directory) to the set
                    align_fqs[fq_key].add(str(fq_path.path.parent))
        align_fqs = {fq_key: tuple(fq_paths)
                     for fq_key, fq_paths in align_fqs.items()}
    else:
        align_fqs = fastq_args
    # -----------------------------------------------------------------------------------------------------------------------

    ## Alignment: 
    # -----------------------------------------------------------------------------------------------------------------------
    bams = dreem.alignment.run(**align_fqs,
                               phred_enc=phred_enc,
                               refset_file=fasta,
                               out_dir=out_dir,
                               temp_dir=temp_dir,
                               save_temp=save_temp,
                               rerun=rerun,
                               resume=resume,
                               parallel=parallel,
                               max_procs=max_procs,
                               fastqc=fastqc,
                               fastqc_extract=fastqc_extract,
                               trim=trim,
                               trim_minq1=trim_minq1,
                               trim_minq2=trim_minq2,
                               trim_adapt15=trim_adapt15,
                               trim_adapt13=trim_adapt13,
                               trim_adapt25=trim_adapt25,
                               trim_adapt23=trim_adapt23,
                               trim_minover=trim_minover,
                               trim_maxerr=trim_maxerr,
                               trim_indels=trim_indels,
                               trim_nextseq=trim_nextseq,
                               trim_discard_trimmed=trim_discard_trimmed,
                               trim_discard_untrimmed=trim_discard_untrimmed,
                               trim_minlen=trim_minlen,
                               align_local=align_local,
                               align_unal=align_unal,
                               align_disc=align_disc,
                               align_mixed=align_mixed,
                               align_dove=align_dove,
                               align_cont=align_cont,
                               align_score=align_score,
                               align_minl=align_minl,
                               align_maxl=align_maxl,
                               align_gbar=align_gbar,
                               align_lseed=align_lseed,
                               align_iseed=align_iseed,
                               align_exten=align_exten,
                               align_reseed=align_reseed,
                               align_pad=align_pad,
                               align_orient=align_orient)
    # -----------------------------------------------------------------------------------------------------------------------

    ## Vectoring
    # -----------------------------------------------------------------------------------------------------------------------
    verbose_print('\nvectoring \n------------------')
    vector_reports = dreem.vectoring.run(out_dir=out_dir,
                                         temp_dir=temp_dir,
                                         bams=bams,
                                         fasta=fasta,
                                         coords=coords,
                                         primers=primers,
                                         primer_gap=primer_gap,
                                         spanall=spanall,
                                         library=library,
                                         parallel=parallel,
                                         max_procs=max_procs,
                                         phred_enc=phred_enc,
                                         min_phred=min_phred,
                                         rerun=rerun,
                                         resume=resume,
                                         save_temp=save_temp)
    # -----------------------------------------------------------------------------------------------------------------------

    ## Clustering (optional)
    # -----------------------------------------------------------------------------------------------------------------------
    if cluster:
        verbose_print('\nclustering \n------------------')

        report_files = tuple(report.pathstr for report in vector_reports)
        dreem.clustering.run(
            report_files=report_files,
            out_dir=top_dir,
            max_clusters=max_clusters,
            min_iter=min_iter,
            signal_thresh=signal_thresh,
            info_thresh=info_thresh,
            include_del=include_del,
            include_gu=include_gu,
            min_reads=min_reads,
            convergence_cutoff=convergence_cutoff,
            num_runs=num_runs,
            n_cpus=max_procs
        )
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
        dreem.aggregation.run(
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
            fill=spanall)

    verbose_print('Done!')
    # -----------------------------------------------------------------------------------------------------------------------
