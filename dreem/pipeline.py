from collections import defaultdict

from dreem import util
from dreem.util import path
import dreem  # import all the macros from the cli.py file
from dreem.util.cli import *
from dreem.util.reads import FastqUnit


def run(
    top_dir= TOP_DIR,
    fastq1= FASTQ1,
    fastq2= FASTQ2,
    fasta= FASTA,
    library= LIBRARY,
    samples= SAMPLES,
    demultiplexing= True,
    clustering= True,
    primers= PRIMERS,
    coords= COORDS,
    fill= FILL,
    parallel= PARALLEL,
    interleaved= DEFAULT_INTERLEAVED_INPUT,
    barcode_start= BARCODE_START,
    barcode_length= BARCODE_LENGTH,
    max_barcode_mismatches= MAX_BARCODE_MISMATCHES,
    max_clusters= MAX_CLUSTERS,
    signal_thresh = SIGNAL_THRESH,
    info_thresh = INFO_THRESH,
    include_g_u= INCLUDE_G_U,
    include_del = INCLUDE_DEL,
    min_reads= MIN_READS,
    convergence_cutoff= CONVERGENCE_CUTOFF,
    min_iterations= MIN_ITER,
    num_runs= NUM_RUNS,
    max_cpus= MAX_CPUS,
    rnastructure_path= RNASTRUCTURE_PATH,
    rnastructure_temperature= RNASTRUCTURE_TEMPERATURE,
    rnastructure_fold_args= RNASTRUCTURE_FOLD_ARGS,
    rnastructure_dms= RNASTRUCTURE_DMS,
    rnastructure_dms_min_unpaired_value= RNASTRUCTURE_DMS_MIN_UNPAIRED_VALUE,
    rnastructure_dms_max_paired_value= RNASTRUCTURE_DMS_MAX_PAIRED_VALUE,
    rnastructure_partition= RNASTRUCTURE_PARTITION,
    rnastructure_probability= RNASTRUCTURE_PROBABILITY,
    verbose= VERBOSE    
    ):
    """Run DREEM. The input arguments are parsed from the command line. They correspond to the parameters of the functions in the other modules.

    Parameters
    ----------
    top_dir: str
        Path to the output folder.
    fastq1: str
        Path to the fastq file 1.
    fastq2: str
        Path to the fastq file 2.
    fasta: str
        Path to the fasta file.
    library: str
        Path to the library file.
    samples: str
        Path to the samples file.
    
    demultiplexing: bool
        Run demultiplexing.
    clustering: bool
        Run clustering.
        
    primers: str
        coordinates for reference: '-c ref-name first last'
    coords: str
        primers for reference: '-p ref-name fwd rev'
    fill: str
        Fill in coordinates of reference sequences for which neither coordinates nor primers were given (default: no).
    parallel: str
        "Parallelize the processing of mutational PROFILES or READS within each profile, turn parallelization OFF, or AUTO matically choose the parallelization method (default: auto).
    interleaved: bool
        Fastq files are interleaved.
    
    barcode_start: int
        Start position of the barcode.
    barcode_length: int
        Length of the barcode.
    max_barcode_mismatches: int
        Maximum number of mismatches in the barcode.
    
    max_clusters: int
        Maximum number of clusters.
    signal_thresh: float
        Signal threshold.
    info_thresh: float
        Information threshold.
    include_g_u: bool
        Include G-U pairs.
    include_del: bool
        Include deletions.
    min_reads: int
        Minimum number of reads to cluster a reference.
    convergence_cutoff: float
        Convergence cutoff.
    min_iterations: int
        Minimum number of iterations.
    num_runs: int
        Number of runs.
    max_cpus: int
        Number of CPUs.
    
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
    
    verbose: bool
        Verbose output.
    
    """

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

    if demultiplexing:
        verbose_print('\ndemultiplexing \n------------------')

        demult_fqs = dreem.demultiplexing.run(
            fastqs=fastqs, fastqi=fastqi, fastq1=fastq1, fastq2=fastq2,
            phred_enc=phred_enc, fasta=fasta, top_dir=top_dir, library=library,
            max_barcode_mismatches=max_barcode_mismatches)
        verbose_print('demultiplexing done')

        # Rearrange the outputs of demultiplexing
        # from {sample1: {ref1: FastqUnit1, ref2: FastqUnit2, ...},3
        #       sample2: {ref1: FastqUnit3, ref2: FastqUnit4, ...},
        #       ...}
        # to {fastqs_dir:  (SampleDir1, SampleDir2, ...),
        #     fastqi_dir:  (SampleDir3, SampleDir4, ...),
        #     fastq12_dir: (SampleDir5, SampleDir6, ...)}
        # to match the format of align_fqs without demultiplexing.
        align_fqs = defaultdict(set)
        for sample, references in demult_fqs.items():
            for ref, fq_unit in references.items():
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
    bams = dreem.alignment.run(parallel=parallel,
                               max_cpus=max_cpus,
                               top_dir=top_dir,
                               fasta=fasta,
                               phred_enc=phred_enc,
                               **align_fqs)
    # -----------------------------------------------------------------------------------------------------------------------

    ## Vectoring
    # -----------------------------------------------------------------------------------------------------------------------
    verbose_print('\nvectoring \n------------------')
    vector_reports = dreem.vectoring.run(top_dir=top_dir,
                                         bam_dirs=[bam.path.parent
                                                   for bam in bams],
                                         fasta=fasta,
                                         coords=coords,
                                         primers=primers,
                                         fill=fill,
                                         library=library,
                                         parallel=parallel,
                                         max_cpus=max_cpus,
                                         phred_enc=phred_enc,
                                         min_phred=min_phred,
                                         rerun=rerun)
    # -----------------------------------------------------------------------------------------------------------------------

    ## Clustering (optional)
    # -----------------------------------------------------------------------------------------------------------------------
    if clustering:
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
            include_g_u=include_g_u,
            min_reads=min_reads,
            convergence_cutoff=convergence_cutoff,
            num_runs=num_runs,
            n_cpus=max_cpus
        )
    # -----------------------------------------------------------------------------------------------------------------------

    ## Aggregate
    # -----------------------------------------------------------------------------------------------------------------------
    verbose_print('\naggregation \n------------------')
    for sample in samples_names:
        clustering_file = os.path.join(top_dir, 'output', 'clustering', sample + '.json') if clustering else None
        vect_out_dir = os.path.join(top_dir, 'output', 'vectoring', sample)
        references = [c for c in os.listdir(vect_out_dir) if os.path.isdir(os.path.join(vect_out_dir, c))]
        sections = [
            [s for s in os.listdir(os.path.join(vect_out_dir, c)) if os.path.isdir(os.path.join(vect_out_dir, c, s))]
            for c in references]
        dreem.aggregation.run(
            bv_files=[os.path.join(top_dir, 'output', 'vectoring', sample, c, s) for c in references for s in
                      sections[references.index(c)]],
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
            fill=fill)

    verbose_print('Done!')
    # -----------------------------------------------------------------------------------------------------------------------
