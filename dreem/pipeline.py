import yaml
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from dreem import util
import dreem # import all the macros from the cli.py file
from dreem.util.cli import *


def run(fasta:str, fastq1:str, fastq2:str=FASTQ2, library:str=LIBRARY, samples:str=SAMPLES, out_dir:str=OUT_DIR, 
        demultiplexing:bool=DEFAULT_DEMULTIPLEXED, clustering:bool=CLUSTERING,
        primers:str=PRIMERS, coords:str=COORDS, fill:bool=FILL, parallel:bool=PARALLEL, interleaved:bool=DEFAULT_INTERLEAVED_INPUT,
        barcode_start:int=BARCODE_START, barcode_length:int=BARCODE_LENGTH, max_barcode_mismatches:int=MAX_BARCODE_MISMATCHES,
        max_clusters:int=MAX_CLUSTERS, signal_thresh:float=SIGNAL_THRESH, info_thresh:float=INFO_THRESH, include_g_u:bool=INCLUDE_G_U, include_del:bool=INCLUDE_DEL, min_reads:int=MIN_READS, convergence_cutoff:float=CONVERGENCE_CUTOFF, min_iter:int=MIN_ITER, num_runs:int=NUM_RUNS, n_cpus:int=N_CPUS,
        rnastructure_path:str=RNASTRUCTURE_PATH, rnastructure_temperature:bool=RNASTRUCTURE_TEMPERATURE, rnastructure_fold_args:str=RNASTRUCTURE_FOLD_ARGS, rnastructure_dms:bool=RNASTRUCTURE_DMS, rnastructure_dms_min_unpaired_value:int=RNASTRUCTURE_DMS_MIN_UNPAIRED_VALUE, rnastructure_dms_max_paired_value:int=RNASTRUCTURE_DMS_MAX_PAIRED_VALUE, rnastructure_partition:bool=RNASTRUCTURE_PARTITION, rnastructure_probability:bool=RNASTRUCTURE_PROBABILITY, poisson:bool=POISSON,verbose=VERBOSE):
    """Run DREEM. The input arguments are parsed from the command line. They correspond to the parameters of the functions in the other modules.

    Parameters
    ----------
    out_dir: str
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
    inlcude_del: bool
        Include deletions.
    min_reads: int
        Minimum number of reads to cluster a reference.
    convergence_cutoff: float
        Convergence cutoff.
    min_iterations: int
        Minimum number of iterations.
    num_runs: int
        Number of runs.
    n_cpus: int
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
    poisson: bool
        Predict Poisson confidence intervals.
    
    verbose: bool
        Verbose output.
    
    """
    def verbose_print(*args):
        print(*args) if verbose else None

    # make output and temp folders
    for folder in ['output', 'temp']:
        os.makedirs(os.path.join(out_dir, folder), exist_ok=True)
        
    # sort fast pairs
    fastq1, fastq2, samples_names = util.util.sort_fastq_pairs(fastq1, fastq2)
    
    # Run DREEM
    verbose_print("""

    ========================================

                RUNNING   DREEM

    ========================================

    """)
    # -----------------------------------------------------------------------------------------------------------------------
    
    ## Demultiplexing (optional)
    # -----------------------------------------------------------------------------------------------------------------------
    if demultiplexing:
        verbose_print('\ndemultiplexing \n------------------')
        
        dreem.demultiplexing.run(out_dir = os.path.join(out_dir, 'output', 'demultiplexing'),
                           fastq1 = fastq1, 
                           fastq2 = fastq2,
                           fasta = fasta,
                           library=library,
                           interleaved=interleaved,
                           barcode_length = barcode_length,
                           barcode_start = barcode_start,
                           max_barcode_mismatches = max_barcode_mismatches)
        verbose_print('demultiplexing done')
    # -----------------------------------------------------------------------------------------------------------------------
    
    ## Alignment: 
    # -----------------------------------------------------------------------------------------------------------------------
    verbose_print('\nalignment \n----------------')
    if demultiplexing:
        fastq1, fastq2, samples_names = [], [], []
        for sample in os.listdir(os.path.join(out_dir, 'output', 'demultiplexing')):
            path = os.path.join(os.path.join(out_dir, 'output', 'demultiplexing'), sample)
            fastq1 = fastq1 + [os.path.join(path, f) for f in os.listdir(path) if f.endswith('_R1.fastq')]
            fastq2 = fastq2 + [os.path.join(path, f) for f in os.listdir(path) if f.endswith('_R2.fastq')]
            samples_names.append(samples_names)

    for idx,(f1, f2, sample) in enumerate(zip(fastq1, fastq2, samples_names)):
        verbose_print('Aligning this fastq pair: ', '\n   ',f1, '\n   ',f2)
        dreem.alignment.run(
                        out_dir=os.path.join(out_dir),#, 'output','alignment'),
                        fasta=fasta,
                        fastq=f1,
                        fastq2=f2,
                        demultiplexed=demultiplexing
                        )
    # -----------------------------------------------------------------------------------------------------------------------
    
    ## Vectoring
    # -----------------------------------------------------------------------------------------------------------------------
    verbose_print('\nvectoring \n------------------')
    for sample in samples_names:
        path_to_bam = os.path.join(out_dir, 'output', 'alignment', sample)
        bam_files = [os.path.join(path_to_bam, f) for f in os.listdir(path_to_bam) if f.endswith('.bam')]
        dreem.vectoring.run(
                       out_dir= os.path.join(out_dir, 'output'), #TODO
                       bam_dirs= bam_files,
                       fasta=fasta,
                       coords=coords,
                       primers=primers,
                       fill=fill,
                       library=library,
                       parallel=parallel
                       )
    # -----------------------------------------------------------------------------------------------------------------------

    ## Clustering (optional)
    # -----------------------------------------------------------------------------------------------------------------------
    if clustering:
        verbose_print('\nclustering \n------------------')
        
        for sample in samples_names:
            dreem.clustering.run(
                           input_dir = os.path.join(out_dir, 'output', 'vectoring'), 
                           out_dir = os.path.join(out_dir, 'output', 'clustering'),
                           max_clusters = max_clusters,
                           min_iter = min_iter,
                           signal_thresh = signal_thresh,
                           info_thresh = info_thresh,
                           include_del = include_del,
                           include_g_u = include_g_u,
                           min_reads = min_reads,
                           convergence_cutoff = convergence_cutoff,
                           num_runs = num_runs,
                           n_cpus = n_cpus
                           )
    # -----------------------------------------------------------------------------------------------------------------------

    ## Aggregate
    # -----------------------------------------------------------------------------------------------------------------------
    verbose_print('\naggregation \n------------------')
    for sample in samples_names:
        clustering_file = os.path.join(out_dir, 'output', 'clustering', sample+'.json') if clustering else None
        vect_out_dir = os.path.join(out_dir, 'output', 'vectoring', sample)
        constructs = [c for c in os.listdir(vect_out_dir) if os.path.isdir(os.path.join(vect_out_dir, c))]
        sections = [[s for s in os.listdir(os.path.join(vect_out_dir, c)) if os.path.isdir(os.path.join(vect_out_dir, c, s))] for c in constructs]
        dreem.aggregation.run(
                        bv_files= [os.path.join(out_dir,'output','vectoring', sample, c, s) for c in constructs for s in sections[constructs.index(c)]],
                        fasta=fasta,
                        library= library,
                        sample = sample, 
                        samples= samples,
                        clustering_file=clustering_file,
                        out_dir= os.path.join(out_dir, 'output', 'aggregation'),
                        rnastructure_path=rnastructure_path,
                        rnastructure_temperature=rnastructure_temperature,
                        rnastructure_fold_args=rnastructure_fold_args,
                        rnastructure_dms=rnastructure_dms,
                        rnastructure_dms_max_paired_value=rnastructure_dms_max_paired_value,
                        rnastructure_dms_min_unpaired_value=rnastructure_dms_min_unpaired_value,
                        rnastructure_partition=rnastructure_partition,
                        rnastructure_probability=rnastructure_probability,
                        poisson=poisson,
                        coords=coords,
                        primers=primers,
                        fill=fill)

    verbose_print('Done!')
    # -----------------------------------------------------------------------------------------------------------------------
