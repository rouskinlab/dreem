import logging
from datetime import datetime
from enum import StrEnum
import os

import click

# System information
CWD = os.getcwd()
if (NUM_CPUS := os.cpu_count()) is None:
    logging.warning("Failed to determine CPU count: defaulting to 1")
    NUM_CPUS = 1


class ParallelOption(StrEnum):
    """ Options of parallelization.

    BROAD: Process all profiles simultaneously, in parallel.
    DEEP: Process profiles serially, and parallelize within each.
    AUTO: Automatically choose "broad" or "deep" parallelization.
    """
    AUTO = "auto"
    BROAD = "broad"
    DEEP = "deep"


class MateOrientationOption(StrEnum):
    """ Options of mate orientation for alignment with Bowtie2.

    See Bowtie2 manual for full documentation:
    https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
    """
    FR = "fr"
    RF = "rf"
    FF = "ff"


# Input/output options
opt_out_dir = click.option('--out_dir', type=click.Path(file_okay=False),
                           default=os.path.join(CWD, "output"))
opt_temp_dir = click.option('--temp_dir', type=click.Path(file_okay=False),
                            default=os.path.join(CWD, "temp"))

# Resource usage options
opt_parallel = click.option('--parallel',
                            type=click.Choice(tuple(ParallelOption),
                                              case_sensitive=False),
                            default=ParallelOption.AUTO)
opt_max_cpus = click.option('--max_cpus', type=int, default=NUM_CPUS)

# Experiment and analysis setup options
opt_library = click.option('--library',
                           type=click.Path(exists=True, dir_okay=False))
opt_samples = click.option('--samples',
                           type=click.Path(exists=True, dir_okay=False))
opt_rerun = click.option('--rerun/--no-rerun', default=False, type=bool)
opt_resume = click.option('--resume/--no-resume', default=False, type=bool)

# Reference sequence (FASTA) files
arg_fasta = click.argument("fasta",
                           type=click.Path(exists=True, dir_okay=False))

# Sequencing read (FASTQ) files
opt_fastqs = click.option("--fastqs", type=click.Path(), multiple=True,
                          help="FASTQ file of single-end reads")
opt_fastqi = click.option("--fastqi", type=click.Path(), multiple=True,
                          help="FASTQ file of interleaved paired reads")
opt_fastq1 = click.option("--fastq1", type=click.Path(), multiple=True,
                          help="FASTQ file of mate 1 paired-end reads")
opt_fastq2 = click.option("--fastq2", type=click.Path(), multiple=True,
                          help="FASTQ file of mate 2 paired-end reads")

# Sequencing read (FASTQ/BAM) options
opt_phred_enc = click.option("--phred_enc", "-e", type=int, default=33)
opt_min_phred = click.option("--min_phred", "-q", type=int, default=25)

# Demultiplexing options
opt_demultiplex = click.option('--demultiplex', '-dx', type=bool,
                               default=False)
opt_barcode_start = click.option('--barcode_start', '-bs', type=int,
                                 default=0)
opt_barcode_length = click.option('--barcode_length', '-bl', type=int,
                                  default=0)
opt_max_barcode_mismatches = click.option('--max_barcode_mismatches', type=int,
                                          default=1)

# Demultiplexed sequencing read (FASTQ) directories
opt_fastqs_dir = click.option("--fastqs_dir",
                              type=click.Path(exists=True, file_okay=False),
                              multiple=True,
                              help="Directory of single-end FASTQ files")
opt_fastqi_dir = click.option("--fastqi_dir",
                              type=click.Path(exists=True, file_okay=False),
                              multiple=True,
                              help="Directory of interleaved FASTQ files")
opt_fastq12_dir = click.option("--fastq12_dir",
                               type=click.Path(exists=True, file_okay=False),
                               multiple=True,
                               help="Directory of paired-end FASTQ file pairs")

# Alignment map (BAM) files
arg_bams = click.argument("bams", nargs=-1, type=click.Path(exists=True, dir_okay=False))  # path to one or more BAM files

# Adapter trimming options with Cutadapt
opt_trim = click.option("--trim/--no_trim", type=bool, default=True)
opt_trim_minq1 = click.option("--trim_minq1", type=int, default=25)
opt_trim_minq2 = click.option("--trim_minq2", type=int, default=25)
opt_trim_adapt15 = click.option("--trim_adapt15", type=str, multiple=True,
                                default=())
opt_trim_adapt13 = click.option("--trim_adapt13", type=str, multiple=True,
                                default=("AGATCGGAAGAGC",))
opt_trim_adapt25 = click.option("--trim_adapt25", type=str, multiple=True,
                                default=())
opt_trim_adapt23 = click.option("--trim_adapt23", type=str, multiple=True,
                                default=("AGATCGGAAGAGC",))
opt_trim_minover = click.option("--trim_minover", type=int, default=6)
opt_trim_maxerr = click.option("--trim_maxerr", type=float, default=0.1)
opt_trim_indels = click.option("--trim_indels/--trim_no_indels",
                               type=bool, default=True)
opt_trim_nextseq = click.option("--trim_nextseq/--no_trim_nextseq",
                                type=bool, default=False)
opt_trim_xtrim = click.option("--trim_discard_trimmed/--trim_keep_trimmed",
                              type=bool, default=False)
opt_trim_xuntrim = click.option("--trim_discard_untrimmed/--cutadapt_keep_untrimmed",
                                type=bool, default=False)
opt_trim_minlen = click.option("--trim_minlen", type=int, default=20)

# Alignment options with Bowtie2
opt_align_local = click.option("--align_local/--align_e2e", type=bool, default=True)
opt_align_unal = click.option("--align_unal/--align_no_unal", type=bool, default=False)
opt_align_disc = click.option("--align_disc/--align_no_disc", type=bool, default=False)
opt_align_mixed = click.option("--align_mixed/--align_no_mixed", type=bool, default=False)
opt_align_dove = click.option("--align_dove/--align_no_dove", type=bool, default=False)
opt_align_cont = click.option("--align_cont/--align_no_cont", type=bool, default=True)
opt_align_minl = click.option("--align_minl", type=int, default=0)
opt_align_maxl = click.option("--align_maxl", type=int, default=600)
opt_align_score = click.option("--align_score", type=str, default="L,4,0.8")
opt_align_sint = click.option("--align_sint", type=str, default="L,1,0.1")
opt_align_slen = click.option("--align_slen", type=int, default=12)
opt_align_gbar = click.option("--align_gbar", type=int, default=4)
opt_align_exten = click.option("--align_exten", type=int, default=4)
opt_align_reseed = click.option("--align_reseed", type=int, default=0)
opt_align_pad = click.option("--align_pad", type=int, default=2)
opt_align_orient = click.option('--align_orient',
                                type=click.Choice(tuple(MateOrientationOption),
                                                  case_sensitive=False),
                                default=MateOrientationOption.FR)

# Reference region specification options
opt_coords = click.option('--coords', '-c', type=(str, int, int), multiple=True)
opt_primers = click.option('--primers', '-p', type=(str, int, int), multiple=True)
opt_spanall = click.option('--spanall/--no-spanall', type=bool, default=False)

# Clustering options
opt_cluster = click.option('--cluster/--no-cluster', type=bool, default=False)
opt_max_clusters = click.option('--max_clusters', type=int, default=3)
opt_min_iter = click.option('--min_iter', type=int, default=100)
opt_signal_thresh = click.option('--signal_thresh', type=float, default=0.005)
opt_info_thresh = click.option('--info_thresh', type=float, default=0.05)
opt_include_gu = click.option('--include_gu/--exclude_gu', type=bool, default=False)
opt_include_del = click.option('--include_del/--exclude_del', type=bool, default=False)
opt_min_reads = click.option('--min_reads', type=int, default=1000)
opt_convergence_cutoff = click.option('--convergence_cutoff', type=float, default=0.5)
opt_num_runs = click.option('--num_runs', type=int, default=10)

# Aggregation
RNASTRUCTURE_PATH = None
RNASTRUCTURE_TEMPERATURE = False
RNASTRUCTURE_FOLD_ARGS = None
RNASTRUCTURE_DMS = False
RNASTRUCTURE_DMS_MIN_UNPAIRED_VALUE = 0.04
RNASTRUCTURE_DMS_MAX_PAIRED_VALUE = 0.01
RNASTRUCTURE_PARTITION = False
RNASTRUCTURE_PROBABILITY = False

rnastructure_path = click.option('--rnastructure_path', '-rs', type=click.Path(exists=True))
rnastructure_temperature = click.option('--rnastructure_temperature', '-rst',
                                        type=int, default=310)
rnastructure_fold_args = click.option('--rnastructure_fold_args', '-rsa', type=str)
rnastructure_dms = click.option('--rnastructure_dms', '-rsd', type=bool, help='Use the DMS signal to make predictions with RNAstructure', default=   RNASTRUCTURE_DMS)
rnastructure_dms_min_unpaired_value = click.option('--rnastructure_dms_min_unpaired_value', '-rsdmin', type=int, help='Minimum unpaired value for using the dms signal as an input for RNAstructure', default=RNASTRUCTURE_DMS_MIN_UNPAIRED_VALUE)
rnastructure_dms_max_paired_value = click.option('--rnastructure_dms_max_paired_value', '-rsdmax', type=int, help='Maximum paired value for using the dms signal as an input for RNAstructure', default=RNASTRUCTURE_DMS_MAX_PAIRED_VALUE)
rnastructure_partition = click.option('--rnastructure_partition', '-rspa', type=bool, help='Use RNAstructure partition function to predict free energy', default=RNASTRUCTURE_PARTITION)
rnastructure_probability = click.option('--rnastructure_probability', '-rspr', type=bool, help='Use RNAstructure partition function to predict per-base mutation probability', default=RNASTRUCTURE_PROBABILITY)

# Logging options
opt_verbose = click.option('--verbose', '-v', count=True)
opt_quiet = click.option('--quiet', '-q', count=True)
opt_logfile = click.option("--log",
                           type=click.Path(exists=False, dir_okay=False),
                           default=os.path.join(CWD, datetime.now().strftime(
                               "dreem-log_%Y-%m-%d_%H:%M:%S"
                           )))
