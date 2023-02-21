import logging
from datetime import datetime
from enum import IntEnum, StrEnum
import os

import click

from dreem.vector.vector import (MATCH_INT, DELET_INT, INS_5_INT,
                                 INS_3_INT, SUB_N_INT, SUB_A_INT,
                                 SUB_C_INT, SUB_G_INT, SUB_T_INT)
COVER_INT = 255


# System information
CWD = os.getcwd()
if (NUM_CPUS := os.cpu_count()) is None:
    logging.warning("Failed to determine CPU count: defaulting to 1")
    NUM_CPUS = 1


class MateOrientationOption(StrEnum):
    """ Options of mate orientation for alignment with Bowtie2.

    See Bowtie2 manual for full documentation:
    https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
    """
    FR = "fr"
    RF = "rf"
    FF = "ff"


class CountOption(StrEnum):
    """ Options for count-based statistics """
    COVERAGE = "v"
    MATCHES = "w"
    DELETIONS = "d"
    ALL_INSERTIONS = "i"
    INSERTIONS_5PRIME = "i5"
    INSERTIONS_3PRIME = "i3"
    ALL_SUBSTITUTIONS = "s"
    SUBS_TO_ADENINE = "a"
    SUBS_TO_CYTOSINE = "c"
    SUBS_TO_GUANINE = "g"
    SUBS_TO_THYMINE = "t"
    ALL_MUTATIONS = "m"


class CountOptionValue(IntEnum):
    """ Values for count-based statistics """
    COVERAGE = COVER_INT
    MATCHES = MATCH_INT
    DELETIONS = DELET_INT
    ALL_INSERTIONS = INS_5_INT + INS_3_INT
    INSERTIONS_5PRIME = INS_5_INT
    INSERTIONS_3PRIME = INS_3_INT
    ALL_SUBSTITUTIONS = SUB_N_INT
    SUBS_TO_ADENINE = SUB_A_INT
    SUBS_TO_CYTOSINE = SUB_C_INT
    SUBS_TO_GUANINE = SUB_G_INT
    SUBS_TO_THYMINE = SUB_T_INT
    ALL_MUTATIONS = COVER_INT - MATCH_INT



# Input/output options
opt_out_dir = click.option("--out-dir", type=click.Path(file_okay=False),
                           default=os.path.join(CWD, "output"))
opt_temp_dir = click.option("--temp-dir", type=click.Path(file_okay=False),
                            default=os.path.join(CWD, "temp"))
opt_save_temp = click.option("--save-temp/--no-save-temp", type=bool,
                             default=False)

# Resource usage options
opt_parallel = click.option("--parallel/--serial", type=bool, default=True)
opt_max_procs = click.option("--max-procs", type=int, default=NUM_CPUS)

# Experiment and analysis setup options
opt_library = click.option("--library", "-L",
                           type=click.Path(exists=True, dir_okay=False))
opt_samples = click.option("--samples", "-S",
                           type=click.Path(exists=True, dir_okay=False))
opt_rerun = click.option("--rerun/--no-rerun", default=False, type=bool)
opt_resume = click.option("--resume/--no-resume", default=False, type=bool)

# Reference sequence (FASTA) files
arg_fasta = click.argument("fasta",
                           type=click.Path(exists=True, dir_okay=False))

# Sequencing read (FASTQ) files
opt_fastqs = click.option("--fastqs", type=click.Path(dir_okay=False),
                          multiple=True,
                          help="FASTQ file of single-end reads")
opt_fastqi = click.option("--fastqi", type=click.Path(dir_okay=False),
                          multiple=True,
                          help="FASTQ file of interleaved paired reads")
opt_fastq1 = click.option("--fastq1", type=click.Path(dir_okay=False),
                          multiple=True,
                          help="FASTQ file of mate 1 paired-end reads")
opt_fastq2 = click.option("--fastq2", type=click.Path(dir_okay=False),
                          multiple=True,
                          help="FASTQ file of mate 2 paired-end reads")

# Sequencing read (FASTQ/BAM) options
opt_phred_enc = click.option("--phred-enc", "-e", type=int, default=33)
opt_min_phred = click.option("--min-phred", "-q", type=int, default=25)
opt_fastqc = click.option("--fastqc/--no-fastqc", type=bool, default=True)
opt_fastqc_extract = click.option("--fastqc-extract/--no-fastqc-extract",
                                  type=bool, default=True)

# Demultiplexing options
opt_demultiplex = click.option("--demultiplex", "-dx", type=bool,
                               default=False)
opt_barcode_start = click.option("--barcode-start", "-bs", type=int,
                                 default=0)
opt_barcode_length = click.option("--barcode_length", "-bl", type=int,
                                  default=0)
opt_max_barcode_mismatches = click.option("--max_barcode_mismatches", type=int,
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
opt_align_local = click.option("--align-local/--align-e2e", type=bool,
                               default=True)
opt_align_unal = click.option("--align-unal/--align-no_unal", type=bool,
                              default=False)
opt_align_disc = click.option("--align-disc/--align-no_disc", type=bool,
                              default=False)
opt_align_mixed = click.option("--align-mixed/--align-no_mixed", type=bool,
                               default=False)
opt_align_dove = click.option("--align-dove/--align-no_dove", type=bool,
                              default=False)
opt_align_cont = click.option("--align-cont/--align-no_cont", type=bool,
                              default=True)
opt_align_minl = click.option("--align-minl", type=int, default=0)
opt_align_maxl = click.option("--align-maxl", type=int, default=300)
opt_align_score = click.option("--align-score", type=str, default="L,0,0.5")
opt_align_iseed = click.option("--align-iseed", type=str, default="L,1,0.1")
opt_align_lseed = click.option("--align-lseed", type=int, default=12)
opt_align_gbar = click.option("--align-gbar", type=int, default=4)
opt_align_exten = click.option("--align-exten", type=int, default=4)
opt_align_reseed = click.option("--align-reseed", type=int, default=2)
opt_align_pad = click.option("--align-pad", type=int, default=2)
opt_align_orient = click.option("--align-orient",
                                type=click.Choice(tuple(MateOrientationOption),
                                                  case_sensitive=False),
                                default=MateOrientationOption.FR)

# Reference region specification options
opt_coords = click.option("--coords", "-c", type=(str, int, int),
                          multiple=True)
opt_primers = click.option("--primers", "-p", type=(str, str, str),
                           multiple=True)
opt_primer_gap = click.option("--primer_gap", type=int, default=2)
opt_spanall = click.option("--spanall/--no-spanall", type=bool,
                           default=False)

# Mutational profile report files
arg_report = click.argument("report", nargs=-1,
                            type=click.Path(exists=True, dir_okay=False))

# Clustering options
opt_cluster = click.option("--cluster/--no-cluster", type=bool, default=False)
opt_max_clusters = click.option("--max_clusters", type=int, default=3)
opt_min_iter = click.option("--min_iter", type=int, default=100)
opt_signal_thresh = click.option("--signal_thresh", type=float, default=0.005)
opt_info_thresh = click.option("--info_thresh", type=float, default=0.05)
opt_include_gu = click.option("--include_gu/--exclude_gu", type=bool, default=False)
opt_include_del = click.option("--include_del/--exclude_del", type=bool, default=False)
opt_min_reads = click.option("--min_reads", type=int, default=1000)
opt_convergence_cutoff = click.option("--convergence_cutoff", type=float, default=0.5)
opt_num_runs = click.option("--num_runs", type=int, default=10)


# Statistics options
opt_stats_count = click.option("--count", "-c",
                               type=click.Choice(tuple(CountOption),
                                                 case_sensitive=False),
                               multiple=True)
opt_stats_frac = click.option("--frac", "-f",
                              type=click.Choice(tuple(CountOption),
                                                case_sensitive=False),
                              multiple=True)


# Aggregation
RNASTRUCTURE_PATH = None
RNASTRUCTURE_TEMPERATURE = False
RNASTRUCTURE_FOLD_ARGS = None
RNASTRUCTURE_DMS = False
RNASTRUCTURE_DMS_MIN_UNPAIRED_VALUE = 0.04
RNASTRUCTURE_DMS_MAX_PAIRED_VALUE = 0.01
RNASTRUCTURE_PARTITION = False
RNASTRUCTURE_PROBABILITY = False

rnastructure_path = click.option("--rnastructure_path", "-rs", type=click.Path(exists=True))
rnastructure_temperature = click.option("--rnastructure_temperature", "-rst",
                                        type=int, default=310)
rnastructure_fold_args = click.option("--rnastructure_fold_args", "-rsa", type=str)
rnastructure_dms = click.option("--rnastructure_dms", "-rsd", type=bool, help="Use the DMS signal to make predictions with RNAstructure", default=   RNASTRUCTURE_DMS)
rnastructure_dms_min_unpaired_value = click.option("--rnastructure_dms_min_unpaired_value", "-rsdmin", type=int, help="Minimum unpaired value for using the dms signal as an input for RNAstructure", default=RNASTRUCTURE_DMS_MIN_UNPAIRED_VALUE)
rnastructure_dms_max_paired_value = click.option("--rnastructure_dms_max_paired_value", "-rsdmax", type=int, help="Maximum paired value for using the dms signal as an input for RNAstructure", default=RNASTRUCTURE_DMS_MAX_PAIRED_VALUE)
rnastructure_partition = click.option("--rnastructure_partition", "-rspa", type=bool, help="Use RNAstructure partition function to predict free energy", default=RNASTRUCTURE_PARTITION)
rnastructure_probability = click.option("--rnastructure_probability", "-rspr", type=bool, help="Use RNAstructure partition function to predict per-base mutation probability", default=RNASTRUCTURE_PROBABILITY)

# Logging options
opt_verbose = click.option("--verbose", "-v", count=True)
opt_quiet = click.option("--quiet", "-q", count=True)
opt_logfile = click.option("--log",
                           type=click.Path(exists=False, dir_okay=False),
                           default=os.path.join(CWD, datetime.now().strftime(
                               "dreem-log_%Y-%m-%d_%H:%M:%S"
                           )))
