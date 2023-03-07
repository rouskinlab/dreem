from datetime import datetime
from enum import Enum, IntEnum
import logging
import os
from typing import Any, Callable, Iterable

from click import Choice, Option, Parameter, Path

from .seq import (MATCH_INT, DELET_INT, INS_5_INT, INS_3_INT,
                  SUB_N_INT, SUB_A_INT, SUB_C_INT, SUB_G_INT, SUB_T_INT,
                  AMBIG_INT)

# System information
CWD = os.getcwd()
if (NUM_CPUS := os.cpu_count()) is None:
    logging.warning("Failed to determine CPU count: defaulting to 1")
    NUM_CPUS = 1

DEFAULT_PHRED_ENC = 33
DEFAULT_MIN_PHRED = 25


class DreemCommandName(Enum):
    """ Commands for DREEM """
    TEST = "test"
    DEMULTIPLEX = "demultiplex"
    ALIGN = "align"
    VECTOR = "vector"
    CLUSTER = "cluster"
    AGGREGATE = "aggregate"
    DRAW = "draw"


class MateOrientationOption(Enum):
    """ Options of mate orientation for alignment with Bowtie2.

    See Bowtie2 manual for full documentation:
    https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
    """
    FR = "fr"
    RF = "rf"
    FF = "ff"


class CountOption(Enum):
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
    COVERAGE = AMBIG_INT
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
    ALL_MUTATIONS = AMBIG_INT - MATCH_INT


class AdapterSequence(Enum):
    """ Adapter sequences """
    ILLUMINA_3P = "AGATCGGAAGAGC"


# Input/output options
opt_out_dir = Option(("--out-dir", "-o"),
                     type=Path(file_okay=False),
                     default=os.path.join(CWD, "output"),
                     help="Where to output all finished files")
opt_temp_dir = Option(("--temp-dir", "-t"),
                      type=Path(file_okay=False),
                      default=os.path.join(CWD, "temp"),
                      help="Where to write all temporary files")
opt_save_temp = Option(("--save-temp/--erase-temp",),
                       type=bool,
                       default=False,
                       help=("Whether to save or erase temporary files "
                             "after the program exits"))

# Resource usage options
opt_parallel = Option(("--parallel/--no-parallel",),
                      type=bool,
                      default=True,
                      help="Whether to run multiple jobs in parallel")
opt_max_procs = Option(("--max-procs",),
                       type=int,
                       default=NUM_CPUS,
                       help="Maximum number of simultaneous processes")

# Experiment and analysis setup options
opt_library = Option(("--library",),
                     type=Path(dir_okay=False),
                     default="",
                     help="Library CSV file")
opt_samples = Option(("--samples",),
                     type=Path(dir_okay=False),
                     default="",
                     help="Samples file")
opt_rerun = Option(("--rerun/--no-rerun",),
                   type=bool,
                   default=False,
                   help="Whether to regenerate files that already exist")
opt_resume = Option(("--resume/--no-resume",),
                    type=bool,
                    default=False,
                    help=("Whether to use any existing temporary files "
                          "to resume a process that terminated"))

# Reference sequence (FASTA) files
opt_fasta = Option(("--fasta", "-f"),
                   type=Path(exists=True, dir_okay=False),
                   help="FASTA file of all reference sequences in the project")

# Sequencing read (FASTQ) files
opt_fastqs = Option(("--fastqs", "-s"),
                    type=Path(exists=True, dir_okay=False),
                    multiple=True,
                    default=(),
                    help="FASTQ file of single-end reads")
opt_fastqi = Option(("--fastqi", "-i"),
                    type=Path(exists=True, dir_okay=False),
                    multiple=True,
                    default=(),
                    help="FASTQ file of interleaved paired reads")
opt_fastq1 = Option(("--fastq1", "-1"),
                    type=Path(exists=True, dir_okay=False),
                    multiple=True,
                    default=(),
                    help="FASTQ file of mate 1 paired-end reads")
opt_fastq2 = Option(("--fastq2", "-2"),
                    type=Path(exists=True, dir_okay=False),
                    multiple=True,
                    default=(),
                    help="FASTQ file of mate 2 paired-end reads")

# Sequencing read (FASTQ/BAM) options
opt_phred_enc = Option(("--phred-enc", "-e"),
                       type=int,
                       default=DEFAULT_PHRED_ENC,
                       help="Phred score encoding in FASTQ/SAM/BAM files")
opt_min_phred = Option(("--min-phred", "-q"),
                       type=int,
                       default=DEFAULT_MIN_PHRED,
                       help="Minimum Phred score to use a base call")
opt_fastqc = Option(("--fastqc/--no-fastqc",),
                    type=bool,
                    default=True,
                    help="Whether to check quality of FASTQ files")
opt_fastqc_extract = Option(("--fastqc-extract/--fastqc-no-extract",),
                            type=bool,
                            default=True,
                            help="Whether to unzip FASTQC reports")

# Demultiplexing options
opt_demultiplex = Option(("--demult-on/--demult-off",),
                         type=bool,
                         default=False,
                         help="Whether to run demultiplexing")
opt_barcode_start = Option(("--barcode-start",),
                           type=int,
                           default=0)
opt_barcode_length = Option(("--barcode-length",),
                            type=int,
                            default=0)
opt_max_barcode_mismatches = Option(("--max_barcode_mismatches",),
                                    type=int,
                                    default=1)

# Demultiplexed sequencing read (FASTQ) directories
opt_fastqs_dir = Option(("--fastqs-dir", "-S"),
                        type=Path(exists=True, file_okay=False),
                        multiple=True,
                        default=(),
                        help="Demultiplexed FASTQ files of single-end reads")
opt_fastqi_dir = Option(("--fastqi-dir", "-I"),
                        type=Path(exists=True, file_okay=False),
                        multiple=True,
                        default=(),
                        help="Demultiplexed FASTQ files of interleaved paired-end reads")
opt_fastq12_dir = Option(("--fastq12-dir", "-P"),
                         type=Path(exists=True, file_okay=False),
                         multiple=True,
                         default=(),
                         help="Demultiplexed pairs of FASTQ files of mate 1 and mate 2 reads")

# Alignment map (BAM) files
opt_bamf = Option(("--bamf", "-b"),
                  type=Path(exists=True, dir_okay=False),
                  multiple=True,
                  default=(),
                  help="BAM file")
opt_bamd = Option(("--bamd", "-B"),
                  type=Path(exists=True, file_okay=False),
                  multiple=True,
                  default=(),
                  help="Directory of BAM files")

# Adapter trimming options with Cutadapt
opt_cutadapt = Option(("--cut/--no-cut",),
                      type=bool,
                      default=True,
                      help="Whether to trim reads with Cutadapt before alignment")
opt_cut_q1 = Option(("--cut-q1",),
                    type=int,
                    default=DEFAULT_MIN_PHRED,
                    help="Phred score for read 1 quality trimming")
opt_cut_q2 = Option(("--cut-q2",),
                    type=int,
                    default=DEFAULT_MIN_PHRED,
                    help="Phred score for read 2 quality trimming")
opt_cut_g1 = Option(("--cut-g1",),
                    type=str,
                    multiple=True,
                    default=(),
                    help="5' adapter for read 1")
opt_cut_a1 = Option(("--cut-a1",),
                    type=str,
                    multiple=True,
                    default=(AdapterSequence.ILLUMINA_3P.value,),
                    help="3' adapter for read 1")
opt_cut_g2 = Option(("--cut-g2",),
                    type=str,
                    multiple=True,
                    default=(),
                    help="5' adapter for read 2")
opt_cut_a2 = Option(("--cut-a2",),
                    type=str,
                    multiple=True,
                    default=(AdapterSequence.ILLUMINA_3P.value,),
                    help="3' adapter for read 2")
opt_cut_o = Option(("--cut-O",),
                   type=int,
                   default=6,
                   help="Minimum overlap of read and adapter")
opt_cut_e = Option(("--cut-e",),
                   type=float,
                   default=0.1,
                   help="Error tolerance for adapters")
opt_cut_indels = Option(("--cut-indels/--cut-no-indels",),
                        type=bool,
                        default=True,
                        help="Whether to allow indels in adapters")
opt_cut_nextseq = Option(("--cut-nextseq/--cut-no-nextseq",),
                         type=bool,
                         default=False,
                         help="Whether to trim high-quality Gs from 3' end")
opt_cut_discard_trimmed = Option(
    ("--cut-discard-trimmed/--cut-keep-trimmed",),
    type=bool,
    default=False,
    help="Whether to discard reads in which an adapter was found")
opt_cut_discard_untrimmed = Option(
    ("--cut-discard-untrimmed/--cut-keep-untrimmed",),
    type=bool,
    default=False,
    help="Whether to discard reads in which no adapter was found")
opt_cut_m = Option(("--cut-m",),
                   type=int,
                   default=20,
                   help="Discard reads shorter than this length after trimming")

# Alignment options with Bowtie2
opt_bt2_local = Option(("--bt2-local/--bt2-end-to-end",),
                       type=bool,
                       default=True,
                       help="Whether to perform local or end-to-end alignment")
opt_bt2_discordant = Option(("--bt2-discordant/--bt2-no-discordant",),
                            type=bool,
                            default=False,
                            help="Whether to output discordant alignments")
opt_bt2_mixed = Option(("--bt2-mixed/--bt2-no-mixed",),
                       type=bool,
                       default=False,
                       help="Whether to align individual mates of unaligned pairs")
opt_bt2_dovetail = Option(("--bt2-dovetail/--bt2-no-dovetail",),
                          type=bool,
                          default=False,
                          help="Whether to treat dovetailed mate pairs as concordant")
opt_bt2_contain = Option(("--bt2-contain/--bt2-no-contain",),
                         type=bool,
                         default=True,
                         help="Whether to treat nested mate pairs as concordant")
opt_bt2_unal = Option(("--bt2-unal/--bt2-no-unal",),
                      type=bool,
                      default=False,
                      help="Whether to output unaligned reads")
opt_bt2_i = Option(("--bt2-I",),
                   type=int,
                   default=0,
                   help="Minimum fragment length for valid paired-end alignments")
opt_bt2_x = Option(("--bt2-X",),
                   type=int,
                   default=600,
                   help="Maximum fragment length for valid paired-end alignments")
opt_bt2_score_min = Option(("--bt2-score-min",),
                           type=str,
                           default="L,0,0.5",
                           help="Minimum score for a valid alignment")
opt_bt2_s = Option(("--bt2-i", "bt2_s"),
                   type=str,
                   default="L,1,0.1",
                   help="Seed interval")
opt_bt2_l = Option(("--bt2-L",),
                   type=int,
                   default=12,
                   help="Seed length")
opt_bt2_gbar = Option(("--bt2-gbar",),
                      type=int,
                      default=4,
                      help="Minimum distance of a gap from end of a read")
opt_bt2_d = Option(("--bt2-D",),
                   type=int,
                   default=4,
                   help="Maximum number of failed seed extensions")
opt_bt2_r = Option(("--bt2-R",),
                   type=int,
                   default=2,
                   help="Maximum number of times to re-seed")
opt_bt2_dpad = Option(("--bt2-dpad",),
                      type=int,
                      default=2,
                      help="Width of padding on alignment matrix, to allow gaps")
opt_bt2_orient = Option(("--bt2-orient",),
                        type=Choice(tuple(op.value for op
                                          in MateOrientationOption),
                                    case_sensitive=False),
                        default=MateOrientationOption.FR.value,
                        help="Valid orientations of paired-end mates")

# Other alignment options
opt_rem_buffer = Option(("--rem-buffer",),
                        type=int,
                        default=65536,
                        help=("Maximum number of reads to hold in memory when "
                              "removing reads with multiple equal alignments"))

# Reference region specification options
opt_coords = Option(("--coords", "-c"),
                    type=(str, int, int),
                    multiple=True,
                    default=(),
                    help=("Reference name, 5' end, and 3' end of a region; "
                          "coordinates are 1-indexed and include both ends"))
opt_primers = Option(("--primers", "-p"),
                     type=(str, str, str),
                     multiple=True,
                     default=(),
                     help=("Reference name, forward primer, and reverse primer "
                           "of a region; reverse primer must be given 5' to 3'"))
opt_primer_gap = Option(("--primer-gap",),
                        type=int,
                        default=2,
                        help=("Number of bases to leave as a gap between the "
                              "end of a primer and the end of the region"))
opt_cfill = Option(("--cfill/--no-cfill",),
                   type=bool,
                   default=False,
                   help=("Whether, for every reference that was not explicitly "
                         "given at least one region (by --initial_coords or "
                         "--primers), to generate coordinates covering the "
                         "entire reference sequence automatically"))

# Vectoring options
opt_batch_size = Option(("--batch-size", "-z"),
                        type=float,
                        default=32.0,
                        help=("Maximum size of each batch of mut_vectors, "
                              "in millions of base calls"))
opt_strict_pairs = Option(("--strict-pairs/--no-strict-pairs",),
                          type=bool,
                          default=True,
                          help=("Whether to require that every paired "
                                "read that maps to a region also have "
                                "a mate that maps to the region"))
opt_ambid = Option(("--ambid/--no-ambid",),
                   type=bool,
                   default=True,
                   help=("Whether to find and label all ambiguous "
                         "insertions and deletions (improves accuracy "
                         "but runs slower)"))

# Mutational profile report files
opt_report = Option(("--mp-report", "-r"),
                    type=Path(exists=True, dir_okay=False),
                    multiple=True,
                    default=())

# Clustering options
opt_cluster = Option(("--cluster-on/--cluster-off",),
                     type=bool,
                     default=False,
                     help="Whether to run clustering")
opt_max_clusters = Option(("--max-clusters", "-k"),
                          type=int,
                          default=3)
opt_min_iter = Option(("--min-iter", "-m"),
                      type=int,
                      default=100)
opt_max_iter = Option(("--max-iter", "-x"),
                      type=int,
                      default=500)
opt_signal_thresh = Option(("--signal-thresh", "-s"),
                           type=float,
                           default=0.005)
opt_info_thresh = Option(("--info-thresh", "-i"),
                         type=float,
                         default=0.05)
opt_include_gu = Option(("--include-gu/--exclude-gu",),
                        type=bool,
                        default=False)
opt_include_del = Option(("--include-del/--exclude-del",),
                         type=bool,
                         default=False)
opt_min_reads = Option(("--min-reads", "-v"),
                       type=int,
                       default=1000)
opt_convergence_cutoff = Option(("--convergence-cutoff", "-g"),
                                type=float,
                                default=0.5)
opt_num_runs = Option(("--num-runs", "-n"),
                      type=int,
                      default=10)

# Statistics options
opt_stats_count = Option(("--count", "-c"),
                         type=Choice(tuple(CountOption),
                                     case_sensitive=False),
                         multiple=True,
                         default=(), )
opt_stats_frac = Option(("--frac", "-f"),
                        type=Choice(tuple(CountOption),
                                    case_sensitive=False),
                        multiple=True,
                        default=(), )

# Aggregation
RNASTRUCTURE_PATH = ''
RNASTRUCTURE_TEMPERATURE = False
RNASTRUCTURE_FOLD_ARGS = ''
RNASTRUCTURE_DMS = False
RNASTRUCTURE_DMS_MIN_UNPAIRED_VALUE = 0.04
RNASTRUCTURE_DMS_MAX_PAIRED_VALUE = 0.01
RNASTRUCTURE_PARTITION = False
RNASTRUCTURE_PROBABILITY = False

sample = Option(("--sample", "-s"),
                type=str,
                default='output')

clustering_file = Option(("--clustering_file", "-cf"),
                            type=Path(exists=True, dir_okay=False),
                            default='')

rnastructure_path = Option(("--rnastructure_path", "-rs"),
                           type=Path(exists=True))
rnastructure_temperature = Option(("--rnastructure_temperature", "-rst"),
                                  type=bool, default=False)
rnastructure_fold_args = Option(("--rnastructure_fold_args", "-rsa"),
                                type=str,
                                default=RNASTRUCTURE_FOLD_ARGS )
rnastructure_dms = Option(("--rnastructure_dms", "-rsd"),
                          type=bool,
                          default=RNASTRUCTURE_DMS,
                          help="Use the DMS signal to make predictions with RNAstructure")
rnastructure_dms_min_unpaired_value = Option(("--rnastructure_dms_min_unpaired_value", "-rsdmin"),
                                             type=int,
                                             default=RNASTRUCTURE_DMS_MIN_UNPAIRED_VALUE,
                                             help="Minimum unpaired value for using the dms signal as an input for RNAstructure")
rnastructure_dms_max_paired_value = Option(("--rnastructure_dms_max_paired_value", "-rsdmax"),
                                           type=int,
                                           default=RNASTRUCTURE_DMS_MAX_PAIRED_VALUE,
                                           help="Maximum paired value for using the dms signal as an input for RNAstructure")
rnastructure_partition = Option(("--rnastructure_partition", "-rspa"),
                                type=bool,
                                default=RNASTRUCTURE_PARTITION,
                                help="Use RNAstructure partition function to predict free energy")
rnastructure_probability = Option(("--rnastructure_probability", "-rspr"),
                                  type=bool,
                                  default=RNASTRUCTURE_PROBABILITY,
                                  help="Use RNAstructure partition function to predict per-base mutation probability")

# Logging options
opt_verbose = Option(("--verbose", "-v"),
                     count=True)
opt_quiet = Option(("--quiet", "-q"),
                   count=True)
opt_logfile = Option(("--log",),
                     type=Path(exists=False, dir_okay=False),
                     default=os.path.join(CWD, datetime.now().strftime(
                         "dreem-log_%Y-%m-%d_%H:%M:%S"
                     )))


def merge_params(*param_lists: list[Parameter]):
    """ Merge lists of Click parameters, dropping duplicates. """
    params = list()
    names = set()
    for param_list in param_lists:
        for param in param_list:
            if param.name not in names:
                params.append(param)
                names.add(param.name)
    return params


class DreemCommand(object):
    """

    """

    def __init__(self,
                 cli_func: Callable,
                 imports: tuple[str, ...],
                 exports: tuple[str, ...],
                 result_key: str | None):
        """
        Parameters
        ----------
        cli_func: callable
            Command line function to wrap
        imports: tuple[str]
            Key(s) to import from the context object into ```kwargs```
            before calling ```cli_func(**kwargs)```. Imported keys
            override existing keys in ```kwargs```.
        exports: tuple[str]
            Key(s) to export from ```kwargs``` to the context object.
            Exported keys override existing keys in ```ctx_obj```.
        result_key: str | None
            Key under which to export the result of the function, or
            None to discard the result.
        """
        self._cli_func = cli_func
        self._exports = exports
        self._imports = imports
        self._result_key = result_key

    @staticmethod
    def _update_keys(updated: dict[str, Any],
                     updater: dict[str, Any],
                     keys: Iterable[str]):
        """ For every key in ```keys```, if ```updater``` contains the
        key, then add its key-value pair to ```updated```, overriding
        any existing key-value pairs with the same key names in
        ```updated```; otherwise, ignore the key. """
        for key in keys:
            try:
                updated[key] = updater[key]
            except KeyError:
                pass

    def __call__(self, ctx_obj: dict[str, Any], *args: Any, **kwargs: Any):
        """
        Wrapper that calls ```cli_func(*args, **kwargs)``` after adding
        keyword arguments from ```ctx_obj``` to ```kwargs```.

        Parameters
        ----------
        ctx_obj: dict[str, Any]
            Object storing attributes of the Click context
        *args: Any:
            Positional arguments to pass to ```cli_func```
        **kwargs: Any:
            Keyword arguments to pass to ```cli_func```

        Return
        ------
        Any
            Return value of ```cli_func(**kwargs)```
        """
        # Make shallow copy of kwargs so the outer scope can assume that
        # kwargs is not modified.
        kwargs = kwargs.copy()
        # Import selected keyword arguments from the context object into
        # the dictionary of keyword arguments.
        self._update_keys(kwargs, ctx_obj, self._imports)
        # Call the function and optionally store its return value(s).
        result = self._cli_func(*args, **kwargs)
        if self._result_key is not None:
            kwargs[self._result_key] = result
        # Export all keyword arguments to the context object so that
        # subsequent commands can import the values.
        self._update_keys(ctx_obj, kwargs, self._exports)
        return result


def dreem_command(imports: tuple[str, ...] = (),
                  exports: tuple[str, ...] = (),
                  result_key: str | None = None):
    def command_decorator(func: Callable):
        return DreemCommand(func, imports, exports, result_key)

    return command_decorator
