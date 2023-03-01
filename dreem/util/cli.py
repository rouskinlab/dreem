import logging
from datetime import datetime
from enum import Enum, IntEnum
import os
from typing import Any, Callable, Iterable

import click

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
opt_out_dir = click.option("--out-dir", "-o", type=click.Path(file_okay=False),
                           default=os.path.join(CWD, "output"),
                           help="Where to output all finished files")
opt_temp_dir = click.option("--temp-dir", "-t", type=click.Path(file_okay=False),
                            default=os.path.join(CWD, "temp"),
                            help="Where to write all temporary files")
opt_save_temp = click.option("--save-temp/--erase-temp", type=bool,
                             default=False,
                             help=("Whether to save or erase temporary files "
                                   "after the program exits"))

# Resource usage options
opt_parallel = click.option("--parallel/--no-parallel", type=bool, default=True,
                            help="Whether to run multiple jobs in parallel")
opt_max_procs = click.option("--max-procs", type=int, default=NUM_CPUS,
                             help="Maximum number of simultaneous processes")

# Experiment and analysis setup options
opt_library = click.option("--library",
                           type=click.Path(exists=True, dir_okay=False))
opt_samples = click.option("--samples",
                           type=click.Path(exists=True, dir_okay=False))
opt_rerun = click.option("--rerun/--no-rerun", default=False, type=bool,
                         help="Whether to regenerate files that already exist")
opt_resume = click.option("--resume/--no-resume", default=False, type=bool,
                          help=("Whether to use any existing temporary files "
                                "to resume a process that terminated"))

# Reference sequence (FASTA) files
opt_fasta = click.option("--fasta", "-r",
                         type=click.Path(exists=True, dir_okay=False),
                         help="FASTA file of reference sequences")

# Sequencing read (FASTQ) files
opt_fastqs = click.option("--fastqs", "-s",
                          type=click.Path(exists=True, dir_okay=False),
                          multiple=True,
                          help="FASTQ file of single-end reads")
opt_fastqi = click.option("--fastqi", "-i",
                          type=click.Path(exists=True, dir_okay=False),
                          multiple=True,
                          help="FASTQ file of interleaved paired reads")
opt_fastq1 = click.option("--fastq1", "-1",
                          type=click.Path(exists=True, dir_okay=False),
                          multiple=True,
                          help="FASTQ file of mate 1 paired-end reads")
opt_fastq2 = click.option("--fastq2", "-2",
                          type=click.Path(exists=True, dir_okay=False),
                          multiple=True,
                          help="FASTQ file of mate 2 paired-end reads")

# Sequencing read (FASTQ/BAM) options
opt_phred_enc = click.option("--phred-enc", "-e", type=int,
                             default=DEFAULT_PHRED_ENC,
                             help="Phred score encoding in FASTQ/SAM/BAM files")
opt_min_phred = click.option("--min-phred", "-q", type=int,
                             default=DEFAULT_MIN_PHRED,
                             help="Minimum Phred score to use a base call")
opt_fastqc = click.option("--fastqc/--no-fastqc", type=bool, default=True,
                          help="Whether to check quality of FASTQ files")
opt_fastqc_extract = click.option("--fastqc-extract/--fastqc-no-extract",
                                  type=bool, default=True,
                                  help="Whether to unzip FASTQC reports")

# Demultiplexing options
opt_demultiplex = click.option("--demultiplex", "-dx", type=bool,
                               default=False)
opt_barcode_start = click.option("--barcode-start", "-bs", type=int,
                                 default=0)
opt_barcode_length = click.option("--barcode-length", "-bl", type=int,
                                  default=0)
opt_max_barcode_mismatches = click.option("--max_barcode_mismatches", type=int,
                                          default=1)

# Demultiplexed sequencing read (FASTQ) directories
opt_fastqs_dir = click.option("--fastqs-dir", "-S",
                              type=click.Path(exists=True, file_okay=False),
                              multiple=True,
                              help="Demultiplexed FASTQ files of single-end reads")
opt_fastqi_dir = click.option("--fastqi-dir", "-I",
                              type=click.Path(exists=True, file_okay=False),
                              multiple=True,
                              help="Demultiplexed FASTQ files of interleaved paired-end reads")
opt_fastq12_dir = click.option("--fastq12-dir", "-P",
                               type=click.Path(exists=True, file_okay=False),
                               multiple=True,
                               help="Demultiplexed pairs of FASTQ files of mate 1 and mate 2 reads")

# Alignment map (BAM) files
opt_bamf = click.option("--bamf", "-b",
                        type=click.Path(exists=True, dir_okay=False),
                        multiple=True,
                        help="BAM file")
opt_bamd = click.option("--bamd", "-B",
                        type=click.Path(exists=True, file_okay=False),
                        multiple=True,
                        help="Directory of BAM files")

# Adapter trimming options with Cutadapt
opt_cutadapt = click.option("--cut/--no-cut", type=bool, default=True,
                            help="Whether to trim reads with Cutadapt before alignment")
opt_cut_q1 = click.option("--cut-q1", type=int, default=DEFAULT_MIN_PHRED,
                          help="Phred score for read 1 quality trimming")
opt_cut_q2 = click.option("--cut-q2", type=int, default=DEFAULT_MIN_PHRED,
                          help="Phred score for read 2 quality trimming")
opt_cut_g1 = click.option("--cut-g1", type=str, multiple=True,
                          default=(),
                          help="5' adapter for read 1")
opt_cut_a1 = click.option("--cut-a1", type=str, multiple=True,
                          default=(AdapterSequence.ILLUMINA_3P.value,),
                          help="3' adapter for read 1")
opt_cut_g2 = click.option("--cut-g2", type=str, multiple=True,
                          default=(),
                          help="5' adapter for read 2")
opt_cut_a2 = click.option("--cut-a2", type=str, multiple=True,
                          default=(AdapterSequence.ILLUMINA_3P.value,),
                          help="3' adapter for read 2")
opt_cut_o = click.option("--cut-O", type=int, default=6,
                         help="Minimum overlap of read and adapter")
opt_cut_e = click.option("--cut-e", type=float, default=0.1,
                         help="Error tolerance for adapters")
opt_cut_indels = click.option("--cut-indels/--cut-no-indels",
                              type=bool, default=True,
                              help="Whether to allow indels in adapters")
opt_cut_nextseq = click.option("--cut-nextseq/--cut-no-nextseq",
                               type=bool, default=False,
                               help="Whether to trim high-quality Gs from 3' end")
opt_cut_discard_trimmed = click.option(
    "--cut-discard-trimmed/--cut-keep-trimmed",
    type=bool, default=False,
    help="Whether to discard reads in which an adapter was found")
opt_cut_discard_untrimmed = click.option(
    "--cut-discard-untrimmed/--cut-keep-untrimmed",
    type=bool, default=False,
    help="Whether to discard reads in which no adapter was found")
opt_cut_m = click.option("--cut-m", type=int, default=20,
                         help="Discard reads shorter than this length after trimming")

# Alignment options with Bowtie2
opt_bt2_local = click.option("--bt2-local/--bt2-end-to-end", type=bool,
                             default=True,
                             help="Whether to perform local or end-to-end alignment")
opt_bt2_discordant = click.option("--bt2-discordant/--bt2-no-discordant", type=bool,
                                  default=False,
                                  help="Whether to output discordant alignments")
opt_bt2_mixed = click.option("--bt2-mixed/--bt2-no-mixed", type=bool,
                             default=False,
                             help="Whether to align individual mates of unaligned pairs")
opt_bt2_dovetail = click.option("--bt2-dovetail/--bt2-no-dovetail", type=bool,
                                default=False,
                                help="Whether to treat dovetailed mate pairs as concordant")
opt_bt2_contain = click.option("--bt2-contain/--bt2-no-contain", type=bool,
                               default=True,
                               help="Whether to treat nested mate pairs as concordant")
opt_bt2_i = click.option("--bt2-I", type=int, default=0,
                         help="Minimum fragment length for valid paired-end alignments")
opt_bt2_x = click.option("--bt2-X", type=int, default=600,
                         help="Maximum fragment length for valid paired-end alignments")
opt_bt2_score_min = click.option("--bt2-score-min", type=str, default="L,0,0.5",
                                 help="Minimum score for a valid alignment")
opt_bt2_s = click.option("--bt2-i", "bt2_s", type=str, default="L,1,0.1",
                         help="Seed interval")
opt_bt2_l = click.option("--bt2-L", type=int, default=12,
                         help="Seed length")
opt_bt2_gbar = click.option("--bt2-gbar", type=int, default=4,
                            help="Minimum distance of a gap from end of a read")
opt_bt2_d = click.option("--bt2-D", type=int, default=4,
                         help="Maximum number of failed seed extensions")
opt_bt2_r = click.option("--bt2-R", type=int, default=2,
                         help="Maximum number of times to re-seed")
opt_bt2_dpad = click.option("--bt2-dpad", type=int, default=2,
                            help="Width of padding on alignment matrix, to allow gaps")
opt_bt2_orient = click.option("--bt2-orient",
                              type=click.Choice(tuple(op.value for op
                                                      in MateOrientationOption),
                                                case_sensitive=False),
                              default=MateOrientationOption.FR.value,
                              help="Valid orientations of paired-end mates")

# Reference region specification options
opt_coords = click.option("--coords", "-c", type=(str, int, int),
                          multiple=True)
opt_primers = click.option("--primers", "-p", type=(str, str, str),
                           multiple=True)
opt_primer_gap = click.option("--primer-gap", type=int, default=2)
opt_cfill = click.option("--cfill/--no-cfill", type=bool,
                         default=False)

# Vectoring options
opt_batch_size = click.option("--batch-size", "-z", type=float, default=32.0,
                              help=("Maximum size of each batch of mut_vectors, "
                                    "in millions of base calls"))
opt_strict_pairs = click.option("--strict-pairs/--no-strict-pairs", type=bool,
                                default=True,
                                help=("Whether to require that every paired "
                                      "read that maps to a region also have "
                                      "a mate that maps to the region"))
opt_ambindel = click.option("--ambindel/--no-ambindel", type=bool, default=True,
                            help=("Whether to find and label all ambiguous "
                                  "insertions and deletions (improves accuracy "
                                  "but runs slower)"))

# Mutational profile report files
opt_report = click.option("--mp-report", "-r",
                          type=click.Path(exists=True, dir_okay=False),
                          multiple=True)

# Clustering options
opt_cluster = click.option("--cluster/--no-cluster", type=bool, default=False)
opt_max_clusters = click.option("--max-clusters", "-k", type=int, default=3)
opt_min_iter = click.option("--min-iter", "-m", type=int, default=100)
opt_max_iter = click.option("--max-iter", "-x", type=int, default=500)
opt_signal_thresh = click.option("--signal-thresh", "-s", type=float, default=0.005)
opt_info_thresh = click.option("--info-thresh", "-i", type=float, default=0.05)
opt_include_gu = click.option("--include-gu/--exclude-gu", type=bool, default=False)
opt_include_del = click.option("--include-del/--exclude-del", type=bool, default=False)
opt_min_reads = click.option("--min-reads", "-v", type=int, default=1000)
opt_convergence_cutoff = click.option("--convergence-cutoff", "-g", type=float, default=0.5)
opt_num_runs = click.option("--num-runs", "-n", type=int, default=10)

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
rnastructure_dms = click.option("--rnastructure_dms", "-rsd", type=bool,
                                help="Use the DMS signal to make predictions with RNAstructure",
                                default=RNASTRUCTURE_DMS)
rnastructure_dms_min_unpaired_value = click.option("--rnastructure_dms_min_unpaired_value", "-rsdmin", type=int,
                                                   help="Minimum unpaired value for using the dms signal as an input for RNAstructure",
                                                   default=RNASTRUCTURE_DMS_MIN_UNPAIRED_VALUE)
rnastructure_dms_max_paired_value = click.option("--rnastructure_dms_max_paired_value", "-rsdmax", type=int,
                                                 help="Maximum paired value for using the dms signal as an input for RNAstructure",
                                                 default=RNASTRUCTURE_DMS_MAX_PAIRED_VALUE)
rnastructure_partition = click.option("--rnastructure_partition", "-rspa", type=bool,
                                      help="Use RNAstructure partition function to predict free energy",
                                      default=RNASTRUCTURE_PARTITION)
rnastructure_probability = click.option("--rnastructure_probability", "-rspr", type=bool,
                                        help="Use RNAstructure partition function to predict per-base mutation probability",
                                        default=RNASTRUCTURE_PROBABILITY)

# Logging options
opt_verbose = click.option("--verbose", "-v", count=True)
opt_quiet = click.option("--quiet", "-q", count=True)
opt_logfile = click.option("--log",
                           type=click.Path(exists=False, dir_okay=False),
                           default=os.path.join(CWD, datetime.now().strftime(
                               "dreem-log_%Y-%m-%d_%H:%M:%S"
                           )))


class DreemCommand(object):
    """

    """

    def __init__(self,
                 cli_func: Callable,
                 imports: tuple[str, ...],
                 exports: None | str | tuple[str, ...]):
        """
        Parameters
        ----------
        cli_func: callable
            Command line function to wrap
        imports: tuple[str]
            Key(s) to import from the context object into ```kwargs```
            before calling ```cli_func(**kwargs)```. Imported keys
            override existing keys in ```kwargs```.
        exports: None | str | tuple[str]
            Key(s) under which the return value(s) of ```cli_func``` are
            exported to the context object, or None to discard them.
        """
        self._cli_func = cli_func
        self._exports = exports
        self._imports = imports

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

    def _export_result(self, result: Any, kwargs: dict[str, Any]):
        """ Export ```result``` to ```kwargs```. """
        if isinstance(self._exports, str):
            # If _exports is a str, then use it as a key for the result.
            kwargs[self._exports] = result
        elif isinstance(self._exports, tuple):
            # If _exports is a tuple of str, then use each item as a key
            # for an item in result; _exports and result must have the
            # same number of items.
            for key, res in zip(self._exports, result, strict=True):
                kwargs[key] = res
        elif self._exports is not None:
            # The only other valid value of _exports is None, in which
            # case the result is not stored in kwargs.
            raise TypeError(self._exports)

    def __call__(self, ctx_obj: dict[str, Any], **kwargs: Any):
        """
        Wrapper that calls ```cli_func(**kwargs)``` after adding keyword
        arguments from ```ctx_obj``` to ```kwargs```,

        Parameters
        ----------
        ctx_obj: dict[str, Any]
            Object storing attributes of the Click context
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
        result = self._cli_func(**kwargs)
        self._export_result(result, kwargs)
        # Export all keyword arguments to the context object so that
        # subsequent commands can import the values.
        ctx_obj.update(kwargs)
        return result


def dreem_command(imports: tuple[str, ...] = (),
                  exports: None | str | tuple[str, ...] = None):
    def command_decorator(func: Callable):
        return DreemCommand(func, imports, exports)

    return command_decorator
