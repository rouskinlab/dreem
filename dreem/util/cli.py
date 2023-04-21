from datetime import datetime
from enum import Enum
import logging
import os

from click import Choice, Option, Parameter, Path

# System information
CWD = os.getcwd()
if (NUM_CPUS := os.cpu_count()) is None:
    logging.warning("Failed to determine CPU count: defaulting to 1")
    NUM_CPUS = 1

DEFAULT_PHRED_ENC = 33
DEFAULT_MIN_PHRED = 25


class MateOrientationOption(Enum):
    """ Options of mate orientation for alignment with Bowtie2.

    See Bowtie2 manual for full documentation:
    https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
    """
    FR = "fr"
    RF = "rf"
    FF = "ff"


class AdapterSequence(Enum):
    """ Adapter sequences """
    ILLUMINA_3P = "AGATCGGAAGAGC"


# Input/output options
opt_out_dir = Option(("--out-dir",),
                     type=Path(file_okay=False),
                     default=os.path.join(".", "output"),
                     help="Where to output all finished files")
opt_temp_dir = Option(("--temp-dir",),
                      type=Path(file_okay=False),
                      default=os.path.join(".", "temp"),
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

# Reference sequence (FASTA) files
opt_fasta = Option(("--fasta",),
                   type=Path(exists=True, dir_okay=False),
                   help="FASTA file of all reference sequences in the project")

# Sequencing read (FASTQ) files
opt_fastqs = Option(("--fastqs",),
                    type=Path(exists=True, dir_okay=False),
                    multiple=True,
                    default=(),
                    help="FASTQ files of single-end reads")
opt_fastqi = Option(("--fastqi",),
                    type=Path(exists=True, dir_okay=False),
                    multiple=True,
                    default=(),
                    help="FASTQ files of interleaved paired reads")
opt_fastq1 = Option(("--fastq1",),
                    type=Path(exists=True, dir_okay=False),
                    multiple=True,
                    default=(),
                    help="FASTQ files of mate 1 paired-end reads")
opt_fastq2 = Option(("--fastq2",),
                    type=Path(exists=True, dir_okay=False),
                    multiple=True,
                    default=(),
                    help="FASTQ files of mate 2 paired-end reads")

# Sequencing read (FASTQ/BAM) options
opt_phred_enc = Option(("--phred-enc",),
                       type=int,
                       default=DEFAULT_PHRED_ENC,
                       help="Phred score encoding in FASTQ/SAM/BAM files")
opt_min_phred = Option(("--min-phred",),
                       type=int,
                       default=DEFAULT_MIN_PHRED,
                       help="Minimum Phred score to use a base call")
opt_fastqc = Option(("--fastqc/--no-fastqc",),
                    type=bool,
                    default=True,
                    help="Whether to check quality of FASTQ files")
opt_qc_extract = Option(("--qc-extract/--qc-no-extract",),
                        type=bool,
                        default=False,
                        help="Whether to unzip FASTQC reports")

# Demultiplexing options

opt_demultiplex = Option(("--demult-on/--demult-off",),
                         type=bool,
                         default=False,
                         help="Whether to run demultiplexing")

opt_parallel_demultiplexing = Option(("--parallel-demultiplexing",),
                                     type=bool,
                                     default=False,
                                     help="Whether to run demultiplexing at maximum speed by submitting multithreaded grep functions")

opt_clipped_demultiplexing = Option(("--clipped",),
                                    type=int,
                                    default=0,
                                    help="Designates the amount of clipped patterns to search for in the sample, will raise compution time")

opt_mismatch_tolerence = Option(("--mismatch-tolerence",),
                                type=int,
                                default=0,
                                help="Designates the allowable amount of mismatches allowed in a string and still be considered a valid pattern find. \
                            will increase non-parallel computation at a factorial rate. use caution going above 2 mismatches. does not apply to clipped sequences.")

opt_index_tolerence = Option(("--index-tolerance",),
                             type=int,
                             default=0,
                             help="Designates the allowable amount of distance you allow the pattern to be found in a read from the reference index")

opt_barcode_start = Option(("--barcode-start",),
                           type=int,
                           default=0,
                           help="index of start of barcode")
opt_barcode_length = Option(("--barcode-length",),
                            type=int,
                            default=0,
                            help="length of barcode")
opt_demulti_overwrite = Option(("--demulti-overwrite",),
                               type=bool,
                               default=False,
                               help="desiginates whether to overwrite the grepped fastq. should only be used if changing setting on the same sample")

# Demultiplexed sequencing read (FASTQ) directories
opt_fastqs_dir = Option(("--fastqs-dir",),
                        type=Path(exists=True, file_okay=False),
                        multiple=True,
                        default=(),
                        help="Directory containing demultiplexed FASTQ files of single-end reads from one sample")
opt_fastqi_dir = Option(("--fastqi-dir",),
                        type=Path(exists=True, file_okay=False),
                        multiple=True,
                        default=(),
                        help="Directory containing demultiplexed FASTQ files of interleaved paired-end reads from one sample")
opt_fastq12_dir = Option(("--fastq12-dir",),
                         type=Path(exists=True, file_okay=False),
                         multiple=True,
                         default=(),
                         help="Directory containing demultiplexed pairs of FASTQ files of mate 1 and mate 2 reads from one sample")

# Alignment map (BAM) files
opt_bamf = Option(("--bamf",),
                  type=Path(exists=True, dir_okay=False),
                  multiple=True,
                  default=(),
                  help="BAM file")
opt_bamd = Option(("--bamd",),
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
                           default="L,1,0.5",
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

# Reference section specification options
opt_coords = Option(("--coords", "-c"),
                    type=(str, int, int),
                    multiple=True,
                    default=(),
                    help=("Reference name, 5' end, and 3' end of a section; "
                          "coordinates are 1-indexed and include both ends"))
opt_primers = Option(("--primers", "-p"),
                     type=(str, str, str),
                     multiple=True,
                     default=(),
                     help=("Reference name, forward primer, and reverse primer "
                           "of a section; reverse primer must be given 5' to 3'"))
opt_primer_gap = Option(("--primer-gap",),
                        type=int,
                        default=2,
                        help=("Number of bases to leave as a gap between the "
                              "end of a primer and the end of the section"))
opt_autosect = Option(("--autosect/--no-autosect",),
                      type=bool,
                      default=False,
                      help=("Whether, for every reference that was not explicitly "
                            "given at least one section (by --initial_coords or "
                            "--primers), to generate coordinates covering the "
                            "entire reference sequence automatically"))

# Vectoring options
opt_batch_size = Option(("--batch-size", "-z"),
                        type=float,
                        default=32.0,
                        help=("Maximum size of each batch of mutation vectors, "
                              "in millions of base calls"))
opt_strict_pairs = Option(("--strict-pairs/--no-strict-pairs",),
                          type=bool,
                          default=True,
                          help=("Whether to require that every paired "
                                "read that maps to a section also have "
                                "a mate that maps to the section"))
opt_ambid = Option(("--ambid/--no-ambid",),
                   type=bool,
                   default=True,
                   help=("Whether to find and label all ambiguous "
                         "insertions and deletions (improves accuracy "
                         "but runs slower)"))

# Mutational profile report files
opt_report = Option(("--mp-report",),
                    type=Path(exists=True, dir_okay=False),
                    multiple=True,
                    default=(),
                    help="Path to the bit vector folder or list of paths to the bit vector folders.")

# Clustering options
opt_cluster = Option(("--clust/--cluster-off",),
                     type=bool,
                     default=False,
                     help="Whether to run clustering")
opt_max_clusters = Option(("--max-clusters",),
                          type=int,
                          default=3,
                          help='Maximum number of clusters.')
opt_min_iter = Option(("--min-iter",),
                      type=int,
                      default=100,
                      help='Minimum number of iteration before checking convergence of EM.')
opt_max_iter = Option(("--max-iter",),
                      type=int,
                      default=500,
                      help='Maximum number of iteration before stopping EM.')
opt_signal_thresh = Option(("--signal-thresh",),
                           type=float,
                           default=0.005,
                           help='Minimum Mutation fraction to keep a base.')
opt_info_thresh = Option(("--info-thresh",),
                         type=float,
                         default=0.05)
opt_include_gu = Option(("--include-gu/--exclude-gu",),
                        type=bool,
                        default=False,
                        help='Whether to include G and U bases in reads.')
opt_polya_max = Option(("--polya-max",),
                       type=int,
                       default=4,
                       help='Maximum length of poly(A) sequences to include.')
opt_include_del = Option(("--include-del/--exclude-del",),
                         type=bool,
                         default=False,
                         help='Whether to include deletions in reads.')
opt_min_reads = Option(("--min-reads",),
                       type=int,
                       default=1000,
                       help='Minimum number of reads to start clustering.')
opt_convergence_cutoff = Option(("--convergence-cutoff",),
                                type=float,
                                default=0.5,
                                help='Minimum difference between the log-likelihood of two consecutive iterations to stop EM.')
opt_num_runs = Option(("--num-runs", "-n"),
                      type=int,
                      default=10,
                      help='Number of time to run the clustering algorithm.')

# Aggregation

opt_bv_files = Option(("--bv-files", "-bv"),
                      type=Path(exists=True, dir_okay=True, file_okay=False),
                      multiple=True,
                      default=(),
                      help="Tuple of paths. Give the path to the sample folder to process every section. Give the path to a report to process a single section.")

opt_clustering_file = Option(("--clustering-file", "-cf"),
                             # type=Path(exists=True, dir_okay=False),
                             default='',
                             help="Path to the json clustering file from dreem clustering.")

opt_rnastructure_path = Option(("--rnastructure-path", "-rs"),
                               type=Path(),
                               help='Path to the RNAstructure executable folder (e.g. /home/user/RNAstructure/exe/). Use this option if RNAstructure is not in your PATH.',
                               default='')
opt_rnastructure_use_temp = Option(("--rnastructure-use-temp", "-rst"),
                                   type=bool, default=False,
                                   help='Use the temperature signal to make predictions with RNAstructure')
opt_rnastructure_fold_args = Option(("--rnastructure-fold-args", "-rsa"),
                                    type=str,
                                    default='',
                                    help="Additional arguments to pass to RNAstructure's Fold command")

opt_rnastructure_use_dms = Option(("--rnastructure-use-dms", "-rsd"),
                                  type=bool,
                                  default=False,
                                  help="Use the DMS signal to make predictions with RNAstructure")
opt_rnastructure_dms_min_unpaired_value = Option(("--rnastructure-dms-min-unpaired-value", "-rsdmin"),
                                                 type=int,
                                                 default=0.07,
                                                 help="Minimum unpaired value for using the dms signal as an input for RNAstructure")
opt_rnastructure_dms_max_paired_value = Option(("--rnastructure-dms-max-paired-value", "-rsdmax"),
                                               type=int,
                                               default=0.01,
                                               help="Maximum paired value for using the dms signal as an input for RNAstructure")
opt_rnastructure_deltag_ensemble = Option(("--rnastructure-deltag-ensemble", "-rspa"),
                                          type=bool,
                                          default=False,
                                          help="Use RNAstructure partition function to predict free energy")
opt_rnastructure_probability = Option(("--rnastructure_probability", "-rspr"),
                                      type=bool,
                                      default=False,
                                      help="Use RNAstructure partition function to predict per-base pairing probability")

# Drawing options

opt_draw_input = Option(("--inpt",),
                        multiple=True,
                        type=Path(exists=True, dir_okay=False),
                        default=(),
                        help="Path to a dreem output format file. Can be specified multiple times.")

opt_section = Option(("--section",),
                        multiple=True,
                        default=('full',),
                        help="Section to draw. Can be specified multiple times.")

opt_flat = Option(("--flat/--no-flat",),
                  is_flag=True,
                  default=True,
                  help="Flatten the output folder structure. This names your files [reference]__[section]__[plot_name].html")

opt_mutation_fraction = Option(("--mutation_fraction", "-mf"),
                               is_flag=True,
                               default=True,
                               help="Plot mutation_fraction plot. See Plots/gallery.")

opt_mutation_fraction_identity = Option(("--mutation_fraction_identity", "-mfi"),
                                        is_flag=True,
                                        default=True,
                                        help="Plot mutation_fraction_identity plot. See Plots/gallery.")

opt_base_coverage = Option(("--base_coverage", "-bc"),
                           is_flag=True,
                           default=True,
                           help="Plot base_coverage plot. See Plots/gallery.")

opt_mutations_in_barcodes = Option(("--mutations_in_barcodes", "-mib"),
                                   is_flag=True,
                                   default=False,
                                   help="Plot mutations_in_barcodes plot. See Plots/gallery.")

opt_mutations_per_read_per_sample = Option(("--mutation_per_read_per_reference", "-mprps"),
                                           is_flag=True,
                                           default=True,
                                           help="Plot mutation_per_read_per_reference plot. See Plots/gallery.")

# Logging options
opt_verbose = Option(("--verbose", "-v"),
                     count=True,
                     help="Print info or info+debug messages to stdout")
opt_quiet = Option(("--quiet",),
                   count=True,
                   help="Suppress warnings or warnings+errors to stdout")
opt_log = Option(("--log",),
                 type=Path(exists=False, dir_okay=False),
                 default=os.path.join(CWD, datetime.now().strftime(
                     "dreem_%Y-%m-%d_%H:%M:%S.log")))
opt_profile = Option(("--profile",),
                     type=Path(exists=False, dir_okay=False),
                     default="",
                     help="Profile code performance and log results to file")

# Misc
opt_version = Option(("--version",),
                     is_flag=True,
                     help="Show the version and exit.")


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
