
from click_option_group import optgroup
import click
from dreem import util, pipeline


@click.command()
@optgroup.group('Files and folders paths')
@optgroup.option('--fasta', '-fa', type=click.Path(exists=True), help='Path to the fasta file', required=True)
@optgroup.option('--fastq1', '-fq1', help='Paths to the fastq1 file (forward primer). Enter multiple times for multiple files', multiple=True, type=click.Path(exists=True), required=True)
@optgroup.option('--fastq2', '-fq2', help='Paths to the fastq2 file (reverse primer). Enter multiple times for multiple files', multiple=True, type=click.Path(exists=True))
@optgroup.option('--samples', '-s', type=click.Path(exists=True), help='Path to the samples.csv file')
@optgroup.option('--library', '-l', type=click.Path(exists=True), help='Path to the library.csv file')
@optgroup.option('--out_dir', '-o', default=os.getcwd(), type=click.Path(exists=True), help='Main directory to run the pipeline in')
@optgroup.option('--clear_directories', '-cd', type=bool, default=False, help='Remove temp and output folders before running')

@optgroup.group('Construct selection')
@optgroup.option("-c", "--coords", type=(str, int, int), multiple=True, help="coordinates for reference: '-c ref-name first last'")
@optgroup.option("-p", "--primers", type=(str, util.Primer, util.Primer), multiple=True, help="primers for reference: '-p ref-name fwd rev'")
@optgroup.option("--fill/--no-fill", default=False,
              help="Fill in coordinates of reference sequences for which "
                   "neither coordinates nor primers were given (default: no).")


@optgroup.group('Demultiplexing parameters')
@optgroup.option('--demultiplexing', '-dx', type=bool, help='Use demultiplexing', default=False)
@optgroup.option('--barcode_start', '-bs', type=int, help='Start position of the barcode in the read')
@optgroup.option('--barcode_end', '-be', type=int, help='End position of the barcode in the read')

#@optgroup.group('Bowtie parameters') #TODO

#@optgroup.group('Cutadapt parameters') #TODO 

@optgroup.group('Vectoring parameters')
@optgroup.option("-P", "--parallel",
              type=click.Choice(["profiles", "reads", "off", "auto"],
                                case_sensitive=False),
              default="auto",
              help="Parallelize the processing of mutational PROFILES or "
              "READS within each profile, turn parallelization OFF, or AUTO"
              "matically choose the parallelization method (default: auto).")

@optgroup.group('Clustering parameters')
@optgroup.option('--clusters', '-cl', type=bool, help='Use clustering', default=False)
@optgroup.option('--n_clusters', '-nc', type=int, help='Number of clusters', default=None)
@optgroup.option('--max_clusters', '-mc', type=int, help='Maximum number of clusters', default=None)
@optgroup.option('--signal_thresh', '-st', type=float, help='Signal threshold', default=None)
@optgroup.option('--info_thresh', '-it', type=float, help='Information threshold', default=None)
@optgroup.option('--include_g_u', '-igu', type=bool, help='Include G and U', default=None)
@optgroup.option('--include_del', '-id', type=bool, help='Include deletions', default=None)
@optgroup.option('--min_reads', '-mr', type=int, help='Minimum number of reads', default=None)
@optgroup.option('--convergence_cutoff', '-cc', type=float, help='Convergence cutoff', default=None)
@optgroup.option('--num_runs', '-nr', type=int, help='Number of runs', default=None)

@optgroup.group('Aggregate parameters')
@optgroup.option('--rnastructure_path', '-rs', type=click.Path(exists=True), help='Path to RNAstructure, to predict structure and free energy', default=None)
@optgroup.option('--rnastructure_temperature', '-rst', type=bool, help='Use sample.csv temperature values for RNAstructure', default=False)
@optgroup.option('--rnastructure_fold_args', '-rsa', type=str, help='Arguments to pass to RNAstructure fold', default=None)
@optgroup.option('--rnastructure_dms', '-rsd', type=bool, help='Use the DMS signal to amke predictions with RNAstructure', default=False)
@optgroup.option('--rnastructure_dms_min_unpaired_value', '-rsdmin', type=int, help='Minimum unpaired value for using the dms signal as an input for RNAstructure', default=0.01)
@optgroup.option('--rnastructure_dms_max_paired_value', '-rsdmax', type=int, help='Maximum paired value for using the dms signal as an input for RNAstructure', default=0.05)
@optgroup.option('--rnastructure_partition', '-rspa', type=bool, help='Use RNAstructure partition function to predict free energy', default=False)
@optgroup.option('--rnastructure_probability', '-rspr', type=bool, help='Use RNAstructure partition function to predict per-base mutation probability', default=False)
@optgroup.option('--poisson', '-po', type=bool, help='Predict Poisson confidence intervals', default=True)

@optgroup.group('Misc. parameters')
@optgroup.option('--verbose', '-v', type=bool, help='Verbose output', default=False)


def cli(**args):
    pipeline.run(**args)

if __name__ == "__main__":
    cli()
