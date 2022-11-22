import yaml, sys, os

from click_option_group import optgroup
import click
from dreem import util


@click.command()
@optgroup.group('Files and folders paths')
@optgroup.option('--fasta', '-fa', type=click.Path(exists=True), help='Path to the fasta file', required=True)
@optgroup.option('--fastq1', '-fq1', help='Paths to the fastq1 file (forward primer). Enter multiple times for multiple files', multiple=True, type=click.Path(exists=True), required=True)
@optgroup.option('--fastq2', '-fq2', help='Paths to the fastq2 file (reverse primer). Enter multiple times for multiple files', multiple=True, type=click.Path(exists=True))
@optgroup.option('--samples', '-s', type=click.Path(exists=True), help='Path to the samples.csv file')
@optgroup.option('--library', '-l', type=click.Path(exists=True), help='Path to the library.csv file')
@optgroup.option('--output', '-o', default=os.getcwd(), type=click.Path(exists=True), help='Main directory to run the pipeline in')
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
@optgroup.option('--clustering', '-cl', type=bool, help='Use clustering', default=False)
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


def run(**args):
    """Run DREEM.

     Args:
         args (_type_): _description_
    """

    # make output and temp folders
    for folder in ['output', 'temp']:
        if args['clear_directories']:
            util.clear_folder(os.path.join(args['output'], folder))
        util.make_folder(os.path.join(args['output'], folder))
    
    # Run DREEM
    print("""

    ========================================

                RUNNING   DREEM

    ========================================

    """)
    if args['demultiplexing']:
        print('\ndemultiplexing \n------------------')
        cmd = util.make_cmd({k:v for k,v in args.items() if k in ['output','fastq1','fastq2','library','barcode_start','barcode_stop']}, module='demultiplexing')
        print(util.run_cmd(cmd))
        print('demultiplexing done')
    
    ## Alignment
    print('\nalignment \n----------------')
    fastq1, fastq2 = [], []
    if args['demultiplexing']:
        for sample in os.listdir(os.path.join(args['output'], 'output', 'demultiplexing')):
            path = os.path.join(os.path.join(args['output'], 'output', 'demultiplexing'), sample)
            fastq1 = fastq1 + [os.path.join(path, f) for f in os.listdir(path) if f.endswith('_R1.fastq')]
            fastq2 = fastq2 + [os.path.join(path, f) for f in os.listdir(path) if f.endswith('_R2.fastq')]
    else:
        fastq1, fastq2 = args['fastq1'], args['fastq2']

    for f1 in fastq1:
        for f2 in fastq2:
            if f1[:-len('_R1.fastq')] == f2[:-len('_R2.fastq')]:
                print('Aligning this fastq pair: ', '\n   ',f1, '\n   ',f2)
                sample = f1.split('/')[-2]
                print(util.run_cmd('dreem-alignment -fa {} -fq1 {} -fq2 {} -o {} -sd {}'.format(args['fasta'], f1, f2, args['output'], sample)))

    ## Vectoring
    print('\nvectoring \n------------------')
    path_alignment = os.path.join(args['output'], 'output', 'alignment')
    cmd = util.make_cmd({k:v for k,v in args.items() if k in ['output','fasta','coords','primers','parallel']}, module='vectoring') \
            + (' --fill ' if args['fill'] else ' --no-fill ')\
            + ' --bam_dir ' + ' --bam_dir '.join([os.path.join(path_alignment, s) for s in os.listdir(path_alignment)]) 
    print(util.run_cmd(cmd))

    ## Clustering
    if args['clustering']:
        print('\nclustering \n------------------')
        path_vectoring = os.path.join(args['output'], 'output', 'vectoring')
        for sample in os.listdir(path_vectoring):
            cmd = util.make_cmd({k:v for k,v in args.items() if k in ['fasta','library','output','N_clusters','max_clusters','signal_thresh','info_thresh','include_G_U','include_del','min_reads','convergence_cutoff','num_runs']}, module='clustering') \
                    + ' --bv_dir ' + os.path.join(path_vectoring, sample) \
                    + ' --name ' + sample
            print(util.run_cmd(cmd))

    ## Aggregate
    print('\naggregating \n------------------')
    path_vectoring = os.path.join(args['output'], 'output', 'vectoring')
    for sample in os.listdir(path_vectoring):
        path_clustering = os.path.join(args['output'], 'output', 'clustering', sample+'.json') if args['clustering'] else None
        bv_dir = os.path.join(path_vectoring, sample)
        cmd = util.make_cmd({k:v for k,v in args.items() if k in ['output','samples','library','rnastructure_path','rnastructure_temperature','rnastructure_fold_args','rnastructure_dms','rnastructure_dms_min_unpaired_value','rnastructure_dms_max_paired_value','rnastructure_partition','rnastructure_probability','poisson','verbose']}, module='aggregate') \
                + ' --bv_dir ' + bv_dir \
                + ' --sample ' + sample \
                + (' --clustering ' + path_clustering if args['clustering'] is not None else '')
        print(util.run_cmd(cmd))

    print('Done!')