import dreem
import click
import os, sys
import pandas as pd
import json
from dreem.clustering.clustering import cluster_likelihood
from dreem import util


@click.command()
@click.option('--fasta', '-fa', type=click.Path(exists=True), help='Path to the fasta file', required=True)
@click.option('--input_dir', '-id', help='Path to the bit vector folder or list of paths to the bit vector folders.', type=click.Path(exists=True), multiple=True)
@click.option('--out_dir', '-o', default=os.path.join(os.getcwd()), type=click.Path(exists=True), multiple = True, help='Where to output files')
@click.option('--library', '-l', type=click.Path(exists=True), help='Path to the library.csv file')
@click.option('--n_clusters', '-nc', type=int, help='Number of clusters', default=None)
@click.option('--max_clusters', '-mc', type=int, help='Maximum number of clusters', default=None)
@click.option('--signal_thresh', '-st', type=float, help='Signal threshold', default=None)
@click.option('--info_thresh', '-it', type=float, help='Information threshold', default=None)
@click.option('--include_g_u', '-igu', type=bool, help='Include G and U', default=None)
@click.option('--include_del', '-idel', type=bool, help='Include deletions', default=None)
@click.option('--min_reads', '-mr', type=int, help='Minimum number of reads', default=None)
@click.option('--convergence_cutoff', '-cc', type=float, help='Convergence cutoff', default=None)
@click.option('--num_runs', '-nr', type=int, help='Number of runs', default=None)


def run(**args):
    """Run the clustering module.

    Clusters the reads of all given bitvectors and outputs the likelihoods of the clusters as `name`.json in the directory `output_path`, using `temp_path` as a temp directory.
    Each bitvector is a file containing the reads of a construct. Bitvectors from the same sample should be grouped in a folder and the path to the folder should be given as `bv_dir`.
    `name` is the name of the output file, and should be the sample name.

    Parameters from args:
    ---------------------

    fasta: str
        Path to the reference FASTA file.
    input_dir: str
        Path to the bit vector folder or list of paths to the bit vector folders.
    out_dir: str
        Path to the output folder.
    library: str
        Path to the library.csv file.
    N_clusters: int
        Number of clusters
    max_clusters: int
        Maximum number of clusters
    signal_thresh: float
        Signal threshold
    info_thresh: float
        Float from 0 to 1, where 1 means that all bases are unvalid and 0 means that all bases are valid (valid:= just 0s and 1s in the bitvector). If info_thresh of a read is above this threshold, it is considered unvalid and isn't used.
    include_G_U: bool
        Include G and U
    include_del: bool
        Include deletions
    min_reads: int
        Minimum number of reads
    convergence_cutoff: float
        Convergence cutoff
    num_runs: int
        Number of runs
    """

    # Extract the arguments
    samples = [os.path.basename(os.path.normpath(sample)) for sample in args['input_dir']] if isinstance(args['input_dir'], list) else [os.path.basename(os.path.normpath(args['input_dir']))]
    constructs, bitvectors_path = {}, {}
    for sample in samples:
        constructs[sample] = [construct[:-len('.orc')] for construct in os.listdir(os.path.join(args['input_dir'], sample)) if construct.endswith('.orc')]
        bitvectors_path[sample] = [os.path.join(args['input_dir'], sample, construct + '.orc') for construct in constructs[sample]]
    fasta = args['fasta']
    root = args['out_dir']
    library = pd.read_csv(args['library'])
    temp_folder = os.path.join(root,'temp','clustering') 
    output_folder = os.path.join(root,'output','clustering') 
    output_file = os.path.join(output_folder, 'clustering.json')
    N_clusters = args['n_clusters']
    max_clusters = args['max_clusters']
    signal_thresh = args['signal_thresh']
    info_thresh = args['info_thresh']
    include_G_U = args['include_g_u']
    include_del = args['include_del']
    min_reads = args['min_reads']
    convergence_cutoff = args['convergence_cutoff']
    num_runs = args['num_runs']

    # Create the output folder
    util.make_folder(output_folder)
    util.make_folder(temp_folder)

    cluster_likelihoods = {}

    # Run the clustering
    for sample in samples:
        for construct, bv_path in zip(constructs[sample], bitvectors_path[sample]):
            cluster_likelihoods[construct] = {}
            assert construct in library['construct'].values, 'Construct {} not in library'.format(construct)
            library_for_this_construct = library[library['construct'] == construct]
            for _, row in library_for_this_construct.iterrows():
                section, section_start, section_end = row['section'], row['section_start'], row['section_end']
                cluster_likelihoods[construct][section] = cluster_likelihood(bv_path, fasta, section_start, section_end, temp_folder, N_clusters, max_clusters, signal_thresh, info_thresh, include_G_U, include_del, min_reads, convergence_cutoff, num_runs)
                assert cluster_likelihoods[construct][section] is not None, 'Clustering failed for construct {} and section {}'.format(construct, section)
    # Save the cluster likelihoods
    with open(output_file, 'w') as f:
        json.dump(cluster_likelihoods, f)

    return 1

if __name__ == '__main__':
    run()
