import numpy as np

import dreem
import click
import os, sys
import pandas as pd
import json
from dreem import util

def cluster_likelihood(bit_vector_path, fasta, section_start, section_end, temp_folder, N_clusters, max_clusters, signal_thresh, info_thresh, include_G_U, include_del, min_reads, convergence_cutoff, num_runs):
    """
    Run the clustering algorithm for a given section of a construct.

    Parameters
    ----------
    bit_vector_path: str
        Path to the bit vector.
    fasta: str
        Path to the fasta file.
    section_start: int
        Start of the section.
    section_end: int
        End of the section.
    temp_folder: str
        Path to the temporary folder.
    N_clusters: int
        Number of clusters
    max_clusters: int
        Maximum number of clusters
    signal_thresh: float
        Signal threshold
    info_thresh: float
        Information threshold
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
    # Run the clustering
    try:
        cluster_likelihoods = __cluster_likelihood(bit_vector_path, fasta, section_start, section_end, temp_folder, N_clusters, max_clusters, signal_thresh, info_thresh, include_G_U, include_del, min_reads, convergence_cutoff, num_runs)
        return cluster_likelihoods
    except:
        return None

def __cluster_likelihood(bit_vector_path, fasta, section_start, section_end, temp_folder, N_clusters, max_clusters, signal_thresh, info_thresh, include_G_U, include_del, min_reads, convergence_cutoff, num_runs):

    placeholder = {
        'cluster_1': np.random.random(),
        'cluster_2': np.random.random(),
        'cluster_3': np.random.random(),
    }

    return placeholder

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
    input_dir = args['input_dir']
    fasta = args['fasta']
    root = args['out_dir']
    sample = args['sample']
    library = pd.read_csv(args['library'])
    temp_folder = os.path.join(root,'temp','clustering') 
    output_folder = os.path.join(root,'output','clustering') 
    output_file = os.path.join(output_folder, sample+'.json')
    N_clusters = args['n_clusters']
    max_clusters = args['max_clusters']
    signal_thresh = args['signal_thresh']
    info_thresh = args['info_thresh']
    include_G_U = args['include_g_u']
    include_del = args['include_del']
    min_reads = args['min_reads']
    convergence_cutoff = args['convergence_cutoff']
    num_runs = args['num_runs']
    samples

    # Create the output folder
    util.make_folder(output_folder)
    util.make_folder(temp_folder)

    cluster_likelihoods = {}
    
    # Remove this
    raise NotImplementedError('This module is not implemented yet')

    # Run the clustering
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

