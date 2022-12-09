import numpy as np

import dreem
import click
import os, sys
import pandas as pd
import json
from dreem import util


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
        
    Returns
    -------
    
    clustering: dict
        The structure is the following:
        {construct_1}: # name of the bitvector file (from fasta file)
            {section_1}: # name of the clustered section (from library)
                {read_1}:  # read number (from bitvector)
                    K1_1: likelihood that this read belongs to cluster 1 when using 1 cluster
                    K2_1: likelihood that this read belongs to cluster 1 when using 2 clusters
                    K2_2: likelihood that this read belongs to cluster 2 when using 2 clusters
                    K3_1: likelihood that this read belongs to cluster 1 when using 3 clusters
                    K3_2: likelihood that this read belongs to cluster 2 when using 3 clusters
                    K3_3: likelihood that this read belongs to cluster 3 when using 3 clusters   
    
    """

    # Create the output folder
    os.makedirs(args['out_dir'], exist_ok=True)
    os.makedirs(os.path.join(args['out_dir'], 'temp'), exist_ok=True)
    
    cluster_likelihoods = {}
    
    # Remove this
    raise NotImplementedError('This module is not implemented yet')

    # Run the clustering 
    cluster_likelihoods = {}
    
    # Save the cluster likelihoods
    with open(output_file, 'w') as f:
        json.dump(cluster_likelihoods, f)

    return 1

