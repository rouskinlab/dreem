import dreem
import click
import os, sys
import pandas as pd
import json
from dreem.clustering import cluster_likelihood
from main import run


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


def cli(**args):
    run(**args)

if __name__ == '__main__':
    run()
