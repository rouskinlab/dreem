from logging import getLogger
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from dreem.vector.bits import BitVector
from .emalgo import EmClustering
from ..util import path

logger = getLogger(__name__)

IDX_NCLUSTERS = "NumClusters"
IDX_CLUSTER = "Cluster"
IDXS_CLUSTERS = IDX_NCLUSTERS, IDX_CLUSTER


def write_results(bvec: BitVector,
                  clusters: dict[int, EmClustering],
                  out_dir: Path):
    """ Write CSV files of the proportions, mutation rates, and read
    responsibilities for each cluster. Return the responsibilities. """
    # Define the cluster numbers.
    nc_pairs = pd.MultiIndex.from_tuples([(nclust, clust) for nclust in clusters
                                          for clust in range(1, nclust + 1)],
                                         names=IDXS_CLUSTERS)
    # Proportions: proportion of each cluster in the ensemble
    write_table(build_table(clusters, EmClustering.output_props, nc_pairs),
                bvec, out_dir, path.CLUST_PROP_TABLE)
    # Mutation rates: fraction of mutated bits at each position
    write_table(build_table(clusters, EmClustering.output_mus, nc_pairs),
                bvec, out_dir, path.CLUST_MUT_TABLE)
    # Responsibilities: likelihood that each read came from each cluster
    resps_file = write_table(build_table(clusters, EmClustering.output_resps,
                                         nc_pairs),
                             bvec, out_dir, path.CLUST_RESP_TABLE, gzip=True)
    return resps_file


def build_table(clusters: dict[int, EmClustering],
                output_func: Callable[[EmClustering], pd.DataFrame],
                columns: pd.MultiIndex):
    """ Build a DataFrame of one attribute of one or more clusters. """
    table = None
    for k, cluster in clusters.items():
        for i, column in output_func(cluster).items():
            if table is None:
                table = pd.DataFrame(index=column.index, columns=columns,
                                     dtype=float)
            table.loc[:, (k, i)] = column
    return table


def write_table(table: pd.DataFrame, bvec: BitVector,
                out_dir: Path, name: str, gzip: bool = False):
    """ Write a DataFrame of one clustering attribute to a CSV file. """
    file = path.build(path.ModSeg, path.SampSeg, path.RefSeg,
                      path.SectSeg, path.ClustTabSeg,
                      top=out_dir, module=path.MOD_CLUST,
                      sample=bvec.loader.sample, ref=bvec.loader.ref,
                      end5=bvec.section.end5, end3=bvec.section.end3,
                      table=name, ext=(path.CSVZIP_EXT if gzip
                                       else path.CSV_EXT))
    table.to_csv(file, header=True, index=True)
    logger.info(f"Wrote {name} of {bvec} to {file}")
    return file


def load_cluster_resps(cluster_resps_file: Path):
    """ Load the responsibilities of the reads for each cluster. """
    return pd.read_csv(cluster_resps_file,
                       index_col=[0],
                       header=list(range(len(IDXS_CLUSTERS))))
