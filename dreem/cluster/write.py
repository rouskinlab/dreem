from logging import getLogger
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from .bvec import BitVector
from .emalgo import EmClustering
from ..util import path

logger = getLogger(__name__)


def write_results(bitvectors: BitVector,
                  clusters: dict[int, EmClustering],
                  out_dir: Path):
    """ Write CSV files of the proportions, mutation rates, and read
    responsibilities for each cluster. Return the responsibilities. """
    # Define the cluster numbers.
    ki_pairs = pd.MultiIndex.from_tuples([(k, i) for k in clusters
                                          for i in range(1, k + 1)],
                                         names=["k", "i"])
    # Proportions: proportion of each cluster in the ensemble
    write_table(build_table(clusters, EmClustering.output_props,
                            EmClustering.props_index(), ki_pairs),
                bitvectors, out_dir, path.CLUST_PROP_TABLE)
    # Mutation rates: fraction of mutated bits at each position
    write_table(build_table(clusters, EmClustering.output_mus,
                            bitvectors.positions, ki_pairs),
                bitvectors, out_dir, path.CLUST_MUT_TABLE)
    # Responsibilities: likelihood that each read came from each cluster
    return write_table(build_table(clusters, EmClustering.output_resps,
                                   bitvectors.read_names, ki_pairs),
                       bitvectors, out_dir, path.CLUST_RESP_TABLE, gzip=True)


def build_table(clusters: dict[int, EmClustering],
                output_func: Callable[[EmClustering], pd.DataFrame],
                index: list | np.ndarray | pd.Index,
                columns: pd.MultiIndex):
    """ Build a DataFrame of one attribute of one or more clusters. """
    table = pd.DataFrame(dtype=float, index=index, columns=columns)
    for k, cluster in clusters.items():
        for i, column in output_func(cluster).items():
            table.loc[:, (k, i)] = column
    return table


def write_table(table: pd.DataFrame, bvec: BitVector,
                out_dir: Path, name: str, gzip: bool = False):
    """ Write a DataFrame of one clustering attribute to a CSV file. """
    file = path.build(path.ModSeg, path.SampSeg, path.RefSeg,
                      path.SectSeg, path.ClustTabSeg,
                      top=out_dir, module=path.MOD_CLUST, sample=bvec.sample,
                      ref=bvec.ref, end5=bvec.end5, end3=bvec.end3, table=name,
                      ext=(path.CSVZIP_EXT if gzip else path.CSV_EXT))
    table.to_csv(file, header=True, index=True)
    logger.info(f"Wrote {name} of {bvec} to {file}")
    return file
