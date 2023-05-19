from itertools import chain
from logging import getLogger
from pathlib import Path
from typing import Callable, Iterable

import pandas as pd

from .emalgo import EmClustering
from ..bit.vect import BitVector
from ..util import path
from ..util.epu import get_common_attrib
from ..util.report import DECIMAL_PRECISION

logger = getLogger(__name__)

IDX_NCLUSTERS = "NumClusters"
IDX_CLUSTER = "Cluster"
IDXS_CLUSTERS = IDX_NCLUSTERS, IDX_CLUSTER


def kc_pairs(ks: Iterable[int]):
    """
    Return the multi-index for a data frame with each number of clusters
    in ```ks```. The first level of the multi-index is the total number
    of clusters, ```k```; the second level is the number of the cluster,
    ```c```. For example, ```(3, 2)``` signifies cluster number 2 from
    the best run of EM clustering with ```k``` = 3 clusters in total.
    Note that the length of the multi-index equals ```sum(ks)```.

    Examples
    --------
    >>> kc_pairs([1, 2, 3]).tolist()
    [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3)]
    """
    return pd.MultiIndex.from_tuples([(k, c) for k in ks
                                      for c in range(1, k + 1)],
                                     names=IDXS_CLUSTERS)


def write_results(clusters: dict[int, list[EmClustering]],
                  out_dir: Path):
    """ Write CSV files of the proportions, mutation rates, and read
    responsibilities for each cluster. Return the file paths. """
    bvec = get_common_attrib(chain(*clusters.values()), "bvec", key=id)
    # Proportions: proportion of each cluster in the ensemble
    props = write_tables(build_tables(clusters, EmClustering.output_props),
                         bvec, out_dir, path.CLUST_PROP_TABLE)
    # Mutation rates: fraction of mutated bits at each position
    mus = write_tables(build_tables(clusters, EmClustering.output_mus),
                       bvec, out_dir, path.CLUST_MUS_TABLE)
    # Responsibilities: likelihood that each read came from each cluster
    resps = write_tables(build_tables(clusters, EmClustering.output_resps),
                         bvec, out_dir, path.CLUST_RESP_TABLE, gzip=True)
    return props, mus, resps


def build_tables(clusters: dict[int, list[EmClustering]],
                 output_func: Callable[[EmClustering], pd.DataFrame]):
    """ Build a DataFrame of one attribute of one or more clusters. """
    tables: dict[int, pd.DataFrame] = dict()
    for k, runs in clusters.items():
        for r, run in enumerate(runs):
            ktable = output_func(run)
            if r not in tables:
                tables[r] = pd.DataFrame(index=ktable.index,
                                         columns=kc_pairs(clusters),
                                         dtype=float)
            for c, column in ktable.items():
                tables[r].loc[:, (k, c)] = column
    return tables


def write_tables(dfs: dict[int, pd.DataFrame], bvec: BitVector, out_dir: Path,
                 name: str, gzip: bool = False) -> dict[int, Path]:
    return {run: write_table(df, bvec, out_dir, name, run, gzip)
            for run, df in dfs.items()}


def table_path(out_dir: Path, sample: str, ref: str, sect: str,
               table: str, run: int, gzip: bool):
    return path.buildpar(path.ModSeg, path.SampSeg, path.RefSeg, path.SectSeg,
                         path.ClustTabSeg, top=out_dir, module=path.MOD_CLUST,
                         sample=sample, ref=ref, sect=sect, table=table,
                         run=run, ext=(path.CSVZIP_EXT if gzip
                                       else path.CSV_EXT))


def write_table(df: pd.DataFrame, bvec: BitVector, out_dir: Path,
                table: str, run: int, gzip: bool = False):
    """ Write a DataFrame of one clustering attribute to a CSV file. """
    file = table_path(out_dir, bvec.loader.sample, bvec.loader.ref,
                      bvec.section.name, table, run, gzip)
    df.round(DECIMAL_PRECISION).to_csv(file, header=True, index=True)
    logger.info(f"Wrote {table} run {run} of {bvec} to {file}")
    return file


def load_cluster_resps(cluster_resps_file: Path):
    """ Load the responsibilities of the reads for each cluster. """
    return pd.read_csv(cluster_resps_file,
                       index_col=[0],
                       header=list(range(len(IDXS_CLUSTERS))))
