from logging import getLogger
from pathlib import Path
from typing import Callable, Iterable

import pandas as pd

from .emalgo import EmClustering
from ..call.load import BitVecLoader
from ..core import path
from ..core.report import DECIMAL_PRECISION

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


def write_results(loader: BitVecLoader, k_runs: dict[int, list[EmClustering]]):
    """ Write CSV files of the proportions, mutation rates, counts, and
    read responsibilities, for each run. Return the file paths. """
    for k, runs in k_runs.items():
        for rank, run in enumerate(runs):
            if run.ncls != k:
                logger.error(f"{run} does not have {k} clusters")
                continue
            # Proportions: proportion of each cluster in the ensemble
            write_table(loader, rank, run,
                        EmClustering.output_props, path.CLUST_PROP_RUN_TABLE)
            # Mutation rates: fraction of mutated bits at each position
            write_table(loader, rank, run,
                        EmClustering.output_mus, path.CLUST_MUS_RUN_TAB)
            # Responsibilities: likelihood that a read came from a cluster
            write_table(loader, rank, run,
                        EmClustering.output_resps, path.CLUST_RESP_RUN_TABLE,
                        gzip=True)
            # Counts: observed and expected counts of each bit vector
            write_table(loader, rank, run,
                        EmClustering.output_counts, path.CLUST_COUNT_RUN_TABLE,
                        gzip=True)


def table_path(out_dir: Path, sample: str, ref: str, sect: str,
               table: str, k: int, run: int, gzip: bool):
    return path.buildpar(path.ModSeg, path.SampSeg, path.RefSeg, path.SectSeg,
                         path.ClustTabSeg, top=out_dir, module=path.MOD_CLUST,
                         sample=sample, ref=ref, sect=sect, table=table, k=k,
                         run=run, ext=(path.CSVZIP_EXT if gzip
                                       else path.CSV_EXT))


def write_table(loader: BitVecLoader,
                rank: int,
                run: EmClustering,
                output_func: Callable[[EmClustering], pd.DataFrame],
                table: str, *,
                gzip: bool = False):
    """ Write a DataFrame of one clustering attribute to a CSV file. """
    data = output_func(run)
    file = table_path(loader.out_dir, loader.sample, loader.ref, loader.sect,
                      table, run.ncls, rank, gzip)
    data.round(DECIMAL_PRECISION).to_csv(file, header=True, index=True)
    logger.info(f"Wrote {table} of {run} to {file}")
    return file


def load_cluster_resps(cluster_resps_file: Path):
    """ Load the responsibilities of the reads for each cluster. """
    return pd.read_csv(cluster_resps_file,
                       index_col=[0],
                       header=list(range(len(IDXS_CLUSTERS))))
