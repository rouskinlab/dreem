from logging import getLogger
from pathlib import Path
from typing import Callable

import pandas as pd

from .emalgo import EmClustering
from ..mask.load import MaskLoader
from ..core import path

logger = getLogger(__name__)


FLOAT_PRECISION = 6  # number of digits behind the decimal point


def write_results(loader: MaskLoader, k_runs: dict[int, list[EmClustering]]):
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


def write_table(loader: MaskLoader,
                rank: int,
                run: EmClustering,
                output_func: Callable[[EmClustering], pd.DataFrame],
                table: str, *,
                gzip: bool = False):
    """ Write a DataFrame of one clustering attribute to a CSV file. """
    data = output_func(run)
    file = table_path(loader.out_dir, loader.sample, loader.ref, loader.sect,
                      table, run.ncls, rank, gzip)
    data.round(FLOAT_PRECISION).to_csv(file, header=True, index=True)
    logger.info(f"Wrote {table} of {run} to {file}")
    return file
