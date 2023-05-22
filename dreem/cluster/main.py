from concurrent.futures import Future, ProcessPoolExecutor
from functools import partial
from logging import getLogger
from pathlib import Path

from click import command

from .krun import cluster_filt
from ..filt.load import FilterLoader
from ..util import docdef
from ..util.cli import (opt_filt,
                        opt_max_clusters, opt_em_runs,
                        opt_min_em_iter, opt_max_em_iter, opt_em_thresh,
                        opt_min_fmut_pos,
                        opt_parallel, opt_max_procs)
from ..util.parallel import get_num_parallel


logger = getLogger(__name__)


params = [
    # Input/output paths
    opt_filt,
    # Clustering options
    opt_max_clusters,
    opt_em_runs,
    opt_min_em_iter,
    opt_max_em_iter,
    opt_em_thresh,
    opt_min_fmut_pos,
    # Parallelization
    opt_max_procs,
    opt_parallel,
]


@command("cluster", params=params)
def cli(*args, **kwargs):
    return run(*args, **kwargs)


@docdef.auto()
def run(filt: tuple[str, ...], *,
        max_clusters: int,
        em_runs: int,
        min_em_iter: int,
        max_em_iter: int,
        em_thresh: float,
        min_fmut_pos: float,
        max_procs: int,
        parallel: bool) -> list[Path]:
    """ Run the clustering module. """
    if max_clusters == 0:
        # Exit immediately if the maximum number of clusters is 0.
        return list()
    loaders = list(map(FilterLoader.open, map(Path, filt)))
    # Determine how to parallelize clustering.
    n_tasks, n_procs_per_task = get_num_parallel(len(filt), max_procs,
                                                 parallel=parallel,
                                                 hybrid=False)
    # Run EM clustering on every section of every set of bit vectors.
    cluster_func = partial(cluster_filt,
                           max_clusters=max_clusters,
                           n_runs=em_runs,
                           min_iter=min_em_iter,
                           max_iter=max_em_iter,
                           conv_thresh=em_thresh,
                           max_procs=n_procs_per_task)
    if n_tasks > 1:
        # Run multiple sections of bit vectors in parallel.
        with ProcessPoolExecutor(max_workers=n_tasks) as pool:
            futures: list[Future] = list()
            for loader in loaders:
                futures.append(pool.submit(cluster_func, loader))
                logger.debug(f"Submitted EM clustering for {loader}")
            cluster_files: list[Path] = [future.result() for future in futures]
    else:
        # Run each section of bit vectors one at a time.
        cluster_files: list[Path] = list()
        for loader in loaders:
            cluster_files.append(cluster_func(loader))
    # Remove any None values (indicating that clustering failed).
    return list(filter(None, cluster_files))
