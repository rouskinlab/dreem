from logging import getLogger
from pathlib import Path

from click import command

from .krun import cluster
from ..core import docdef, path
from ..core.cli import (opt_call, opt_max_clusters, opt_em_runs,
                        opt_min_em_iter, opt_max_em_iter, opt_em_thresh,
                        opt_parallel, opt_max_procs)
from ..core.parallel import dispatch

logger = getLogger(__name__)

params = [
    # Input files
    opt_call,
    # Clustering options
    opt_max_clusters,
    opt_em_runs,
    opt_min_em_iter,
    opt_max_em_iter,
    opt_em_thresh,
    # Parallelization
    opt_max_procs,
    opt_parallel,
]


@command(path.MOD_CLUST, params=params)
def cli(*args, **kwargs):
    return run(*args, **kwargs)


@docdef.auto()
def run(call: tuple[str, ...], *,
        max_clusters: int,
        em_runs: int,
        min_em_iter: int,
        max_em_iter: int,
        em_thresh: float,
        max_procs: int,
        parallel: bool) -> list[Path]:
    """ Run the clustering module. """
    if max_clusters == 0:
        # Exit immediately if the maximum number of clusters is 0.
        return list()
    # Run clustering on each set of called mutations.
    return dispatch(cluster, max_procs, parallel,
                    args=[(Path(file),) for file in call],
                    kwargs=dict(max_clusters=max_clusters,
                                n_runs=em_runs,
                                min_iter=min_em_iter,
                                max_iter=max_em_iter,
                                conv_thresh=em_thresh))
