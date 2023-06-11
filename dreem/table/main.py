from itertools import chain
from logging import getLogger
from pathlib import Path

from click import command

from .write import write
from ..core import docdef, path
from ..core.cli import (opt_rel, opt_mask, opt_clust,
                        opt_max_procs, opt_parallel, opt_rerun)
from ..core.parallel import as_list_of_tuples, dispatch

logger = getLogger(__name__)

params = [opt_rel, opt_mask, opt_clust,
          opt_max_procs, opt_parallel, opt_rerun]


@command(path.MOD_TABLE, params=params)
def cli(*args, **kwargs):
    """ Tabulate per-base and per-read mutation counts from 'relate' and
    'mask', and mutation rates and read memberships from 'cluster'. """
    return run(*args, **kwargs)


@docdef.auto()
def run(rel: tuple[str, ...], mask: tuple[str, ...], clust: tuple[str, ...],
        max_procs: int, parallel: bool, **kwargs):
    """
    Run the table module.
    """
    rels = path.find_files_multi(map(Path, rel), [path.RelateRepSeg])
    if rels:
        logger.debug(f"Found relate report files: {rels}")
    masks = path.find_files_multi(map(Path, mask), [path.MaskRepSeg])
    if masks:
        logger.debug(f"Found mask report files: {masks}")
    clusts = path.find_files_multi(map(Path, clust), [path.ClustRepSeg])
    if clusts:
        logger.debug(f"Found cluster report files: {clusts}")
    tasks = as_list_of_tuples(chain(rels, masks, clusts))
    return list(chain(*dispatch(write, max_procs, parallel,
                                args=tasks, kwargs=kwargs,
                                pass_n_procs=False)))
