from itertools import chain
from logging import getLogger
from pathlib import Path

from click import command

from .io import write
from ..core import docdef, path
from ..core.cli import (opt_rel, opt_call, opt_clust,
                        opt_max_procs, opt_parallel, opt_rerun)
from ..core.parallel import dispatch

logger = getLogger(__name__)

params = [opt_rel, opt_call, opt_clust,
          opt_max_procs, opt_parallel, opt_rerun]


@command(path.MOD_TABLE, params=params)
def cli(**kwargs):
    return run(**kwargs)


@docdef.auto()
def run(rel: tuple[str, ...],
        call: tuple[str, ...],
        clust: tuple[str, ...],
        max_procs: int,
        parallel: bool,
        **kwargs):
    """
    Run the table module.
    """
    tasks = [(Path(report_file),) for report_file in chain(rel, call, clust)]
    return list(chain(*dispatch(write, max_procs, parallel,
                                args=tasks, kwargs=kwargs,
                                pass_n_procs=False)))
