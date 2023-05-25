from logging import getLogger
from pathlib import Path

from click import command

from .io import (RelVecPosTableWriter, RelVecReadTableWriter,
                 BitVecPosTableWriter, BitVecReadTableWriter,
                 write)
from ..core import docdef, path
from ..core.cli import opt_rel, opt_call, opt_clust, opt_max_procs, opt_parallel
from ..core.parallel import dispatch

logger = getLogger(__name__)

params = [opt_rel, opt_call, opt_max_procs, opt_parallel]


@command(path.MOD_TABLE, params=params)
def cli(**kwargs):
    return run(**kwargs)


@docdef.auto()
def run(rel: tuple[str, ...],
        call: tuple[str, ...],
        clust: tuple[str, ...],
        max_procs: int,
        parallel: bool):
    """
    Run the table module.
    """
    tasks = list()
    for report in rel:
        for writer in (RelVecPosTableWriter, RelVecReadTableWriter):
            tasks.append((Path(report), writer))
    for report in call:
        for writer in (BitVecPosTableWriter, BitVecReadTableWriter):
            tasks.append((Path(report), writer))
    for report in clust:
        pass
    return dispatch(write, max_procs, parallel, args=tasks, pass_n_procs=False)
