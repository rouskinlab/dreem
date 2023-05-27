
from pathlib import Path

from ..core import path


def load(table: Path):
    """ Load an RNA profile from a table. """
    fields = path.parse(path.ModSeg, path.SampSeg, path.RefSeg, path.SectSeg,
                        path.MutTabSeg,
                        )
