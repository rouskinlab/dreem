from abc import ABC, abstractmethod
from pathlib import Path

import pandas as pd

from ..call.load import BitVecLoader
from ..cluster.load import ClusterLoader
from ..core import path


class Table(ABC):
    @classmethod
    @abstractmethod
    def seg(cls) -> path.Segment | None:
        return

    @classmethod
    def path_segs(cls):
        return path.ModSeg, path.SampSeg, path.RefSeg, path.SectSeg, cls.seg()

    @classmethod
    @abstractmethod
    def fields(cls) -> tuple[str, ...]:
        return tuple()




class TableWriter(Table, ABC):
    pass


class CallTable(Table, ABC):
    @classmethod
    def seg(cls):
        return path.MutTabSeg







def call(report_file: Path):
    """ Calculate tabular data for sets of bit vectors. """
    # Load the bit vectors.
    bitvec = BitVecLoader.open(report_file)
    mus = bitvec.mus


def cluster(report_file: Path):
    """ Calculate tabular data for clusters of bit vectors. """
    clusts = ClusterLoader.open(report_file)
