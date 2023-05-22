from __future__ import annotations
from logging import getLogger

from ..util import path
from ..util.report import (Report, calc_seqlen, calc_time_taken,
                           calc_perc_vec, calc_n_batches, calc_speed)

logger = getLogger(__name__)


class MutVecReport(Report):
    __slots__ = ("sample", "ref", "seq", "length",
                 "n_vectors", "n_readerr", "perc_vec", "checksums",
                 "n_batches", "began", "ended", "taken", "speed")

    def __init__(self, /, *,
                 length=calc_seqlen,
                 perc_vec=calc_perc_vec,
                 n_batches=calc_n_batches,
                 taken=calc_time_taken,
                 speed=calc_speed,
                 **kwargs):
        # Note that the named keyword arguments must come after **kwargs
        # because they are calculated using
        super().__init__(**kwargs, length=length, perc_vec=perc_vec,
                         n_batches=n_batches, taken=taken, speed=speed)

    @classmethod
    def path_segs(cls):
        return super().path_segs() + (path.MutVecRepSeg,)

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.MOD: path.MOD_VECT}

    @classmethod
    def batch_seg(cls) -> tuple[path.Segment, str]:
        return path.MutVecBatSeg, path.ORC_EXT
