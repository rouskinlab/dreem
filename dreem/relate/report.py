from __future__ import annotations
from logging import getLogger

from ..core import path
from ..core.report import (Report, calc_seqlen, calc_time_taken,
                           calc_perc_rel_pass, calc_n_batches, calc_speed)

logger = getLogger(__name__)


class RelateReport(Report):
    __slots__ = ("sample", "ref", "seq", "length",
                 "n_reads_rel_pass", "n_reads_rel_fail", "perc_rel_pass", "checksums",
                 "n_batches", "began", "ended", "taken", "speed")

    def __init__(self, /, *,
                 length=calc_seqlen,
                 perc_rel_pass=calc_perc_rel_pass,
                 n_batches=calc_n_batches,
                 taken=calc_time_taken,
                 speed=calc_speed,
                 **kwargs):
        # Note that the named keyword arguments must come after **kwargs
        # because they are calculated using
        super().__init__(**kwargs, length=length, perc_rel_pass=perc_rel_pass,
                         n_batches=n_batches, taken=taken, speed=speed)

    @classmethod
    def path_segs(cls):
        return super().path_segs() + (path.RelateRepSeg,)

    @classmethod
    def auto_fields(cls):
        return {**super().auto_fields(), path.MOD: path.MOD_REL}

    @classmethod
    def batch_seg(cls) -> tuple[path.Segment, str]:
        return path.RelateBatSeg, path.ORC_EXT
