from __future__ import annotations

from abc import ABC
from logging import getLogger
from math import inf, nan
from pathlib import Path

from .batch import get_batch_dir, get_batch_path, BATCH_NUM_START
from ..util import path
from ..util.report import (Report, calc_seqlen, calc_time_taken,
                           ChecksumsF, NumErrF, NumVecF, TimeTakenF)
from ..util.files import digest_file

logger = getLogger(__name__)


def get_report_path(out_dir: Path, sample: str, ref: str):
    report_seg = path.VecRepSeg.build(ext=path.JSON_EXT)
    return get_batch_dir(out_dir, sample, ref).joinpath(report_seg)


# Field calculation functions

def calc_n_batches(report: Report):
    return len(report.get_field(ChecksumsF))


def calc_perc_vec(report: Report) -> float:
    nvecs = report.get_field(NumVecF)
    nerrs = report.get_field(NumErrF)
    try:
        return 100. * nvecs / (nvecs + nerrs)
    except ZeroDivisionError:
        return nan


def calc_speed(report: Report) -> float:
    nvecs = report.get_field(NumVecF)
    taken = report.get_field(TimeTakenF)
    try:
        return nvecs / taken
    except ZeroDivisionError:
        return inf if nvecs > 0 else nan


class VectorReport(Report, ABC):
    __slots__ = ["sample", "ref", "seq", "length",
                 "n_vectors", "n_readerr", "perc_vec", "checksums",
                 "n_batches", "began", "ended", "taken", "speed"]

    def __init__(self, /, *,
                 length=calc_seqlen,
                 perc_vec=calc_perc_vec,
                 n_batches=calc_n_batches,
                 taken=calc_time_taken,
                 speed=calc_speed,
                 **kwargs):
        super().__init__(length=length, perc_vec=perc_vec,
                         n_batches=n_batches, taken=taken, speed=speed,
                         **kwargs)

    def get_path(self, out_dir: Path):
        return get_report_path(out_dir, self.sample, self.ref)

    def find_invalid_batches(self, out_dir: Path):
        """ Return all the batches of mutation vectors that either do
        not exist or do not match their expected checksums. """
        missing = list()
        badsum = list()
        for batch, checksum in enumerate(self.checksums, start=BATCH_NUM_START):
            batch_path = get_batch_path(out_dir, self.sample, self.ref, batch)
            if batch_path.is_file():
                check = digest_file(batch_path)
                if check != checksum:
                    # The file exists does not match the checksum.
                    badsum.append(batch_path)
                    logger.critical(f"Expected {batch_path} to have MD5 "
                                    f"checksum {checksum}, but got {check}")
            else:
                # The batch file does not exist.
                missing.append(batch_path)
        return missing, badsum

    @classmethod
    def open(cls, json_file: Path) -> VectorReport:
        report = super().open(json_file)
        report_path = path.parse(json_file, path.ModSeg, path.SampSeg,
                                 path.RefSeg, path.VecRepSeg)
        missing, badsum = report.find_invalid_batches(report_path[path.TOP])
        if missing:
            raise FileNotFoundError(", ".join(map(str, missing)))
        if badsum:
            raise ValueError(f"Bad MD5 sums: {', '.join(map(str, badsum))}")
        return report
