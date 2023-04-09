from __future__ import annotations
from datetime import datetime
from inspect import signature
import json
from logging import getLogger
from pathlib import Path
from typing import Any

from .batch import get_batch_dir, get_batch_path, BATCH_NUM_START
from ..util import path
from ..util.seq import DNA
from ..util.util import digest_file


logger = getLogger(__name__)


def _get_report_path(out_dir: Path, sample: str, ref: str):
    report_seg = path.VecRepSeg.build({path.EXT: path.JSON_EXT})
    return get_batch_dir(out_dir, sample, ref).joinpath(report_seg)


class VectorReport(object):
    """
    Read and write a report about a mutational profile, including:
    - the sample, reference, and section
    - number of mutation vectors
    - number of mutation vector batch files and their checksums
    - beginning and ending time, duration, and speed of vectoring

    'GTATTTTTACAACAATTACC'
    """

    class AbstractField(object):
        dtype: type
        key: str
        alias: str

        def __init__(self, value: Any):
            if not isinstance(value, self.dtype):
                raise TypeError(f"{self.__class__.__name__} expected a value "
                                f"of type '{self.dtype.__name__}' but got "
                                f"type '{type(value).__name__}'")
            self.value: Any = value

        @classmethod
        def parse(cls, value: str):
            # For the base class, just cast value to the expected type.
            return cls(cls.dtype(value))

        def __str__(self):
            return str(self.value)

    class AbstractStrField(AbstractField):
        """ String field """
        dtype = str

    class AbstractNumField(AbstractField):
        """ Numeric field with optional minimal and maximal values """
        dtype: int | float
        minimum: int | float | None = None
        maximum: int | float | None = None

        def __init__(self, value: int | float):
            super().__init__(value)
            # Verify bounds. Note that if value is NaN, the < and >
            # comparisons always return False, so the checks will pass,
            # which is the intended behavior. Do not check for
            # if not value >= self.minimum and
            # if not value <= self.maximum
            # because then the checks will fail for NaN values.
            if self.minimum is not None:
                if value < self.minimum:
                    raise ValueError(f"{self.alias} must be ≥ {self.minimum}, "
                                     f"but got {value}")
            if self.maximum is not None:
                if value > self.maximum:
                    raise ValueError(f"{self.alias} must be ≤ {self.maximum}, "
                                     f"but got {value}")

    class AbstractPosIntField(AbstractNumField):
        """ Positive integer """
        dtype, minimum = int, 1

    class AbstractNonNegIntField(AbstractNumField):
        """ Non-negative integer """
        dtype, minimum = int, 0

    class AbstractNonNegFloatField(AbstractNumField):
        """ Non-negative real number """
        dtype, minimum = float, 0.0

    class SampleField(AbstractStrField):
        key, alias = "sample", "Sample Name"

    class RefField(AbstractStrField):
        key, alias = "ref", "Reference Name"

    class LengthField(AbstractPosIntField):
        key, alias = "length", "Section Length"

    class AbstractDnaField(AbstractField):
        dtype = DNA

        @classmethod
        def parse(cls, value: str):
            # Need to encode the value to bytes before casting.
            return cls(cls.dtype(value.encode()))

        def __str__(self):
            # Need to decode the value from DNA to str.
            return self.value.decode()

    class SeqField(AbstractDnaField):
        key, alias = "seq", "Section Sequence"

    class BoolField(AbstractField):
        dtype = bool

        @classmethod
        def parse(cls, value: str):
            lower = value.strip().lower()
            if lower == "true":
                return cls(True)
            if lower == "false":
                return cls(False)
            raise ValueError(f"Cannot parse '{value}' as {cls.dtype.__name__}")

    class NumVectorsField(AbstractNonNegIntField):
        key, alias = "n_vectors", "Reads Vectorized"

    class NumReadErrorsField(AbstractNonNegIntField):
        key, alias = "n_readerr", "Reads with Errors"

    class FracVectorsField(AbstractNonNegFloatField):
        key, alias = "f_success", "Fraction Vectorized"

    class NumBatchesField(AbstractNonNegIntField):
        key, alias = "n_batches", "Batches"

    class ListStrField(AbstractField):
        dtype, delim = list, ", "

        @classmethod
        def parse(cls, value: str):
            return cls(value.strip().split(cls.delim) if value else list())

        def __str__(self):
            return self.delim.join(self.value)

    class ChecksumsField(ListStrField):
        key, alias = "checksums", "MD5 Checksums"

    class AbstractDateTimeField(AbstractField):
        dtype, dtfmt = datetime, "%Y-%m-%d %H:%M:%S.%f"

        @classmethod
        def parse(cls, value: str):
            return cls(cls.dtype.strptime(value, cls.dtfmt))

        def __str__(self):
            return self.value.strftime(self.dtfmt)

    class TimeBeganField(AbstractDateTimeField):
        key, alias = "began", "Began"

    class TimeEndedField(AbstractDateTimeField):
        key, alias = "ended", "Ended"

    class TimeTakenField(AbstractNonNegFloatField):
        key, alias = "taken", "Time taken (s)"

    class SpeedField(AbstractNonNegFloatField):
        key, alias = "speed", "Speed (1/s)"

    fields: list[type[AbstractField]] = [SampleField, RefField, LengthField,
                                         SeqField, NumVectorsField,
                                         NumReadErrorsField, FracVectorsField,
                                         NumBatchesField, ChecksumsField,
                                         TimeBeganField, TimeEndedField,
                                         TimeTakenField, SpeedField]

    fields_by_key: dict[str, type[AbstractField]] = {field.key: field
                                                     for field in fields}

    fields_by_alias: dict[str, type[AbstractField]] = {field.alias: field
                                                       for field in fields}

    @staticmethod
    def compute_length(seq: DNA):
        return len(seq)

    @staticmethod
    def compute_frac_vect(n_vectors: int, n_readerr: int):
        try:
            frac_vect = n_vectors / (n_vectors + n_readerr)
        except ZeroDivisionError:
            frac_vect = float("nan")
        return round(frac_vect, 5)

    @staticmethod
    def compute_n_batches(checksums: list[str]):
        return len(checksums)

    @staticmethod
    def compute_time_taken(began: datetime, ended: datetime):
        """ Compute time taken in seconds. """
        dtime = ended - began
        time_taken = dtime.seconds + dtime.microseconds / 1e6
        return round(time_taken, 3)

    @staticmethod
    def compute_speed(n_vectors: int, time_taken: float):
        try:
            speed = n_vectors / time_taken
        except ZeroDivisionError:
            speed = float("inf" if n_vectors else "nan")
        return round(speed, 1)

    def __init__(self, /, *, sample: str, ref: str, seq: DNA,
                 n_vectors: int, n_readerr: int, checksums: list[str],
                 began: datetime, ended: datetime):
        # Fill fields using keyword arguments.
        self[self.SampleField] = sample
        self[self.RefField] = ref
        self[self.SeqField] = seq
        self[self.NumVectorsField] = n_vectors
        self[self.NumReadErrorsField] = n_readerr
        self[self.ChecksumsField] = checksums
        self[self.TimeBeganField] = began
        self[self.TimeEndedField] = ended
        # Compute computable fields.
        self[self.LengthField] = len(seq)
        self[self.FracVectorsField] = self.compute_frac_vect(n_vectors,
                                                             n_readerr)
        self[self.NumBatchesField] = len(checksums)
        self[self.TimeTakenField] = self.compute_time_taken(began, ended)
        self[self.SpeedField] = self.compute_speed(n_vectors,
                                                   self[self.TimeTakenField])

    init_keys = [key for key, param in signature(__init__).parameters.items()
                 if param.kind == param.KEYWORD_ONLY]

    def __getitem__(self, field: type[AbstractField]):
        return self.__getattribute__(field.key)

    def __setitem__(self, field: type[AbstractField], value: Any):
        self.__setattr__(field.key, field(value).value)

    def find_invalid_batches(self, out_dir: Path, validate_checksums: bool):
        """ Return all the batches of mutation vectors that either do
        not exist or do not match their expected checksums. """
        missing = list()
        badsum = list()
        for batch, checksum in enumerate(self[self.ChecksumsField],
                                         start=BATCH_NUM_START):
            batch_path = get_batch_path(out_dir,
                                        self[self.SampleField],
                                        self[self.RefField],
                                        batch)
            if batch_path.is_file():
                if validate_checksums:
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

    def save(self, out_dir: Path):
        data = {alias: str(field(self[field]))
                for alias, field in self.fields_by_alias.items()}
        report_path = _get_report_path(out_dir,
                                       self[self.SampleField],
                                       self[self.RefField])
        report_path.parent.mkdir(parents=True, exist_ok=True)
        with open(report_path, "w") as f:
            json.dump(data, f)
        return report_path

    @classmethod
    def load(cls,
             report_file: Path,
             validate_checksums: bool = True) -> VectorReport:
        """ Load a mutation vector report from a file. """
        kwargs = dict()
        with open(report_file, "r") as f:
            for alias, value in json.load(f).items():
                field = cls.fields_by_alias[alias]
                if field.key in cls.init_keys:
                    kwargs[field.key] = field.parse(value).value
        report = cls(**kwargs)
        report_path = path.parse(report_file, path.ModSeg, path.SampSeg,
                                 path.RefSeg, path.VecRepSeg)
        if report_path[path.SAMP] != report[cls.SampleField]:
            raise ValueError(f"Samples from report ({report[cls.SampleField]}) "
                             f"and path ({report_path[path.SAMP]}) disagree")
        if report_path[path.REF] != report[cls.RefField]:
            raise ValueError(f"References from report ({report[cls.RefField]}) "
                             f"and path ({report_path[path.REF]}) disagree")
        missing, badsum = report.find_invalid_batches(report_path[path.TOP],
                                                      validate_checksums)
        if missing:
            raise FileNotFoundError(f"Batches: {', '.join(map(str, missing))}")
        if badsum:
            raise ValueError(f"Bad MD5 sums: {', '.join(map(str, badsum))}")
        return report
