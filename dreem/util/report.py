from __future__ import annotations

from abc import ABC, abstractmethod
from datetime import datetime
from inspect import getmembers
import json
from functools import cache
from logging import getLogger
from math import isclose, isnan, inf, nan
from pathlib import Path
import re
import sys
from typing import Any, Hashable, Callable, Iterable

from .path import STR_CHARS_SET
from .seq import DNA

logger = getLogger(__name__)


# Field class

def identity(value: Any):
    """ Identity function. """
    return value


def truthy(*_: Any):
    """ Truthy function. """
    return True


class Field(object):
    __slots__ = ["key", "title", "dtype", "iconv", "oconv",
                 "check_val", "check_rep_val"]

    def __init__(self, /, key: str, title: str, dtype: type, *,
                 iconv: Callable[[Any], Any] | None = None,
                 oconv: Callable[[Any], Any] | None = None,
                 check_val: Callable[[Any], bool] | None = None,
                 check_rep_val: Callable[[Report, Any], bool] | None = None):
        """
        Parameters
        ----------
        key: str
            Key under which the field is stored as an attribute of the
            Report instance.
        title: str
            Title by which the field is identified in the report file.
        dtype: type
            Data type of the field
        iconv: Callable[[Any], Any] | None = None
            Convert an input value to the right type for instantiation.
            If omitted, just cast the value to the proper data type.
        oconv: Callable[[Any], Any] | None = None
            Convert an output value to the right type for exporting.
            If omitted, export the value with no change.
        check_val: Callable[[Any], bool] | None = None
            Validate the value of the field (after validating its type)
            upon assigning the value to an attribute of a Report object.
            Must accept one argument (the value) and return True if the
            value is valid and False otherwise. If None, the value is
            not validated -- only its data type.
        check_rep_val: Callable[[Report, Any], bool] | None = None
            Validate the value of the field (after validating its type)
            upon assigning the value to an attribute of a Report object.
            Must accept two arguments (the Report object and the value)
            and return True if the value is valid and False otherwise.
            If None, the value is not validated -- only its data type.
        """
        self.key = key
        self.title = title
        self.dtype = dtype
        self.iconv = self.dtype if iconv is None else iconv
        self.oconv = identity if oconv is None else oconv
        self.check_val = truthy if check_val is None else check_val
        self.check_rep_val = truthy if check_rep_val is None else check_rep_val

    def validate(self, report: Report, value: Any):
        # Validate the type.
        if not isinstance(value, self.dtype):
            raise TypeError(f"{self} expected value to be {self.dtype}, "
                            f"but got {type(value)}")
        # Validate the value.
        if not self.check_val(value):
            raise ValueError(f"{self} got invalid value: {value}")
        if not self.check_rep_val(report, value):
            raise ValueError(f"{self} got invalid value: {value}")

    def __str__(self):
        return f"Report Field '{self.title}' ({self.key})"


# Field calculation functions

# Note that each function takes a single argument: a Report instance.
# So these functions could be implemented as Report instance methods.
# But this implementation could cause confusion because no one Report
# class can work with all these methods.

def calc_time_taken(report: Report) -> float:
    """ Calculate the time taken in minutes. """
    delta = report.get_field(TimeEndedF) - report.get_field(TimeBeganF)
    minutes = (delta.seconds + 1e-6 * delta.microseconds) / 60.
    if minutes < 0.:
        raise ValueError(f"Time taken must be positive, but got {minutes} min")
    return minutes


def calc_seqlen(report: Report):
    return len(report.get_field(SeqF))


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


# Field value checking functions

def agrees(num1: float, num2: float, precision: int | float):
    """ Check if two floats agree within a given decimal precision. """
    return isclose(num1, num2, abs_tol=10. ** -precision, rel_tol=0.)


def nanagrees(num1: float, num2: float, precision: int | float):
    """ Like agrees, but also return True if both floats are NaN. """
    return (isnan(num1) and isnan(num2)) or agrees(num1, num2, precision)


def check_name(name: str):
    return bool(name) and not bool(set(name) - STR_CHARS_SET)


def check_nonneg_int(num: int):
    return isinstance(num, int) and num >= 0


def check_pos_int(num: int):
    return isinstance(num, int) and num > 0


def check_nonneg_float(num: float):
    return isinstance(num, float) and num >= 0.


def check_nannonneg_float(num: float):
    return isinstance(num, float) and num >= 0. or isnan(num)


def check_bool(x: bool):
    return isinstance(x, bool)


def check_float(num: float):
    return isinstance(num, float) and not isnan(num)


def check_pos_float(num: float):
    return isinstance(num, float) and num > 0.


def check_probability(probability: float):
    return isinstance(probability, float) and 0. <= probability <= 1.


def check_nanprobability(probability: float):
    return isinstance(probability, float) and (0. <= probability <= 1.
                                               or isnan(probability))


def check_percentage(percentage: float):
    return isinstance(percentage, float) and 0. <= percentage <= 100.


def check_nanpercentage(percentage: float):
    return isinstance(percentage, float) and (0. <= percentage <= 100.
                                              or isnan(percentage))


def check_position(report: Report, pos: int):
    return check_pos_int(pos) and pos <= len(report.get_field(SeqF))


def check_positions(report: Report, positions: list[int]):
    return (isinstance(positions, list)
            and all(check_position(report, pos) for pos in positions))


def check_end5(report: Report, pos: int):
    return check_position(report, pos)


def check_end3(report: Report, end3: int):
    return end3 == report.get_field(End5F) + len(report.get_field(SeqF)) - 1


def check_seqlen(report: Report, seqlen: int):
    return check_pos_int(seqlen) and seqlen == len(report.get_field(SeqF))


def check_perc_vec(report: Report, perc_vec: float):
    return (check_nanpercentage(perc_vec) and agrees(perc_vec,
                                                     calc_perc_vec(report),
                                                     PERC_VEC_PRECISION))


def check_checksums(checksums: list[str]):
    checksum_pattern = re.compile("^[0-9A-Fa-f]{32}$")
    return (isinstance(checksums, list)
            and all(isinstance(checksum, str) for checksum in checksums)
            and all(map(checksum_pattern.match, checksums)))


def check_n_batches(report: Report, n_batches):
    return check_nonneg_int(n_batches) and n_batches == calc_n_batches(report)


def check_time_ended(report: Report, ended: datetime):
    return (isinstance(ended, datetime)
            and ended >= report.get_field(TimeBeganF))


def check_time_taken(report: Report, taken: float):
    return (check_nonneg_float(taken) and agrees(taken,
                                                 calc_time_taken(report),
                                                 TIME_TAKEN_PRECISION))


def check_speed(report: Report, speed: float):
    return (check_nannonneg_float(speed) and nanagrees(speed,
                                                       calc_speed(report),
                                                       SPEED_PRECISION))


def check_ints_range(ints: Iterable[int]):
    if not all(isinstance(num, int) for num in ints):
        return False
    ints_sort = sorted(ints)
    return ints_sort == list(range(ints_sort[0], ints_sort[-1] + 1))


def check_cluster_dict(cdict: dict[int, Any]):
    return (isinstance(cdict, dict)
            and check_ints_range(cdict.keys())
            and all(map(check_pos_int, cdict.keys())))


def check_clusts_list_nonneg_int(lis: dict[int, list[int]]):
    return (check_cluster_dict(lis)
            and all(isinstance(li, list) for li in lis.values())
            and all(all(map(check_nonneg_int, li)) for li in lis.values()))


def check_clusts_floats(floats: dict[int, float]):
    return (check_cluster_dict(floats)
            and all(map(check_float, floats.values())))


def check_clusts_list_floats(lfs: dict[int, list[float]]):
    return (check_cluster_dict(lfs)
            and all(isinstance(lf, list) for lf in lfs.values())
            and all(all(map(check_float, lf)) for lf in lfs.values()))


# Field definitions

DATETIME_FORMAT = "%Y-%m-%d at %H:%M:%S"
DECIMAL_PRECISION = 5  # general precision for decimals
PERC_VEC_PRECISION = 2
TIME_TAKEN_PRECISION = 2
SPEED_PRECISION = -2


@cache
def get_oconv_float(precision: int = DECIMAL_PRECISION):
    def oconv_float(num: float):
        return round(num, precision)

    return oconv_float


@cache
def get_oconv_list_float(precision: int = DECIMAL_PRECISION):
    oconv_float = get_oconv_float(precision)

    def oconv_list_float(nums: list[float]) -> list[float]:
        return list(map(oconv_float, nums))
    return oconv_list_float


@cache
def get_oconv_dict_float(precision: int = DECIMAL_PRECISION):
    oconv_float = get_oconv_float(precision)

    def oconv_dict_float(dnum: dict[Hashable, float]) -> dict[Hashable, float]:
        return {d: oconv_float(num) for d, num in dnum.items()}
    return oconv_dict_float


@cache
def get_oconv_dict_list_float(precision: int = DECIMAL_PRECISION):
    oconv_list_float = get_oconv_list_float(precision)

    def oconv_dict_list_float(dnums: dict[Hashable, list[float]]
                              ) -> dict[Hashable, list[float]]:
        return {d: oconv_list_float(nums) for d, nums in dnums.items()}
    return oconv_dict_list_float


def iconv_datetime(text: str):
    return datetime.strptime(text, DATETIME_FORMAT)


def oconv_datetime(dtime: datetime):
    return dtime.strftime(DATETIME_FORMAT)


# General
SampleF = Field("sample", "Name of Sample", str, check_val=check_name)
RefF = Field("ref", "Name of Reference", str, check_val=check_name)
SeqF = Field("seq", "Sequence of Reference/Section", DNA,
             iconv=DNA.parse, oconv=DNA.__str__)
SectF = Field("sect", "Name of Section", str, check_val=check_name)
End5F = Field("end5", "5' end of Section", int, check_rep_val=check_end5)
End3F = Field("end3", "3' end of Section", int, check_rep_val=check_end3)
SeqLenF = Field("length", "Length of Sequence (nt)", int,
                check_rep_val=check_seqlen)
TimeBeganF = Field("began", "Time Began", datetime,
                   iconv=iconv_datetime, oconv=oconv_datetime)
TimeEndedF = Field("ended", "Time Ended", datetime,
                   iconv=iconv_datetime, oconv=oconv_datetime,
                   check_rep_val=check_time_ended)
TimeTakenF = Field("taken", "Time Taken (minutes)", float,
                   oconv=get_oconv_float(TIME_TAKEN_PRECISION),
                   check_rep_val=check_time_taken)

# Mutation vector generation
NumVecF = Field("n_vectors", "Number of Reads Vectorized", int,
                check_val=check_nonneg_int)
NumErrF = Field("n_readerr", "Number of Reads with Errors", int,
                check_val=check_nonneg_int)
PercVecF = Field("perc_vec", "Percent of Reads Vectorized (%)", float,
                 oconv=get_oconv_float(PERC_VEC_PRECISION),
                 check_val=check_percentage)
NumBatchF = Field("n_batches", "Number of Batches of Vectors", int,
                  check_val=check_nonneg_int)
ChecksumsF = Field("checksums", "MD5 Checksums", list,
                   check_val=check_checksums)
SpeedF = Field("speed", "Speed (Vectors per minute)", float,
               oconv=get_oconv_float(SPEED_PRECISION),
               check_rep_val=check_speed)

# Positional filtering
CountDelF = Field("count_del", "Count Deletions as Mutations", bool)
CountInsF = Field("count_ins", "Count Insertions as Mutations", bool)
ExclPolyAF = Field("exclude_polya",
                   "Minimum Length of Poly(A) Sequences to Exclude (nt)",
                   int, check_val=check_nonneg_int)
ExclGUF = Field("exclude_gu", "Exclude G/U Bases", bool)
ExclUserPosF = Field("exclude_pos", "Exclude User-Defined Positions",
                     list, check_rep_val=check_positions)

# Bit vector filtering
MinInfoPosF = Field("min_ninfo_pos",
                    "Minimum Number of Informative Reads per Position",
                    int, check_val=check_nonneg_int)
MinMutPosF = Field("min_fmut_pos",
                   "Minimum Fraction of Mutations per Position",
                   float, oconv=get_oconv_float(),
                   check_val=check_probability)
MaxMutPosF = Field("max_fmut_pos",
                   "Maximum Fraction of Mutations per Position",
                   float, oconv=get_oconv_float(),
                   check_val=check_probability)
MinMutGapF = Field("min_mut_gap",
                   "Minimum Gap Between Mutations (nt)",
                   int, check_val=check_nonneg_int)
MinInfoReadF = Field("min_finfo_read",
                     "Minimum Fraction of Informative Positions per Read",
                     float, oconv=get_oconv_float(),
                     check_val=check_probability)
MaxMutReadF = Field("max_fmut_read",
                    "Maximum Fraction of Mutations per Read",
                    float, oconv=get_oconv_float(),
                    check_val=check_probability)
NumPosInitF = Field("n_pos_init",
                    "Number of Positions in Total",
                    int, check_val=check_nonneg_int)
NumPosCutPolyAF = Field("n_pos_polya",
                        "Number of Positions Cut -- Poly(A) Sequence",
                        int, check_val=check_nonneg_int)
NumPosCutGUF = Field("n_pos_gu",
                     "Number of Positions Cut -- G/U Base",
                     int, check_val=check_nonneg_int)
NumPosCutUserF = Field("n_pos_user",
                       "Number of Positions Cut -- User-Defined",
                       int, check_val=check_nonneg_int)
NumPosCutLoInfoF = Field("n_min_ninfo_pos",
                         "Number of Positions Cut -- Too Few Informative Reads",
                         int, check_val=check_nonneg_int)
NumPosCutLoMutF = Field("n_min_fmut_pos",
                        "Number of Positions Cut -- Too Few Mutations",
                        int, check_val=check_nonneg_int)
NumPosCutHiMutF = Field("n_max_fmut_pos",
                        "Number of Positions Cut -- Too Many Mutations",
                        int, check_val=check_nonneg_int)
NumPosKeptF = Field("n_pos_kept",
                    "Number of Positions Kept",
                    int, check_val=check_nonneg_int)
NumReadInitF = Field("n_reads_init",
                     "Number of Reads in Total",
                     int, check_val=check_nonneg_int)
NumReadLoInfoF = Field("n_min_finfo_read",
                       "Number of Reads Cut -- Too Few Informative Positions",
                       int, check_val=check_nonneg_int)
NumReadHiMutF = Field("n_max_fmut_read",
                      "Number of Reads Cut -- Too Many Mutations",
                      int, check_val=check_nonneg_int)
NumReadCloseMutF = Field("n_min_mut_gap",
                         "Number of Reads Cut -- Mutations Too Close Together",
                         int, check_val=check_nonneg_int)
NumReadKeptF = Field("n_reads_kept",
                     "Number of Reads Kept",
                     int, check_val=check_nonneg_int)
NumUniqReadKeptF = Field("n_uniq_reads_kept",
                         "Number of Unique Reads Kept",
                         int, check_val=check_nonneg_int)

# EM clustering
MinIterClustF = Field("min_iter",
                      "Minimum Number of EM Iterations",
                      int, check_val=check_nonneg_int)
MaxIterClustF = Field("max_iter",
                      "Maximum Number of EM Iterations",
                      int, check_val=check_pos_int)
ClustConvThreshF = Field("conv_thresh",
                         "Convergence Threshold for Log Likelihood",
                         float, oconv=get_oconv_float(),
                         check_val=check_pos_float)
MaxClustsF = Field("max_clust",
                   "Maximum Attempted Number of Clusters",
                   int, check_val=check_pos_int)
ClustNumRunsF = Field("num_runs",
                      "Number of Independent EM Runs",
                      int, check_val=check_pos_int)
NumClustsF = Field("n_clust",
                   "Optimal Number of Clusters",
                   int, check_val=check_pos_int)
ClustsBicF = Field("bic",
                   "Bayesian Information Criterion of Best EM Run",
                   dict, oconv=get_oconv_dict_float(),
                   check_val=check_clusts_floats)
ClustsConvF = Field("converged",
                    "Iterations Until Convergence for Each EM Run",
                    dict, check_val=check_clusts_list_nonneg_int)
ClustsLogLikesF = Field("log_likes",
                        "Log Likelihood of Each EM Run",
                        dict, oconv=get_oconv_dict_list_float(),
                        check_val=check_clusts_list_floats)
ClustsLikeMeanF = Field("log_like_mean",
                        "Log Likelihood Mean Among EM Runs",
                        dict, oconv=get_oconv_dict_float(),
                        check_val=check_clusts_floats)
ClustsLikeStdF = Field("log_like_std",
                       "Log Likelihood Std. Dev. Among EM Runs",
                       dict, oconv=get_oconv_dict_float(),
                       check_val=check_clusts_floats)
ClustsVarInfoF = Field("var_info",
                       "Variation of Information Among EM Runs",
                       dict, oconv=get_oconv_dict_float(),
                       check_val=check_clusts_floats)


# Field managing functions

@cache
def fields() -> list[Field]:
    return [member for _, member in getmembers(sys.modules[__name__])
            if isinstance(member, Field)]


@cache
def field_keys() -> dict[str, Field]:
    return {field.key: field for field in fields()}


@cache
def field_titles() -> dict[str, Field]:
    return {field.title: field for field in fields()}


def lookup_key(key: str) -> Field:
    try:
        return field_keys()[key]
    except KeyError:
        raise KeyError(f"No Report Field is keyed '{key}'")


def lookup_title(title: str) -> Field:
    try:
        return field_titles()[title]
    except KeyError:
        raise KeyError(f"No field is titled '{title}'")


# Report classes

class Report(ABC):
    """ Abstract base class for reports from steps in DREEM. """
    __slots__ = []

    def __init__(self, **kwargs: Any | Callable[[Report], Any]):
        # Get the values of all attributes (specified in __slots__) from
        # the keyword arguments.
        for key in self.__slots__:
            # Get the value of the keyword argument.
            kwarg_val = kwargs.pop(key)
            if callable(kwarg_val):
                # If the value of the keyword argument is callable, then
                # it must accept one argument -- self -- and return the
                # value of the attribute.
                attr_val = kwarg_val(self)
            else:
                # If not, then use the value of the keyword argument.
                attr_val = kwarg_val
            self.__setattr__(key, attr_val)
        # Confirm that no unexpected keywords were given.
        if kwargs:
            raise TypeError("Got unexpected keyword arguments for "
                            f"{self.__class__.__name__}: {kwargs}")

    @abstractmethod
    def get_path(self, out_dir: Path) -> Path:
        raise NotImplementedError

    def get_field(self, field: Field):
        """ Return the value of a field of the report using the field
        instance directly, not its key. """
        return self.__getattribute__(field.key)

    def to_dict(self):
        """ Return a dict of raw values of the fields, keyed by the
        titles of their fields. """
        odata = dict()
        for key in self.__slots__:
            field = lookup_key(key)
            odata[field.title] = field.oconv(self.get_field(field))
        return odata

    def save(self, out_dir: Path):
        report_path = self.get_path(out_dir)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Created directory {report_path.parent}")
        with open(report_path, "w") as f:
            json.dump(self.to_dict(), f, indent=4)
        logger.debug(f"Wrote {self} to {report_path}")
        return report_path

    @classmethod
    def from_dict(cls, odata: dict[str, Any]):
        """ Convert a dict of raw values (keyed by the titles of their
        fields) into a dict of encoded values (keyed by the keys of
        their fields), from which a new Report is instantiated. """
        if not isinstance(odata, dict):
            raise TypeError("Report classmethod from_data expected 'dict',"
                            f"but got '{type(odata).__name__}'")
        # Read every raw value, keyed by the title of its field.
        idata = dict()
        for title, value in odata.items():
            # Get the field corresponding to the title.
            field = lookup_title(title)
            # Cast the value to the input type and key it by the field.
            idata[field.key] = field.iconv(value)
        # Instantiate and return a new Report from the values.
        return cls(**idata)

    @classmethod
    def open(cls, json_file: Path):
        """ Create a new Report instance from a JSON file. """
        with open(json_file) as f:
            return cls.from_dict(json.load(f))

    def __setattr__(self, key, value):
        """ Validate the attribute name and value before setting it. """
        lookup_key(key).validate(self, value)
        super().__setattr__(key, value)
