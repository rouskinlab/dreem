from abc import ABC, abstractmethod
from pathlib import Path

from .report import Report, OutDirF, SampleF, RefF, SectF, End5F, End3F
from .sect import Section
from .seq import DNA


class DataLoader(ABC):
    """ Base class for loading data from a step of DREEM. """

    def __init__(self, report: Report):
        if not isinstance(report, self.report_type()):
            raise TypeError(f"Expected report of type {self.report_type()}, "
                            f"but got {type(report)}")
        self._rep = report

    @classmethod
    @abstractmethod
    def report_type(cls):
        """ Type of the report for this Loader. """
        return Report

    @property
    def out_dir(self) -> Path:
        """ Output directory. """
        return self._rep.get_field(OutDirF)

    @property
    def sample(self) -> str:
        """ Name of sample. """
        return self._rep.get_field(SampleF)

    @property
    def ref(self) -> str:
        """ Name of reference. """
        return self._rep.get_field(RefF)

    @property
    @abstractmethod
    def sect(self) -> str | None:
        """ Name of section, or None if the report has no section. """
        try:
            return self._rep.get_field(SectF)
        except AttributeError:
            return None

    @property
    @abstractmethod
    def section(self):
        """ Section of the reference, or the full reference if the report
        has no section. """
        return Section(ref="", ref_seq=DNA(b""))

    @property
    @abstractmethod
    def seq(self):
        """ Sequence of the section, or of the full reference if the
        report has no section. """
        return DNA(b"")

    @property
    def end5(self):
        """ 5' end of section, or 1 if no section. """
        try:
            return self._rep.get_field(End5F)
        except AttributeError:
            return 1

    @property
    def end3(self):
        """ 3' end of section, or len(seq) if no section. """
        try:
            return self._rep.get_field(End3F)
        except AttributeError:
            return len(self.seq)

    @property
    def positions(self):
        """ Numeric positions of the section. """
        return self.section.positions

    @property
    def index(self):
        """ Base-position indexes of the section. """
        return self.section.index

    @classmethod
    def open(cls, report_file: Path):
        """ Create a new Loader from a report file. """
        return cls(cls.report_type().open(report_file))

    def __str__(self):
        return f"Mutation Vectors from '{self.sample}' aligned to '{self.ref}'"
