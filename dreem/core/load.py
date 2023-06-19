from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path
from typing import Any

import numpy as np

from . import path
from .report import (OutDirF, SampleF, RefF, SectF, End5F, End3F,
                     Report, BatchReport)
from .sect import Section, seq_pos_to_index
from .seq import DNA


class DataLoader(ABC):
    """ Base class for loading data from a step of DREEM. """

    def __init__(self, report: Report | BatchReport):
        if not isinstance(report, self.get_report_type()):
            raise TypeError(f"Expected report type {self.get_report_type()}, "
                            f"but got {type(report)}")
        self._rep: Report | BatchReport = report

    @classmethod
    @abstractmethod
    def get_report_type(cls):
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
    def seq(self):
        """ Sequence of the reference. """
        return DNA(b"")

    @property
    def positions(self):
        """ Numeric positions of the reference sequence. """
        return np.arange(1, len(self.seq) + 1)

    @property
    def index(self):
        """ Base-position indexes of the reference sequence. """
        return seq_pos_to_index(self.seq, self.positions, start=1)

    @abstractmethod
    def load_data(self, *args, **kwargs):
        """ Load a dataset. """

    @classmethod
    def open(cls, report_file: Path):
        """ Create a new DataLoader from a report file. """
        return cls(cls.get_report_type().open(report_file))

    def __str__(self):
        return f"{self.__class__.__name__} for '{self.sample}' on '{self.ref}'"


class BatchLoader(DataLoader, ABC):
    """ Load a dataset that is split into batches. """

    @classmethod
    @abstractmethod
    def get_report_type(cls):
        """ Type of the report for this Loader. """
        return BatchReport

    def build_batch_path(self, batch: int):
        """ Path to the batch with the given number. """
        return self._rep.build_batch_path(self.out_dir, batch,
                                          **self._rep.path_fields())

    def iter_batches(self, *args, **kwargs):
        """ Iterate through all batches of data. """
        for batch_path in self._rep.iter_batch_paths():
            yield self.load_data(batch_path, *args, **kwargs)


class ChainLoader(DataLoader, ABC):
    """ Load data via a DataLoader from a previous step. """

    @classmethod
    @abstractmethod
    def get_upstream_type(cls):
        """ Type of the upstream data loader in the chain. """
        return DataLoader

    @property
    def _upstream_fields(self):
        """ Fields for creating the upstream data loader. """
        return {path.SAMP: self.sample, path.REF: self.ref}

    @cached_property
    def _upload(self) -> DataLoader | BatchLoader:
        """ Upstream data loader in the chain. """
        ups_type = self.get_upstream_type()
        rep_type = ups_type.get_report_type()
        return ups_type.open(rep_type.build_path(self.out_dir,
                                                 **self._upstream_fields))


class BatchChainLoader(ChainLoader, BatchLoader, ABC):
    """ Load data via a BatchLoader from a previous step. """

    @classmethod
    @abstractmethod
    def get_upstream_type(cls):
        """ Type of the upstream data loader in the chain. """
        return BatchLoader

    def _iter_upstream_batches(self, **kwargs):
        """ Iterate through the batches of the upstream BatchLoader. """
        yield from self._upload.iter_batches(**kwargs)

    def iter_report_batches(self, **kwargs):
        """ Iterate through the batches that the report knows about. """
        yield from super().iter_batches(**kwargs)

    @abstractmethod
    def _process_batch(self, upstream_batch: Any, report_batch: Any):
        """ Given one batch of data from the upstream BatchLoader and
        the corresponding batch from the report of this BatchLoader,
        generate a new batch of processed data. """

    def iter_batches(self, *args, **kwargs):
        """ Override `self.iter_batches()` and yield processed batches
        instead of batches straight from the report. """
        upstream_batches = self._iter_upstream_batches(**kwargs)
        report_batches = self.iter_report_batches(**kwargs)
        yield from map(self._process_batch, upstream_batches, report_batches)


class SectLoader(DataLoader, ABC):
    """ Load data from a report that has a section. """

    @property
    def sect(self) -> str:
        """ Name of the section. """
        return self._rep.get_field(SectF)

    @property
    def end5(self) -> int:
        """ 5' end of section. """
        return self._rep.get_field(End5F)

    @property
    def end3(self) -> int:
        """ 3' end of section. """
        return self._rep.get_field(End3F)

    @abstractmethod
    def _get_refseq(self):
        """ Sequence of the full reference. """
        return DNA(b"")

    @cached_property
    def section(self):
        """ Section of the reference. """
        return Section(ref=self.ref, refseq=self._get_refseq(), name=self.sect,
                       end5=self.end5, end3=self.end3)

    @property
    def seq(self):
        """ Sequence of the section. """
        return self.section.seq

    @property
    def positions(self):
        """ Numeric positions of the section. """
        return self.section.positions

    @property
    def index(self):
        """ Base-position indexes of the section. """
        return self.section.index

    def __str__(self):
        return f"{super().__str__()}:{self.sect}"


class SectBatchChainLoader(BatchChainLoader, SectLoader, ABC):

    @property
    def _upstream_fields(self):
        if issubclass(self.get_upstream_type(), SectLoader):
            # Add the section to the fields for the upstream loader.
            return {**super()._upstream_fields, path.SECT: self.sect}
        return super()._upstream_fields
