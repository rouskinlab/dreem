from __future__ import annotations

from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path
from typing import Any, Iterable

import numpy as np

from . import path
from .report import (OutDirF, SampleF, RefF, SectF, End5F, End3F,
                     Report, BatchReport)
from .sect import Section, seq_pos_to_index
from .seq import DNA


class DataLoader(ABC):
    """ Base class for loading data from a step of DREEM. """

    def __init__(self, report: Report | BatchReport, **kwargs):
        if not isinstance(report, self.get_report_type()):
            raise TypeError(f"Expected report type {self.get_report_type()}, "
                            f"but got {type(report)}")
        self._report: Report | BatchReport = report
        self._kwargs = kwargs

    @classmethod
    @abstractmethod
    def get_report_type(cls):
        """ Type of the report for this Loader. """
        return Report

    @classmethod
    def _open_report(cls, report_file: Path):
        """ Open a report of the correct type for this class. """
        return cls.get_report_type().open(report_file)

    @property
    def out_dir(self) -> Path:
        """ Output directory. """
        return self._report.get_field(OutDirF)

    @property
    def sample(self) -> str:
        """ Name of sample. """
        return self._report.get_field(SampleF)

    @property
    def ref(self) -> str:
        """ Name of reference. """
        return self._report.get_field(RefF)

    @abstractmethod
    def get_refseq(self):
        """ Sequence of the reference. """
        return DNA(b"")

    @property
    def seq(self):
        """ Sequence of the reference. """
        return self.get_refseq()

    @property
    def positions(self):
        """ Numeric positions of the reference sequence. """
        return np.arange(1, len(self.seq) + 1)

    @property
    def index(self):
        """ Base-position indexes of the reference sequence. """
        return seq_pos_to_index(self.seq, self.positions, start=1)

    @abstractmethod
    def _load_data_private(self, data_file: Path):
        """ Load a private dataset. """

    def clone(self, **kwargs):
        """ Return a new DataLoader with updated keyword arguments. """
        return self.__class__(self._report, **{**self._kwargs, **kwargs})

    @classmethod
    def open(cls, report_file: Path, **kwargs):
        """ Create a new DataLoader from a report file. """
        return cls(cls._open_report(report_file), **kwargs)

    def __str__(self):
        return f"{self.__class__.__name__} for '{self.sample}' on '{self.ref}'"

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return self._report == other._report and self._kwargs == other._kwargs


class BatchLoader(DataLoader, ABC):
    """ Load a dataset that is split into batches. """

    @classmethod
    @abstractmethod
    def get_report_type(cls):
        """ Type of the report for this Loader. """
        return BatchReport

    def build_batch_path(self, batch: int):
        """ Path to the batch with the given number. """
        return self._report.build_batch_path(self.out_dir, batch,
                                             **self._report.path_fields())

    def _iter_batches_private(self):
        """ Yield each batch of private data. """
        for batch_file in self._report.iter_batch_paths():
            yield self._load_data_private(batch_file)

    def iter_batches_public(self):
        """ Yield every batch of public data. """
        # For the base BatchLoader class, all data are public, so the
        # "private" data are really public data and are yielded from
        # this method. The public/private distinction matters for the
        # BatchChainLoader derived class. """
        yield from self._iter_batches_private()


class ChainLoader(DataLoader, ABC):
    """ Load data via a DataLoader from a previous step. """

    @classmethod
    @abstractmethod
    def get_import_type(cls):
        """ Type of the data loader that is immediately before this data
        loader in the chain and from which data are thus imported. """

    @property
    @abstractmethod
    def import_kwargs(self):
        """ Keyword arguments for creating the imported data loader. """
        return dict()

    @property
    def import_path(self):
        """ Fields for creating the imported data loader. """
        return {path.SAMP: self.sample, path.REF: self.ref}

    @cached_property
    def import_loader(self):
        """ Data loader that is immediately before this data loader in
        the chain and from which data are thus imported. """
        import_type = self.get_import_type()
        report_type = import_type.get_report_type()
        return import_type.open(report_type.build_path(self.out_dir,
                                                       **self.import_path),
                                **self.import_kwargs)


class BatchChainLoader(ChainLoader, BatchLoader, ABC):
    """ Load data via a BatchLoader from a previous step. """

    @abstractmethod
    def _publish_batch(self, private_batch: Any, imported_batch: Any,
                       modifier: Any = None):
        """ Given one batch of imported data and the corresponding batch
        of private data, return one batch of public data. """

    def _publish_batch_mods(self, private_batch: Any, imported_batch: Any,
                            modifiers: Iterable):
        """ Given one batch of imported data, the corresponding batch of
        private data, and an iterable of modifiers for _publish_batch,
        yield one batch of public data for each modifier. """
        for modifier in modifiers:
            yield self._publish_batch(private_batch, imported_batch, modifier)

    def iter_batches_import(self):
        """ Yield every imported batch. """
        yield from self.import_loader.iter_batches_public()

    def _iter_batches_private_import(self):
        """ Yield every pair of private and imported batches. """
        yield from zip(self._iter_batches_private(), self.iter_batches_import(),
                       strict=True)

    def iter_batches_public(self):
        for private, imported in self._iter_batches_private_import():
            yield self._publish_batch(private, imported)

    def iter_batches_modified(self, modifiers: Iterable):
        for private, imported in self._iter_batches_private_import():
            yield self._publish_batch_mods(private, imported, modifiers)


class SectLoader(DataLoader, ABC):
    """ Load data from a report that has a section. """

    @property
    def sect(self) -> str:
        """ Name of the section. """
        return self._report.get_field(SectF)

    @property
    def end5(self) -> int:
        """ 5' end of the section. """
        return self._report.get_field(End5F)

    @property
    def end3(self) -> int:
        """ 3' end of the section. """
        return self._report.get_field(End3F)

    @cached_property
    def section(self):
        """ Section of the reference. """
        return Section(ref=self.ref, refseq=self.get_refseq(),
                       end5=self.end5, end3=self.end3, name=self.sect)

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
    def import_path(self):
        if issubclass(self.get_import_type(), SectLoader):
            # Add the section to the fields for the upstream loader.
            return {**super().import_path, path.SECT: self.sect}
        return super().import_path
