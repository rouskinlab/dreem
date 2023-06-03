from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path
from typing import Iterable

import pandas as pd

from .base import (POS_FIELD, READ_FIELD, POPAVG_TITLE,
                   Table, SectTable, PosTable, PropTable, ReadTable,
                   RelPosTable, RelReadTable,
                   MaskPosTable, MaskReadTable,
                   ClustPosTable, ClustReadTable, ClustPropTable)
from ..cluster.load import CLUST_NAME_IDX
from ..core import path
from ..core.rna import RnaProfile
from ..core.sect import Section


# Table Loader Base Classes ############################################

class TableLoader(Table, ABC):
    """ Load a table from a file. """

    def __init__(self, table_file: Path):
        fields = path.parse(table_file, *self.path_segs())
        self._out_dir = fields[path.TOP]
        self._sample = fields[path.SAMP]
        self._ref = fields[path.REF]
        if not self.path.samefile(table_file):
            raise ValueError(f"Invalid path for '{self.__class__.__name__}': "
                             f"{table_file} (expected {self.path})")

    @property
    def out_dir(self):
        return self._out_dir

    @property
    def sample(self):
        return self._sample

    @property
    def ref(self):
        return self._ref

    @classmethod
    @abstractmethod
    def index_col(cls):
        return ""

    @cached_property
    def data(self):
        return pd.read_csv(self.path, index_col=self.index_col())


class SectTableLoader(TableLoader, SectTable, ABC):
    def __init__(self, table_file: Path):
        fields = path.parse(table_file, *self.path_segs())
        self._sect = fields[path.SECT]
        super().__init__(table_file)

    @property
    def sect(self):
        return self._sect


# Load by Index (position/read/proportion) #############################

class PosTableLoader(TableLoader, PosTable, ABC):
    """ Load data indexed by position. """
    @classmethod
    def index_col(cls):
        return POS_FIELD


class ReadTableLoader(TableLoader, ReadTable, ABC):
    """ Load data indexed by read. """
    @classmethod
    def index_col(cls):
        return READ_FIELD


class PropTableLoader(TableLoader, PropTable, ABC):
    """ Load cluster proportions. """
    @classmethod
    def index_col(cls):
        return CLUST_NAME_IDX


# Load by Source and Index #############################################

class SectPosTableLoader(SectTableLoader, PosTableLoader, ABC):
    """ Load sectioned data indexed by position. """
    @abstractmethod
    def iter_profiles(self, sections: Iterable[Section]):
        """ Yield RNA mutational profiles from the table. """
        for section in sections:
            yield RnaProfile("", section, "", "", pd.Series())


class SectReadTableLoader(SectTableLoader, ReadTableLoader, ABC):
    """ Load sectioned data indexed by read. """


class SectPropTableLoader(SectTableLoader, PropTableLoader, ABC):
    """ Load sectioned data indexed by cluster. """


# Instantiable Table Loaders ###########################################

class RelPosTableLoader(PosTableLoader, RelPosTable):
    """ Load relation data indexed by position. """


class RelReadTableLoader(ReadTableLoader, RelReadTable):
    """ Load relation data indexed by read. """


class MaskPosTableLoader(SectPosTableLoader, MaskPosTable):
    """ Load masked bit vector data indexed by position. """

    def iter_profiles(self, sections: Iterable[Section]):
        for section in sections:
            yield RnaProfile(title=path.fill_whitespace(POPAVG_TITLE),
                             section=section,
                             sample=self.sample,
                             data_sect=self.sect,
                             reacts=self.fm)


class MaskReadTableLoader(SectReadTableLoader, MaskReadTable):
    """ Load masked bit vector data indexed by read. """


class ClusterPosTableLoader(SectPosTableLoader, ClustPosTable):
    """ Load cluster data indexed by position. """

    def iter_profiles(self, sections: Iterable[Section]):
        """ Yield RNA mutational profiles from a table. """
        for section in sections:
            for cluster in self.cluster_names:
                yield RnaProfile(title=path.fill_whitespace(cluster),
                                 section=section,
                                 sample=self.sample,
                                 data_sect=self.sect,
                                 reacts=self.data[cluster])


class ClusterReadTableLoader(SectReadTableLoader, ClustReadTable):
    """ Load cluster data indexed by read. """


class ClusterPropTableLoader(SectPropTableLoader, ClustPropTable):
    """ Load cluster data indexed by cluster. """


# Helper Functions #####################################################

def load(table_file: Path):
    """ Helper function to load a TableLoader from a table file. """
    for loader_type in (RelPosTableLoader,
                        RelReadTableLoader,
                        MaskPosTableLoader,
                        MaskReadTableLoader,
                        ClusterPosTableLoader,
                        ClusterReadTableLoader,
                        ClusterPropTableLoader):
        try:
            # Try to load the table file using the type of TableLoader.
            return loader_type(table_file)
        except (FileNotFoundError, ValueError):
            # The TableLoader was not the correct type for the file.
            pass
    # None of the TableLoader types were able to open the file.
    raise ValueError(f"Failed to open table: {table_file}")
