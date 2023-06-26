from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path
from typing import Iterable

import pandas as pd

from .base import (POS_TITLE, READ_TITLE, POPAVG_TITLE, MUTAT_REL,
                   Table, PosTable, PropTable, ReadTable,
                   RelPosTable, RelReadTable,
                   MaskPosTable, MaskReadTable,
                   ClustPosTable, ClustReadTable, ClustPropTable)
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
        self._sect = fields[path.SECT]
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

    @property
    def sect(self):
        return self._sect

    @classmethod
    @abstractmethod
    def index_col(cls):
        return ""

    @cached_property
    def data(self):
        return pd.read_csv(self.path, index_col=self.index_col())


# Load by Index (position/read/proportion) #############################

class PosTableLoader(TableLoader, PosTable, ABC):
    """ Load data indexed by position. """
    @classmethod
    def index_col(cls):
        return POS_TITLE

    @abstractmethod
    def iter_profiles(self, sections: Iterable[Section]):
        """ Yield RNA mutational profiles from the table. """
        for section in sections:
            yield RnaProfile("", section, "", "", pd.Series())


class ReadTableLoader(TableLoader, ReadTable, ABC):
    """ Load data indexed by read. """
    @classmethod
    def index_col(cls):
        return READ_TITLE


class PropTableLoader(TableLoader, PropTable, ABC):
    """ Load cluster proportions. """
    @classmethod
    def index_col(cls):
        return CLUST_NAME_IDX


# Instantiable Table Loaders ###########################################

class RelPosTableLoader(PosTableLoader, RelPosTable):
    """ Load relation data indexed by position. """

    def iter_profiles(self, sections: Iterable[Section]):
        # Relation table loaders contain unmasked, unfiltered reads and
        # are thus unsuitable for generating RNA profiles.
        yield from ()


class RelReadTableLoader(ReadTableLoader, RelReadTable):
    """ Load relation data indexed by read. """


class MaskPosTableLoader(PosTableLoader, MaskPosTable):
    """ Load masked bit vector data indexed by position. """

    def iter_profiles(self, sections: Iterable[Section]):
        for section in sections:
            yield RnaProfile(title=path.fill_whitespace(POPAVG_TITLE),
                             section=section,
                             sample=self.sample,
                             data_sect=self.sect,
                             reacts=self._get_rel_frac(MUTAT_REL))


class MaskReadTableLoader(ReadTableLoader, MaskReadTable):
    """ Load masked bit vector data indexed by read. """


class ClusterPosTableLoader(PosTableLoader, ClustPosTable):
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


class ClusterReadTableLoader(ReadTableLoader, ClustReadTable):
    """ Load cluster data indexed by read. """


class ClusterPropTableLoader(PropTableLoader, ClustPropTable):
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
