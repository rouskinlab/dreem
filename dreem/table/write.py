from abc import ABC, abstractmethod
from logging import getLogger
from pathlib import Path

import pandas as pd

from .base import (POS_FIELD, READ_FIELD, SEQ_FIELD, Table, SectTable,
                   RelTable, MaskTable, ClustTable,
                   PosTable, ReadTable, PropTable,
                   RelPosTable, RelReadTable, MaskPosTable, MaskReadTable,
                   ClustPosTable, ClustReadTable, ClustPropTable)
from .tabulate import (Tabulator, RelTabulator, MaskTabulator, ClustTabulator,
                       tabulate_loader)
from ..cluster.load import ClustLoader
from ..core import path
from ..mask.load import MaskLoader
from ..relate.load import RelateLoader

logger = getLogger(__name__)


# Table Writer Base Classes ############################################

class TableWriter(Table, ABC):
    """ Write a table to a file. """

    def __init__(self, tabulator: Tabulator):
        self._tab = tabulator

    @property
    def out_dir(self):
        return self._tab.out_dir

    @property
    def sample(self):
        return self._tab.sample

    @property
    def ref(self):
        return self._tab.ref

    @abstractmethod
    def load_data(self):
        """ Load the table's data from a DataLoader. """
        return pd.DataFrame()

    @classmethod
    @abstractmethod
    def write_precision(cls):
        """ Number of digits to keep behind the decimal place. """
        return 0

    @property
    def data(self):
        return self.load_data()

    def write(self, rerun: bool):
        """ Write the table's rounded data to the table's CSV file. """
        # Check if the File exists.
        if rerun or not self.path.is_file():
            self.data.round(self.write_precision()).to_csv(self.path)
        else:
            logger.warning(f"File exists: {self.path}")
        return self.path


class SectTableWriter(TableWriter, SectTable, ABC):
    """ Write a table associated with a section. """

    @property
    def sect(self):
        return self._tab.sect


# Write by Source (relate/mask/cluster) ################################

class RelTableWriter(TableWriter, RelTable, ABC):
    """ Write a table of relation vectors. """

    @classmethod
    def write_precision(cls):
        return 0


class MaskTableWriter(SectTableWriter, MaskTable, ABC):
    """ Write a table of masked bit vectors. """

    @classmethod
    def write_precision(cls):
        return 1


class ClustTableWriter(SectTableWriter, ClustTable, ABC):
    """ Write a table of clusters. """

    @classmethod
    def write_precision(cls):
        return 6


# Write by Index (position/read/cluster) ###############################

class PosTableWriter(TableWriter, PosTable, ABC):

    def load_data(self):
        # Load the data for each position, including excluded positions.
        data = self._tab.tabulate_by_pos().reindex(index=self._tab.index)
        # Replace the base-position formatted index with numeric format.
        data.index = pd.Index(self._tab.positions, name=POS_FIELD)
        # Insert the sequence into the first column of the data frame.
        data.insert(0, SEQ_FIELD, pd.Series(list(self._tab.seq.decode()),
                                            index=data.index))
        return data


class ReadTableWriter(TableWriter, ReadTable, ABC):

    def load_data(self):
        # Load the data for each read.
        data = self._tab.tabulate_by_read()
        # Rename the index.
        data.index.rename(READ_FIELD, inplace=True)
        return data


class PropTableWriter(TableWriter, PropTable, ABC):

    def load_data(self):
        # Load the data for each cluster.
        data = self._tab.tabulate_by_clust()
        # Rename the one column.
        data.columns = [CLUST_PROP_COL]
        return data


# Write by Source and Index ############################################

class SectPosTableLoader(SectTableWriter, PosTableWriter, ABC):
    """ Load sectioned data indexed by position. """


# Instantiable Table Writers ###########################################

class RelPosTableWriter(RelTableWriter, PosTableWriter, RelPosTable):
    pass


class RelReadTableWriter(RelTableWriter, ReadTableWriter, RelReadTable):
    pass


class MaskPosTableWriter(MaskTableWriter, PosTableWriter, MaskPosTable):
    pass


class MaskReadTableWriter(MaskTableWriter, ReadTableWriter, MaskReadTable):
    pass


class ClustPosTableWriter(ClustTableWriter, PosTableWriter, ClustPosTable):

    @classmethod
    def clusters_on_columns(cls):
        return True


class ClustReadTableWriter(ClustTableWriter, ReadTableWriter, ClustReadTable):

    @classmethod
    def clusters_on_columns(cls):
        return True


class ClustPropTableWriter(ClustTableWriter, PropTableWriter, ClustPropTable):

    @classmethod
    def clusters_on_columns(cls):
        return False


# Helper Functions #####################################################

def infer_report_loader_type(report_file: Path):
    """ Given a report file path, infer the type of Loader it needs. """
    if path.RelateRepSeg.ptrn.match(report_file.name):
        return RelateLoader
    if path.MaskRepSeg.ptrn.match(report_file.name):
        return MaskLoader
    if path.ClustRepSeg.ptrn.match(report_file.name):
        return ClustLoader
    raise ValueError(f"Failed to infer loader for {report_file}")


def get_tabulator_writer_types(tabulator: Tabulator):
    if isinstance(tabulator, RelTabulator):
        return RelPosTableWriter, RelReadTableWriter
    if isinstance(tabulator, MaskTabulator):
        return MaskPosTableWriter, MaskReadTableWriter
    if isinstance(tabulator, ClustTabulator):
        return ClustPosTableWriter, ClustReadTableWriter, ClustPropTableWriter
    raise TypeError(f"Invalid tabulator type: {type(tabulator).__name__}")


def get_tabulator_writers(tabulator: Tabulator):
    for writer_type in get_tabulator_writer_types(tabulator):
        yield writer_type(tabulator)


def write(report_file: Path, rerun: bool):
    """ Helper function to write a table from a report file. """
    # Determine the needed type of report loader.
    report_loader_type = infer_report_loader_type(report_file)
    # Load the report.
    report_loader = report_loader_type.open(report_file)
    # Create the tabulator for the report's data.
    tabulator = tabulate_loader(report_loader)
    # For each table associated with this tabulator, create the table,
    # write it, and return the path to the table output file.
    return [table.write(rerun) for table in get_tabulator_writers(tabulator)]
