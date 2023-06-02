from __future__ import annotations

from abc import ABC, abstractmethod
from functools import cached_property
from logging import getLogger
from pathlib import Path

import pandas as pd

from .base import (CLUST_PROP_COL, CLUST_PROP_IDX,
                   DELET_FIELD, INSRT_FIELD, MATCH_FIELD, MUTAT_FIELD,
                   SUB_A_FIELD, SUB_C_FIELD, SUB_G_FIELD, SUB_T_FIELD,
                   SUB_N_FIELD, TOTAL_FIELD,
                   POS_FIELD, READ_FIELD, SEQ_FIELD,
                   Table,
                   CountTable, SectTable, RelTable, MaskTable, ClustTable,
                   PosTable, ReadTable, PropTable,
                   RelPosTable, RelReadTable, MaskPosTable, MaskReadTable,
                   ClustPosTable, ClustReadTable, ClustPropTable)
from ..mask.load import MaskLoader
from ..cluster.load import ClustLoader
from ..core import path
from ..core.bitc import BitCaller, SemiBitCaller
from ..core.bitv import BitCounter
from ..relate.load import RelateLoader

logger = getLogger(__name__)


# Tabulator Classes ####################################################

class Tabulator(ABC):
    """ Base class for tabulating data for multiple tables from a report
    loader. """

    def __init__(self, loader: RelateLoader | MaskLoader | ClustLoader):
        self._load = loader

    @property
    def out_dir(self):
        return self._load.out_dir

    @property
    def sample(self):
        return self._load.sample

    @property
    def ref(self):
        return self._load.ref

    @property
    def sect(self):
        return self._load.sect

    @property
    def seq(self):
        return self._load.seq

    @property
    def positions(self):
        return self._load.positions

    @property
    def index(self):
        return self._load.index

    @abstractmethod
    def tabulate_by_pos(self):
        """ Compute the data for each position. """
        return pd.DataFrame()

    @abstractmethod
    def tabulate_by_read(self):
        """ Compute the data for each read / bit vector. """
        return pd.DataFrame()

    @abstractmethod
    def list_writer_types(self) -> list[type[TableWriter]]:
        """ List the TableWriter classes for this tabulator. """
        return list()

    def get_writers(self):
        """ Yield each TableWriter object for this tabulator. """
        return (writer_type(self) for writer_type in self.list_writer_types())

    @staticmethod
    def from_loader(report_loader: RelateLoader | MaskLoader | ClustLoader):
        """ Return a new DataLoader, choosing the subclass based on the
        type of the argument `report_loader`. """
        if isinstance(report_loader, RelateLoader):
            return RelTabulator(report_loader)
        if isinstance(report_loader, MaskLoader):
            return MaskTabulator(report_loader)
        if isinstance(report_loader, ClustLoader):
            return ClustTabulator(report_loader)
        raise TypeError(f"Invalid loader type: {type(report_loader).__name__}")


class CountTabulator(Tabulator, ABC):
    """ Base class for constructing data for multiple count-based tables
    from a report loader. """

    @classmethod
    def iter_bit_callers(cls):
        """ Yield every BitCaller, one for each type of information to
        include in the table of counts. """
        # Initialize two semi-bit-callers for convenience.
        ref_caller = SemiBitCaller.from_counts(count_ref=True,
                                               cache_all=True)
        mut_caller = SemiBitCaller.from_counts(count_sub=True,
                                               count_del=True,
                                               count_ins=True,
                                               cache_all=True)
        # Count all base calls (everything but the bytes 0 and 255).
        yield TOTAL_FIELD, BitCaller(SemiBitCaller(),
                                     SemiBitCaller.from_counts(count_ref=True,
                                                               count_sub=True,
                                                               count_del=True,
                                                               count_ins=True))
        # Count matches to the reference sequence.
        yield MATCH_FIELD, BitCaller(mut_caller, ref_caller)
        # Count all types of mutations, relative to reference matches.
        yield MUTAT_FIELD, BitCaller(ref_caller, mut_caller)
        # Count each type of mutation, relative to reference matches.
        yield SUB_N_FIELD, BitCaller(ref_caller,
                                     SemiBitCaller.from_counts(count_sub=True))
        yield SUB_A_FIELD, BitCaller(ref_caller,
                                     SemiBitCaller("ca", "ga", "ta"))
        yield SUB_C_FIELD, BitCaller(ref_caller,
                                     SemiBitCaller("ac", "gc", "tc"))
        yield SUB_G_FIELD, BitCaller(ref_caller,
                                     SemiBitCaller("ag", "cg", "tg"))
        yield SUB_T_FIELD, BitCaller(ref_caller,
                                     SemiBitCaller("at", "ct", "gt"))
        yield DELET_FIELD, BitCaller(ref_caller,
                                     SemiBitCaller.from_counts(count_del=True))
        yield INSRT_FIELD, BitCaller(ref_caller,
                                     SemiBitCaller.from_counts(count_ins=True))

    @cached_property
    def bit_counters(self) -> dict[str, BitCounter]:
        """ Return a BitCounter for the batches of relation vectors. """
        # Collect all the bit callers.
        callers = dict(self.iter_bit_callers())
        # Create a bit counter for each bit caller.
        counters = {field: BitCounter() for field in callers}
        # Load each batch of relation vectors.
        for rel_batch in self._load.iter_rel_batches():
            # Use each bit caller to create bit vectors from this batch.
            for field, caller in callers.items():
                bit_batch = caller.call(rel_batch)
                # Count the bits in the batch with its bit counter.
                counters[field].add_batch(bit_batch)
        # Prevent each bit counter from accepting more data.
        for counter in counters.values():
            counter.close()
        return counters

    def tabulate_by_pos(self):
        """ DataFrame of the bit count for each position and caller. """
        return pd.DataFrame.from_dict({field: counter.nmuts_per_pos
                                       for field, counter
                                       in self.bit_counters.items()})

    def tabulate_by_read(self):
        """ DataFrame of the bit count for each read and caller. """
        return pd.DataFrame.from_dict({field: counter.nmuts_per_read
                                       for field, counter
                                       in self.bit_counters.items()})


class RelTabulator(CountTabulator):
    def list_writer_types(self):
        return [RelPosTableWriter, RelReadTableWriter]


class MaskTabulator(CountTabulator):
    def list_writer_types(self):
        return [MaskPosTableWriter, MaskReadTableWriter]


class ClustTabulator(Tabulator):
    def list_writer_types(self):
        return [ClustPosTableWriter, ClustReadTableWriter, ClustPropTableWriter]

    def tabulate_by_pos(self):
        """ Mutation rates of each position in each cluster. """
        return self._load.mus

    def tabulate_by_read(self):
        """ Responsibilities of each read in each cluster. """
        return self._load.resps

    def tabulate_by_clust(self):
        """ Proportion of each cluster. """
        return self._load.props.to_frame()


# Table Writer Base Classes ############################################

class TableWriter(Table, ABC):
    """ Write a table to a file. """

    def __init__(self, tabulator: Tabulator):
        self._out_dir = tabulator.out_dir
        self._sample = tabulator.sample
        self._ref = tabulator.ref
        self._data = self.load_data(tabulator)

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
    def load_data(cls, _: Tabulator):
        return pd.DataFrame()

    @classmethod
    @abstractmethod
    def write_precision(cls):
        """ Number of digits to keep behind the decimal place. """
        return 0

    @property
    def data(self):
        return self._data

    def write(self, rerun: bool):
        """ Write the table's rounded data to the table's CSV file. """
        # Check if the table already exists.
        if rerun or not self.path.is_file():
            self.data.round(self.write_precision()).to_csv(self.path)
        else:
            logger.warning(f"Table already exists: {self.path}")
        return self.path


class CountTableWriter(TableWriter, CountTable, ABC):
    """ Write a table that counts bits of multiple types (base calls,
    matches, and each type of mutation). """


class SectTableWriter(TableWriter, SectTable, ABC):
    """ Write a table associated with a section. """

    def __init__(self, tabulator: Tabulator):
        super().__init__(tabulator)
        self._sect = tabulator.sect

    @property
    def sect(self):
        return self._sect


# Write by Source (relate/mask/cluster) ################################

class RelTableWriter(CountTableWriter, RelTable, ABC):
    """ Write a table of relation vectors. """

    @classmethod
    def write_precision(cls):
        return 0


class MaskTableWriter(SectTableWriter, CountTableWriter, MaskTable, ABC):
    """ Write a table of masked bit vectors. """

    @classmethod
    def write_precision(cls):
        return 2


class ClustTableWriter(SectTableWriter, ClustTable, ABC):
    """ Write a table of clusters. """

    @classmethod
    def write_precision(cls):
        return 5

    @classmethod
    @abstractmethod
    def clusters_on_columns(cls):
        """ Whether the columns are the indexes with the clusters. """
        return False

    @classmethod
    def load_data(cls, tabulator: Tabulator):
        data = super().load_data(tabulator)
        # Replace the cluster index with the cluster names.
        if cls.clusters_on_columns():
            data.columns = pd.Index(cls.format_names(data.columns))
        else:
            data.index = pd.Index(cls.format_names(data.index))
        return data


# Write by Index (position/read/cluster) ###############################

class PosTableWriter(TableWriter, PosTable, ABC):

    @classmethod
    def load_data(cls, tabulator: Tabulator):
        # Load the data for each position, including excluded positions.
        data = tabulator.tabulate_by_pos().reindex(index=tabulator.index)
        # Replace the base-position formatted index with numeric format.
        data.index = pd.Index(tabulator.positions, name=POS_FIELD)
        # Insert the sequence into the first column of the data frame.
        data.insert(0, SEQ_FIELD, pd.Series(list(tabulator.seq.decode()),
                                            index=data.index))
        return data


class ReadTableWriter(TableWriter, ReadTable, ABC):

    @classmethod
    def load_data(cls, tabulator: Tabulator):
        # Load the data for each read.
        data = tabulator.tabulate_by_read()
        # Rename the index.
        data.index.rename(READ_FIELD, inplace=True)
        return data


class PropTableWriter(TableWriter, PropTable, ABC):

    @classmethod
    def load_data(cls, tabulator: ClustTabulator):
        # Load the data for each cluster.
        data = tabulator.tabulate_by_clust()
        # Rename the index.
        data.index.rename(CLUST_PROP_IDX)
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


def write(report_file: Path, rerun: bool):
    """ Helper function to write a table from a report file. """
    # Determine the needed type of report loader.
    report_loader_type = infer_report_loader_type(report_file)
    # Load the report.
    report_loader = report_loader_type.open(report_file)
    # Create the tabulator for the report's data.
    tabulator = Tabulator.from_loader(report_loader)
    # For each table associated with this tabulator, create the table,
    # write it, and return the path to the table output file.
    return [table.write(rerun) for table in tabulator.get_writers()]
