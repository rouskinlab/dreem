from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path
import re
from typing import Any, Iterable

import pandas as pd

from ..call.load import BitVecLoader
from ..cluster.load import ClusterLoader
from ..core import path
from ..core.bit import SemiBitCaller, BitCaller, BitCounter
from ..core.mu import calc_mu_df
from ..core.seq import DNA
from ..relate.load import RelVecLoader

# Definitions ##########################################################

# General fields
POS_FIELD = "Position"
SEQ_FIELD = "Base"
READ_FIELD = "Read Name"

# Relation vector fields
BASES_FIELD = "All Base Calls"
MATCH_FIELD = "Matches"
MUTAT_FIELD = "All Muts"
DELET_FIELD = "Deletions"
INSRT_FIELD = "Insertions"
SUB_N_FIELD = "All Subs"
SUB_A_FIELD = "Subs to A"
SUB_C_FIELD = "Subs to C"
SUB_G_FIELD = "Subs to G"
SUB_T_FIELD = "Subs to T"

# Bit vector fields
NINFO_FIELD = "Num Info"
NMUTS_FIELD = "Num Muts"
FINFO_FIELD = "Frac Info"
FMUTS_FIELD = "Frac Muts"
FMUTO_FIELD = "Frac Muts (Obs)"
FMUTA_FIELD = "Frac Muts (Adj)"

# Cluster fields
CLUST_FORMAT = "Cluster {k}-{c}"
CLUST_PATTERN = re.compile("Cluster ([0-9]+)-([0-9]+)")
CLUST_PROP_IDX = "Cluster"
CLUST_PROP_COL = "Proportion"

# Constants
PRECISION = 6  # number of digits to keep in floating-point numbers


# Table Base Class #####################################################

class Table(ABC):
    @property
    @abstractmethod
    def out_dir(self) -> Path:
        """ Table's output directory. """
        return Path()

    @property
    @abstractmethod
    def sample(self) -> str:
        """ Table's sample name. """
        return ""

    @property
    @abstractmethod
    def ref(self) -> str:
        """ Table's reference name. """
        return ""

    @property
    @abstractmethod
    def sect(self) -> str | None:
        """ Table's section name. """
        return None

    @cached_property
    @abstractmethod
    def data(self) -> pd.DataFrame:
        """ Table's data frame. """
        return pd.DataFrame()

    @classmethod
    @abstractmethod
    def kind(cls) -> str:
        """ Kind of table. """
        return ""

    @classmethod
    def has_sect(cls):
        """ Whether the table is associated with a section. """
        return issubclass(cls, BitVecTable) or issubclass(cls, ClusterTable)

    @classmethod
    def by_read(cls):
        """ Whether the table contains data for each read. """
        return issubclass(cls, ReadTable)

    @classmethod
    def path_segs(cls):
        """ Table's path segments. """
        segs = [path.ModSeg, path.SampSeg, path.RefSeg]
        if cls.has_sect():
            segs.append(path.SectSeg)
        segs.append(path.MutTabSeg)
        return tuple(segs)

    @classmethod
    def gzipped(cls):
        """ Whether the table's file is compressed with gzip. """
        return cls.by_read()

    @classmethod
    def ext(cls):
        """ Table's file extension: either '.csv' or '.csv.gz'. """
        return path.CSVZIP_EXT if cls.gzipped() else path.CSV_EXT

    @property
    def path_fields(self) -> dict[str, Any]:
        """ Table's path fields. """
        fields = {path.TOP: self.out_dir,
                  path.MOD: path.MOD_TABLE,
                  path.SAMP: self.sample,
                  path.REF: self.ref,
                  path.TABLE: self.kind(),
                  path.EXT: self.ext()}
        if self.has_sect():
            fields[path.SECT] = self.sect
        return fields

    @cached_property
    def path(self):
        """ Path of the table's CSV file (possibly gzipped). """
        return path.buildpar(*self.path_segs(), **self.path_fields)


# Type (RelVec/BitVec/Cluster) Base Classes ############################

class RelVecTable(Table, ABC):
    @property
    def bases(self):
        return self.data[BASES_FIELD]

    @property
    def match(self):
        return self.data[MATCH_FIELD]

    @property
    def mutat(self):
        return self.data[MUTAT_FIELD]

    @property
    def delet(self):
        return self.data[DELET_FIELD]

    @property
    def insrt(self):
        return self.data[INSRT_FIELD]

    @property
    def sub_n(self):
        return self.data[SUB_N_FIELD]

    @property
    def sub_a(self):
        return self.data[SUB_A_FIELD]

    @property
    def sub_c(self):
        return self.data[SUB_C_FIELD]

    @property
    def sub_g(self):
        return self.data[SUB_G_FIELD]

    @property
    def sub_t(self):
        return self.data[SUB_T_FIELD]


class BitVecTable(Table, ABC):
    @property
    def ninfo(self):
        return self.data[NINFO_FIELD]

    @property
    def finfo(self):
        return self.data[FINFO_FIELD]

    @property
    def nmuts(self):
        return self.data[NMUTS_FIELD]


class ClusterTable(Table, ABC):
    @classmethod
    def format_names(cls, kc_pairs: Iterable[tuple[int, int]]):
        return [CLUST_FORMAT.format(k=k, c=c) for k, c in kc_pairs]

    @classmethod
    def parse_names(cls, names: Iterable[str]):
        return [(int((kc := CLUST_PATTERN.match(n).groups())[0]), int(kc[1]))
                for n in names]

    @cached_property
    @abstractmethod
    def cluster_names(self) -> list[str]:
        return list()

    @cached_property
    def cluster_tuples(self):
        """ List the clusters as tuples of (n_clusters, cluster). """
        return self.parse_names(self.cluster_names)


# Dimension (Position/Read) Base Classes ###############################

class PosTable(Table, ABC):
    @property
    def pos(self):
        return self.data.index.values

    @property
    def bases(self):
        return self.data[SEQ_FIELD]

    @cached_property
    def seq(self):
        return DNA("".join(self.bases).encode())


class ReadTable(Table, ABC):
    @property
    def reads(self):
        return self.data.index.values


# Operation (Writer/Loader) Base Classes ###############################

class TableWriter(Table, ABC):
    """ Abstract base class for writing data tables to files. """

    def __init__(self, loader: RelVecLoader | BitVecLoader | ClusterLoader):
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
        return self._load.sect if self.has_sect() else None

    @abstractmethod
    def make_data(self):
        return pd.DataFrame()

    @abstractmethod
    def index_data(self, data: pd.DataFrame):
        return data

    @cached_property
    def data(self):
        return self.index_data(self.make_data())

    def write(self):
        """ Write the table's rounded data to the table's CSV file. """
        self.data.round(PRECISION).to_csv(self.path)
        return self.path


class TableLoader(Table, ABC):
    """ Abstract base class for reading data tables from files. """

    def __init__(self, table_file: Path):
        fields = path.parse(table_file, *self.path_segs())
        self._out_dir = fields[path.TOP]
        self._sample = fields[path.SAMP]
        self._ref = fields[path.REF]
        self._sect = fields[path.SECT] if self.has_sect() else None
        if table_file != self.path:
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

    @cached_property
    def data(self):
        return pd.read_csv(self.path)


# Type + Dimension Derived Classes #####################################

class RelVecPosTable(RelVecTable, PosTable, ABC):
    @classmethod
    def kind(cls):
        return path.RELVEC_POS_TAB


class RelVecReadTable(RelVecTable, ReadTable, ABC):
    @classmethod
    def kind(cls):
        return path.RELVEC_READ_TAB


class BitVecPosTable(BitVecTable, PosTable, ABC):
    @classmethod
    def kind(cls):
        return path.BITVEC_POS_TAB

    @property
    def fmuts_obs(self):
        return self.data[FMUTO_FIELD]

    @property
    def fmuts(self):
        return self.data[FMUTA_FIELD]


class BitVecReadTable(BitVecTable, ReadTable, ABC):
    @classmethod
    def kind(cls):
        return path.BITVEC_READ_TAB

    @property
    def fmuts(self):
        return self.data[FMUTS_FIELD]


class ClusterPosTable(ClusterTable, PosTable, ABC):
    @classmethod
    def kind(cls):
        return path.CLUST_MUS_TAB

    @cached_property
    def cluster_names(self):
        return self.data.columns.drop(SEQ_FIELD).to_list()


class ClusterReadTable(ClusterTable, ReadTable, ABC):
    @classmethod
    def kind(cls):
        return path.CLUST_RESP_TAB

    @cached_property
    def cluster_names(self):
        return self.data.columns.to_list()


class ClusterPropTable(ClusterTable, Table, ABC):
    @classmethod
    def kind(cls):
        return path.CLUST_PROP_TAB

    @cached_property
    def cluster_names(self):
        return self.data.index.to_list()


# Type + Operation Derived Classes #####################################

class RelVecTableWriter(RelVecTable, TableWriter, ABC):
    @classmethod
    def iter_mut_callers(cls):
        yield SUB_N_FIELD, SemiBitCaller.from_counts(count_sub=True)
        yield SUB_A_FIELD, SemiBitCaller("ca", "ga", "ta")
        yield SUB_C_FIELD, SemiBitCaller("ac", "gc", "tc")
        yield SUB_G_FIELD, SemiBitCaller("ag", "cg", "tg")
        yield SUB_T_FIELD, SemiBitCaller("at", "ct", "gt")
        yield DELET_FIELD, SemiBitCaller.from_counts(count_del=True)
        yield INSRT_FIELD, SemiBitCaller.from_counts(count_ins=True)

    @classmethod
    def iter_bit_callers(cls):
        ref_caller = SemiBitCaller.from_counts(count_ref=True)
        for field, mut_caller in cls.iter_mut_callers():
            yield field, BitCaller(ref_caller, mut_caller)

    def iter_bit_counters(self):
        for field, caller in self.iter_bit_callers():
            yield field, BitCounter(caller, self._load.iter_batches())

    def iter_counts(self):
        # Count all base calls.
        bc = BitCounter(BitCaller(SemiBitCaller.from_counts(count_ref=True,
                                                            count_sub=True,
                                                            count_del=True,
                                                            count_ins=True),
                                  SemiBitCaller.from_counts()),
                        self._load.iter_batches())
        yield BASES_FIELD, (bc.ninfo_per_vec if self.by_read()
                            else bc.ninfo_per_pos)
        # Count matches and all mutations.
        bc = BitCounter(BitCaller(SemiBitCaller.from_counts(count_ref=True),
                                  SemiBitCaller.from_counts(count_sub=True,
                                                            count_del=True,
                                                            count_ins=True)),
                        self._load.iter_batches())
        if self.by_read():
            yield MATCH_FIELD, bc.ninfo_per_vec - bc.nmuts_per_vec
            yield MUTAT_FIELD, bc.nmuts_per_vec
        else:
            yield MATCH_FIELD, bc.ninfo_per_pos - bc.nmuts_per_pos
            yield MUTAT_FIELD, bc.nmuts_per_pos
        # Count each type of mutation.
        for field, bc in self.iter_bit_counters():
            yield field, (bc.nmuts_per_vec if self.by_read()
                          else bc.nmuts_per_pos)

    def make_data(self):
        return pd.DataFrame.from_dict(dict(self.iter_counts()))


# Dimension + Operation Derived Classes ################################

class PosTableWriter(PosTable, TableWriter, ABC):
    def index_data(self, data: pd.DataFrame):
        # Insert the sequence into the first column of the data frame.
        data.insert(0, SEQ_FIELD, pd.Series(list(self._load.seq.decode()),
                                            index=data.index))
        # Replace the base-position formatted index with numeric format.
        data.index = pd.Index(self._load.positions, name=POS_FIELD)
        return data


class ReadTableWriter(ReadTable, TableWriter, ABC):
    def index_data(self, data: pd.DataFrame):
        # Rename the index.
        data.index.rename(READ_FIELD, inplace=True)
        return data


# Type + Dimension + Operation Concrete Classes ########################

class RelVecPosTableWriter(RelVecPosTable,
                           RelVecTableWriter,
                           PosTableWriter):
    pass


class RelVecReadTableWriter(RelVecReadTable,
                            RelVecTableWriter,
                            ReadTableWriter):
    pass


class BitVecPosTableWriter(BitVecPosTable, PosTableWriter):
    def make_data(self):
        # Assemble the data from the information and mutation counts.
        counter = self._load.get_bit_counter()
        data = pd.DataFrame.from_dict({
            NINFO_FIELD: counter.ninfo_per_pos,
            FINFO_FIELD: counter.finfo_per_pos,
            NMUTS_FIELD: counter.nmuts_per_pos,
            FMUTO_FIELD: counter.fmuts_per_pos,
            FMUTA_FIELD: calc_mu_df(counter.fmuts_per_pos.to_frame(),
                                    self._load.section,
                                    self._load.min_mut_gap).squeeze(axis=1)
        })
        # Add all excluded positions in the section back to the index.
        return data.reindex(index=self._load.section.index)


class BitVecReadTableWriter(BitVecReadTable, ReadTableWriter):
    def make_data(self):
        counter = self._load.get_bit_counter()
        return pd.DataFrame.from_dict({
            NINFO_FIELD: counter.ninfo_per_vec,
            FINFO_FIELD: counter.finfo_per_vec,
            NMUTS_FIELD: counter.nmuts_per_vec,
            FMUTS_FIELD: counter.fmuts_per_vec,
        })


class ClusterPosTableWriter(ClusterPosTable, PosTableWriter):
    def make_data(self):
        # Add all excluded positions in the section back to the index.
        return self._load.mus.reindex(index=self._load.section.index)

    def index_data(self, data: pd.DataFrame):
        # Replace the columns with the cluster names.
        data.columns = pd.Index(self.format_names(data.columns))
        # Complete indexing.
        return super().index_data(data)

    def fmuts(self, n_clusters, cluster) -> pd.Series:
        """ Return the adjusted fraction of mutations for cluster
        `cluster` with `n_clusters` total clusters. """
        return self.data[CLUST_FORMAT.format(k=n_clusters, c=cluster)]


class ClusterReadTableWriter(ClusterReadTable, ReadTableWriter):
    def make_data(self):
        return self._load.resps.T

    def index_data(self, data: pd.DataFrame):
        # Replace the columns with the cluster names.
        data.columns = pd.Index(self.format_names(data.columns))
        # Complete indexing.
        return super().index_data(data)


class ClusterPropTableWriter(ClusterPropTable, TableWriter):
    def make_data(self):
        return self._load.props.to_frame()

    def index_data(self, data: pd.DataFrame):
        # Replace the index with the cluster names.
        data.index = pd.Index(self.format_names(data.index),
                              name=CLUST_PROP_IDX)
        # Replace the column name.
        data.columns = [CLUST_PROP_COL]
        return data


class BitVecPosTableLoader(BitVecPosTable, TableLoader):
    pass


class BitVecReadTableLoader(BitVecReadTable, TableLoader):
    pass


# Helper Functions #####################################################

def infer_report_loader_type(report_file: Path):
    """ Given a report file path, infer the type of Loader it needs. """
    if path.RelateRepSeg.ptrn.match(report_file.name):
        return RelVecLoader
    if path.CallRepSeg.ptrn.match(report_file.name):
        return BitVecLoader
    if path.ClustRepSeg.ptrn.match(report_file.name):
        return ClusterLoader
    raise ValueError(f"Failed to infer loader for {report_file}")


def get_table_types(loader_type: type):
    if loader_type is RelVecLoader:
        return (RelVecPosTableWriter,
                RelVecReadTableWriter)
    if loader_type is BitVecLoader:
        return (BitVecPosTableWriter,
                BitVecReadTableWriter)
    if loader_type is ClusterLoader:
        return (ClusterPosTableWriter,
                ClusterReadTableWriter,
                ClusterPropTableWriter)
    raise TypeError(f"Invalid loader type: {loader_type}")


def write(report_file: Path):
    """ Helper function to write a table from a report file. """
    # Determine the needed type of report loader.
    report_loader_type = infer_report_loader_type(report_file)
    # Load the report.
    report_loader = report_loader_type.open(report_file)
    # For each table type, create the table, write it, and return the
    # path to the table output file.
    return [table_type(report_loader).write()
            for table_type in get_table_types(report_loader_type)]
