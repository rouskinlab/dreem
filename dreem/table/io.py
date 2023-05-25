from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path
from typing import Any

import pandas as pd

from ..call.load import BitVecLoader
from ..cluster.load import ClusterLoader
from ..core import path
from ..core.bit import SemiBitCaller, BitCaller, BitCounter
from ..core.mu import calc_mus
from ..core.seq import DNA
from ..relate.load import RelVecLoader

# Field Definitions ####################################################

POS_FIELD = "Position"
SEQ_FIELD = "Base"
READ_FIELD = "Read Name"
NINFO_FIELD = "Num Info"
NMUTS_FIELD = "Num Muts"
FINFO_FIELD = "Frac Info"
FMUTS_FIELD = "Frac Muts"
FMUTO_FIELD = "Frac Muts (Obs)"
FMUTA_FIELD = "Frac Muts (Adj)"
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
        return cls.kind() in (path.BITVEC_POS_TAB,
                              path.BITVEC_READ_TAB,
                              path.CLUST_MUS_TAB,
                              path.CLUST_PROP_TAB,
                              path.CLUST_RESP_TAB)

    @classmethod
    def by_read(cls):
        """ Whether the table contains data for each read. """
        return cls.kind() in (path.RELVEC_READ_TAB,
                              path.BITVEC_READ_TAB)

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
        return cls.kind() in (path.RELVEC_READ_TAB,
                              path.BITVEC_READ_TAB,
                              path.CLUST_RESP_TAB)

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
    def fmuts_adj(self):
        return self.data[FMUTA_FIELD]


class BitVecReadTable(BitVecTable, ReadTable, ABC):
    @classmethod
    def kind(cls):
        return path.BITVEC_READ_TAB

    @property
    def fmuts(self):
        return self.data[FMUTS_FIELD]


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
            FMUTA_FIELD: calc_mus(counter.fmuts_per_pos.to_frame(),
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


class BitVecPosTableLoader(BitVecPosTable, TableLoader):
    pass


class BitVecReadTableLoader(BitVecReadTable, TableLoader):
    pass


# Helper Functions #####################################################

def infer_loader_type(report_file: Path):
    """ Given a report file path, infer the type of Loader it needs. """
    if path.RelateRepSeg.ptrn.match(report_file.name):
        return RelVecLoader
    if path.CallRepSeg.ptrn.match(report_file.name):
        return BitVecLoader
    if path.ClustRepSeg.ptrn.match(report_file.name):
        return ClusterLoader
    raise ValueError(f"Failed to infer loader for {report_file}")


def infer_loader(report_file: Path):
    """ Given a report file path, return a loader for the file. """
    return infer_loader_type(report_file).open(report_file)


def write(report_file: Path, table_type: type[TableWriter]):
    """ Helper function to write a table from a report file. """
    table = table_type(infer_loader(report_file))
    table.write()
    return table.path
