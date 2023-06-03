from abc import ABC, abstractmethod
from functools import cached_property
from pathlib import Path
from typing import Any

import pandas as pd

from ..core import path
from ..core.seq import DNA


# General fields
POS_FIELD = "Position"
SEQ_FIELD = "Base"
READ_FIELD = "Read Name"
POPAVG_TITLE = "Pop Avg"

# Count fields
TOTAL_FIELD = "Called"
MATCH_FIELD = "Matched"
MUTAT_FIELD = "Mutated"
DELET_FIELD = "Deleted"
INSRT_FIELD = "Inserted"
SUB_N_FIELD = "Subbed"
SUB_A_FIELD = "Sub-A"
SUB_C_FIELD = "Sub-C"
SUB_G_FIELD = "Sub-G"
SUB_T_FIELD = "Sub-T"


# Table Base Classes ###################################################

class Table(ABC):
    """ Table base class. """

    @property
    @abstractmethod
    def out_dir(self):
        """ Path of the table's output directory. """
        return Path()

    @property
    @abstractmethod
    def sample(self):
        """ Name of the table's sample. """
        return ""

    @property
    @abstractmethod
    def ref(self):
        """ Name of the table's reference. """
        return ""

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
    @abstractmethod
    def by_read(cls):
        """ Whether the table contains data for each read. """
        return False

    @classmethod
    @abstractmethod
    def path_segs(cls) -> list[path.Segment]:
        """ Table's path segments. """
        return list()

    @classmethod
    def gzipped(cls):
        """ Whether the table's file is compressed with gzip. """
        return cls.by_read()

    @classmethod
    def ext(cls):
        """ Table's file extension: either '.csv' or '.csv.gz'. """
        return path.CSVZIP_EXT if cls.gzipped() else path.CSV_EXT

    @property
    @abstractmethod
    def path_fields(self) -> dict[str, Any]:
        """ Table's path fields. """
        return {path.TOP: self.out_dir,
                path.MOD: path.MOD_TABLE,
                path.SAMP: self.sample,
                path.REF: self.ref,
                path.TABLE: self.kind(),
                path.EXT: self.ext()}

    @cached_property
    def path(self):
        """ Path of the table's CSV file (possibly gzipped). """
        return path.buildpar(*self.path_segs(), **self.path_fields)


class CountTable(Table, ABC):
    """ Table whose data are counts. """

    @property
    def nb(self):
        """ Number of total base calls. """
        return self.data[TOTAL_FIELD]

    @property
    def nr(self):
        """ Number of reference matches. """
        return self.data[MATCH_FIELD]

    @property
    def nm(self):
        """ Number of mutations. """
        return self.data[MUTAT_FIELD]

    @property
    def nd(self):
        """ Number of deletions. """
        return self.data[DELET_FIELD]

    @property
    def ni(self):
        """ Number of insertions. """
        return self.data[INSRT_FIELD]

    @property
    def ns(self):
        """ Number of substitutions to any nucleotide. """
        return self.data[SUB_N_FIELD]

    @property
    def na(self):
        """ Number of substitutions to A. """
        return self.data[SUB_A_FIELD]

    @property
    def nc(self):
        """ Number of substitutions to C. """
        return self.data[SUB_C_FIELD]

    @property
    def ng(self):
        """ Number of substitutions to G. """
        return self.data[SUB_G_FIELD]

    @property
    def nt(self):
        """ Number of substitutions to T. """
        return self.data[SUB_T_FIELD]

    @cached_property
    def nn(self):
        """ Number of informative mutated or matching bases. """
        return self.nr + self.nm

    @cached_property
    def fn(self):
        """ Fraction of informative (unambiguous) bases. """
        return self.nn / self.nb

    @cached_property
    def fr(self):
        """ Fraction of reference matches. """
        return self.nr / self.nn

    @cached_property
    def fm(self):
        """ Fraction of mutations. """
        return self.nm / self.nn

    @cached_property
    def fd(self):
        """ Fraction of deletions. """
        return self.nd / self.nn

    @cached_property
    def fi(self):
        """ Fraction of insertions. """
        return self.ni / self.nn

    @cached_property
    def fs(self):
        """ Fraction of substitutions to any nucleotide. """
        return self.ns / self.nn

    @cached_property
    def fa(self):
        """ Fraction of substitutions to A. """
        return self.na / self.nn

    @cached_property
    def fc(self):
        """ Fraction of substitutions to C. """
        return self.nc / self.nn

    @cached_property
    def fg(self):
        """ Fraction of substitutions to G. """
        return self.ng / self.nn

    @cached_property
    def ft(self):
        """ Fraction of substitutions to T. """
        return self.nt / self.nn


class SectTable(Table, ABC):
    """ Table with a section of a reference sequence. """

    @property
    @abstractmethod
    def sect(self):
        """ Name of the table's section. """
        return ""

    @classmethod
    def path_segs(cls):
        return [path.ModSeg, path.SampSeg, path.RefSeg, path.SectSeg,
                path.MutTabSeg]

    @property
    def path_fields(self):
        return {**super().path_fields, path.SECT: self.sect}


# Table by Source (relate/mask/cluster) ################################

class RelTable(CountTable, ABC):
    """ Table of relation vectors. """

    @classmethod
    def path_segs(cls):
        return [path.ModSeg, path.SampSeg, path.RefSeg, path.MutTabSeg]

    @property
    def path_fields(self) -> dict[str, Any]:
        return super().path_fields


class MaskTable(SectTable, CountTable, ABC):
    """ Table of masked bit vectors. """


class ClustTable(SectTable, ABC):
    """ Table of clustering results. """


# Table by Index (position/read/proportion) ############################

class PosTable(Table, ABC):
    """ Table indexed by position. """

    @classmethod
    def by_read(cls):
        return False

    @property
    def positions(self):
        return self.data.index.values

    @property
    def end5(self):
        return int(self.positions[0])

    @property
    def end3(self):
        return int(self.positions[-1])

    @property
    def bases(self):
        return self.data[SEQ_FIELD]

    @cached_property
    def seq(self):
        return DNA("".join(self.bases).encode())


class ReadTable(Table, ABC):
    """ Table indexed by read. """

    @classmethod
    def by_read(cls):
        return True

    @property
    def reads(self):
        return self.data.index.values


class PropTable(Table, ABC):
    """ Table indexed by cluster. """

    @classmethod
    def by_read(cls):
        return False

    @property
    def clusters(self):
        return self.data.index.values


# Table by Source and Index ############################################

class RelPosTable(RelTable, PosTable, ABC):
    @classmethod
    def kind(cls):
        return path.RELATE_POS_TAB


class RelReadTable(RelTable, ReadTable, ABC):
    @classmethod
    def kind(cls):
        return path.RELATE_READ_TAB


class MaskPosTable(MaskTable, PosTable, ABC):
    @classmethod
    def kind(cls):
        return path.MASKED_POS_TAB


class MaskReadTable(MaskTable, ReadTable, ABC):
    @classmethod
    def kind(cls):
        return path.MASKED_READ_TAB


class ClustPosTable(ClustTable, PosTable, ABC):
    @classmethod
    def kind(cls):
        return path.CLUST_MUS_TAB

    @cached_property
    def cluster_names(self):
        return self.data.columns.drop(SEQ_FIELD).to_list()


class ClustReadTable(ClustTable, ReadTable, ABC):
    @classmethod
    def kind(cls):
        return path.CLUST_RESP_TAB

    @cached_property
    def cluster_names(self):
        return self.data.columns.to_list()


class ClustPropTable(ClustTable, PropTable, ABC):
    @classmethod
    def kind(cls):
        return path.CLUST_PROP_TAB

    @cached_property
    def cluster_names(self):
        return self.clusters.to_list()
