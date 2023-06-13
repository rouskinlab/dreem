from abc import ABC, abstractmethod
from functools import cache, cached_property
from pathlib import Path
from typing import Any

import pandas as pd

from ..cluster.load import parse_names
from ..core import path
from ..core.seq import DNA

# General fields
POS_FIELD = "Position"
SEQ_FIELD = "Base"
READ_FIELD = "Read Name"
POPAVG_TITLE = "Pop Avg"

# Count fields
TOTAL_FIELD = "Called"
INFOR_FIELD = "Informed"
MATCH_FIELD = "Matched"
MUTAT_FIELD = "Mutated"
DELET_FIELD = "Deleted"
INSRT_FIELD = "Inserted"
SUBST_FIELD = "Subbed"
SUB_A_FIELD = "Subbed-to-A"
SUB_C_FIELD = "Subbed-to-C"
SUB_G_FIELD = "Subbed-to-G"
SUB_T_FIELD = "Subbed-to-T"


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

    FIELD_CODES = {
        'b': TOTAL_FIELD,
        'n': INFOR_FIELD,
        'r': MATCH_FIELD,
        'm': MUTAT_FIELD,
        'd': DELET_FIELD,
        'i': INSRT_FIELD,
        's': SUBST_FIELD,
        'a': SUB_A_FIELD,
        'c': SUB_C_FIELD,
        'g': SUB_G_FIELD,
        't': SUB_T_FIELD,
    }

    @cache
    def _get_field_count(self, field: str):
        # Pull the data field from the table's data frame.
        if field == INFOR_FIELD:
            count = self.data[MATCH_FIELD] + self.data[MUTAT_FIELD]
        else:
            count = self.data[field].copy()
        # Name the series after the field.
        count.rename(field, inplace=True)
        return count

    @cache
    def _get_field_frac(self, field: str):
        # Determine the denominator of the fraction field.
        denom = TOTAL_FIELD if field == INFOR_FIELD else INFOR_FIELD
        # Compute the ratio of the field and its denominator.
        frac = self._get_field_count(field) / self._get_field_count(denom)
        # Name the series after the field.
        frac.rename(field, inplace=True)
        return frac

    def get_field_count(self, code: str):
        """ Count the bits for a field given its field code. """
        return self._get_field_count(self.FIELD_CODES[code])

    def get_field_frac(self, code: str):
        """ Get the fraction for a field given its field code. """
        return self._get_field_frac(self.FIELD_CODES[code])


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

    @property
    def ks(self) -> list[int]:
        """ Numbers of clusters. """
        return sorted({k for k, c in parse_names(self.cluster_names)})


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
