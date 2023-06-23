from abc import ABC, abstractmethod
from functools import cache, cached_property
from pathlib import Path
from typing import Any

import pandas as pd

from ..core import path
from ..core.seq import DNA

# General fields
POS_TITLE = "Position"
SEQ_TITLE = "Base"
READ_TITLE = "Read Name"
REL_NAME = "Relationship"
POPAVG_TITLE = "pop-avg"

# Count relationships
TOTAL_REL = "Called"
INFOR_REL = "Informed"
MATCH_REL = "Matched"
MUTAT_REL = "Mutated"
DELET_REL = "Deleted"
INSRT_REL = "Inserted"
SUBST_REL = "Subbed"
SUB_A_REL = "Subbed-A"
SUB_C_REL = "Subbed-C"
SUB_G_REL = "Subbed-G"
SUB_T_REL = "Subbed-T"


# Table Base Classes ###################################################

class Table(ABC):
    """ Table base class. """

    REL_CODES = {
        'b': TOTAL_REL,
        'n': INFOR_REL,
        'r': MATCH_REL,
        'm': MUTAT_REL,
        'd': DELET_REL,
        'i': INSRT_REL,
        's': SUBST_REL,
        'a': SUB_A_REL,
        'c': SUB_C_REL,
        'g': SUB_G_REL,
        't': SUB_T_REL,
    }

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

    @cache
    def _get_rel_count(self, rel: str):
        # Pull the relationship from the table's data frame.
        if rel == INFOR_REL:
            count = self.data[MATCH_REL] + self.data[MUTAT_REL]
        else:
            count = self.data[rel].copy()
        # Name the series after the relationship.
        count.rename(rel, inplace=True)
        return count

    @cache
    def _get_rel_frac(self, rel: str):
        # Determine the relationship to use as the denominator.
        denom = TOTAL_REL if rel == INFOR_REL else INFOR_REL
        # Compute the ratio of the numerator and the denominator.
        frac = self._get_rel_count(rel) / self._get_rel_count(denom)
        # Name the series after the relationship.
        frac.rename(rel, inplace=True)
        return frac

    def get_rel_count(self, code: str):
        """ Count the bits for a relationship given its code. """
        return self._get_rel_count(self.REL_CODES[code])

    def get_rel_frac(self, code: str):
        """ Compute the fraction for a relationship given its code. """
        return self._get_rel_frac(self.REL_CODES[code])


class SectTable(Table, ABC):
    """ Table with a section of a reference sequence. """

    @property
    @abstractmethod
    def sect(self):
        """ Name of the table's section. """
        return ""

    @classmethod
    def path_segs(cls):
        return (path.ModSeg, path.SampSeg, path.RefSeg, path.SectSeg,
                path.MutTabSeg)

    @property
    def path_fields(self):
        return {**super().path_fields, path.SECT: self.sect}


# Table by Source (relate/mask/cluster) ################################

class RelTable(Table, ABC):
    """ Table of relation vectors. """

    @classmethod
    def path_segs(cls):
        return path.ModSeg, path.SampSeg, path.RefSeg, path.MutTabSeg

    @property
    def path_fields(self) -> dict[str, Any]:
        return super().path_fields


class MaskTable(SectTable, ABC):
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
        return self.data[SEQ_TITLE]

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
        return path.CLUST_POS_TAB

    @cached_property
    def cluster_names(self):
        return self.data.columns.drop(SEQ_TITLE).to_list()


class ClustReadTable(ClustTable, ReadTable, ABC):
    @classmethod
    def kind(cls):
        return path.CLUST_READ_TAB

    @cached_property
    def cluster_names(self):
        return self.data.columns.to_list()


class ClustPropTable(ClustTable, PropTable, ABC):
    @classmethod
    def kind(cls):
        return path.CLUST_PROPS_TAB

    @cached_property
    def cluster_names(self):
        return self.clusters.to_list()
