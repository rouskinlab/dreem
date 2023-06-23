from abc import ABC, abstractmethod
from functools import cached_property
from logging import getLogger

import numpy as np
import pandas as pd

from .base import (TOTAL_REL, DELET_REL, INSRT_REL, MATCH_REL, MUTAT_REL,
                   SUBST_REL, SUB_A_REL, SUB_C_REL, SUB_G_REL, SUB_T_REL,
                   POPAVG_TITLE, REL_NAME)
from ..cluster.indexes import CLS_NAME, ORD_NAME, READ_NAME
from ..cluster.load import ClustLoader
from ..core.bitc import BitCaller, SemiBitCaller
from ..core.bitv import BitCounter, ClustBitCounter, CloseEmptyBitAccumError
from ..core.mu import calc_f_obs_df, calc_mu_adj_df
from ..core.sect import Section
from ..mask.load import MaskLoader
from ..relate.load import RelateLoader

logger = getLogger(__name__)


# Tabulator Classes ####################################################

class Tabulator(ABC):
    """ Base class for tabulating data for multiple tables from a report
    loader. """

    def __init__(self, loader: RelateLoader | MaskLoader | ClustLoader):
        self._loader = loader

    @property
    def out_dir(self):
        return self._loader.out_dir

    @property
    def sample(self):
        return self._loader.sample

    @property
    def ref(self):
        return self._loader.ref

    @property
    def sect(self):
        return self._loader.sect

    @property
    def seq(self):
        """ Sequence of the section covered by the table. """
        return self._loader.seq

    @property
    def positions(self):
        """ Numeric positions of the section covered by the table. """
        return self._loader.positions

    @property
    def index(self):
        """ Index (rows) of the table. """
        return self._loader.index

    @property
    @abstractmethod
    def columns(self):
        """ Columns of the table. """

    @cached_property
    def rels(self):
        """ Relationships counted in the tables. """
        return [rel for rel, _, _ in iter_bit_callers()]

    @cached_property
    @abstractmethod
    def bit_counts(self) -> dict[str, BitCounter]:
        """ Return BitCounters for the batches of relation vectors. """

    @classmethod
    @abstractmethod
    def get_null_value(cls) -> int | float:
        """ The null value for a count: either 0 or NaN. """

    @abstractmethod
    def _tabulate_by_pos(self):
        return pd.DataFrame.from_dict({rel: counter.nyes_per_pos
                                       for rel, counter
                                       in self.bit_counts.items()})

    @abstractmethod
    def _tabulate_by_read(self):
        return pd.DataFrame.from_dict({rel: counter.nyes_per_read
                                       for rel, counter
                                       in self.bit_counts.items()})

    def tabulate_by_pos(self):
        """ DataFrame of the count for each position and bit-caller. """
        try:
            return self._tabulate_by_pos()
        except CloseEmptyBitAccumError:
            # No reads were given.
            return pd.DataFrame(self.get_null_value(),
                                index=self._loader.index,
                                columns=self.columns)

    def tabulate_by_read(self):
        """ DataFrame of the count for each read and bit-caller. """
        try:
            return self._tabulate_by_read()
        except CloseEmptyBitAccumError:
            # No reads were given.
            return pd.DataFrame(self.get_null_value(),
                                index=pd.Index(name=self._loader.index.name),
                                columns=self.columns)


class RelTabulator(Tabulator):

    @classmethod
    def get_null_value(cls):
        return 0

    @property
    def columns(self):
        return pd.Index(self.rels, name=REL_NAME)

    @cached_property
    def bit_counts(self):
        return {
            name: BitCounter(bc.iter(self._loader.iter_batches_processed()))
            for name, bc, _ in iter_bit_callers()
        }

    def _tabulate_by_pos(self):
        return super()._tabulate_by_pos()

    def _tabulate_by_read(self):
        return super()._tabulate_by_read()


class MaskTabulator(Tabulator):

    @classmethod
    def get_null_value(cls):
        return np.nan

    @property
    def columns(self):
        return pd.Index(self.rels, name=REL_NAME)

    @cached_property
    def bit_counts(self):
        return {
            rel: BitCounter(self._loader.iter_batches_processed(bit_caller=bc,
                                                                **kwargs))
            for rel, bc, kwargs in iter_bit_callers()
        }

    def _tabulate_by_pos(self):
        # Count every type of relationship at each position, in the same
        # way as for the superclass.
        counts_obs = super()._tabulate_by_pos()
        # Adjust the counts to correct for observer bias.
        return adjust_counts(counts_obs,
                             self._loader.section,
                             self._loader.min_mut_gap)

    def _tabulate_by_read(self):
        return super()._tabulate_by_read()


class ClustTabulator(Tabulator):

    @classmethod
    def get_null_value(cls):
        return np.nan

    @property
    def clusters(self):
        """ Order and number of each cluster. """
        return self._loader.clusters

    def get_read_names(self):
        return self._loader.import_loader.get_read_names()

    @cached_property
    def columns(self):
        return pd.MultiIndex.from_tuples([(order, cluster, rel)
                                          for order, cluster in self.clusters
                                          for rel in self.rels],
                                         names=[ORD_NAME, CLS_NAME, REL_NAME])

    @cached_property
    def bit_counts(self):
        return dict(
            (rel,
             ClustBitCounter(self._loader.iter_batches_processed(bit_caller=bc,
                                                                 **kwargs)))
            for rel, bc, kwargs in iter_bit_callers()
        )

    def _tabulate_by_pos(self):
        """ DataFrame of the bit count for each position and caller. """
        # Initialize an empty DataFrame.
        counts_obs = pd.DataFrame(self.get_null_value(),
                                  index=self._loader.index,
                                  columns=self.columns)
        # Fill in the DataFrame column by column.
        for col in self.columns:
            nk, k, rel = col
            counts_obs[col] = self.bit_counts[rel].nyes_per_pos.loc[:, (nk, k)]
        # Adjust the counts to correct for observer bias.
        counts_adj = dict()
        for order, k in self.clusters:
            for rel, adj in adjust_counts(counts_obs.loc[:, (order, k)],
                                          self._loader.section,
                                          self._loader.min_mut_gap).items():
                counts_adj[order, k, rel] = adj
        counts_adj = pd.DataFrame.from_dict(counts_adj)
        counts_adj.columns.rename(counts_obs.columns.names, inplace=True)
        return counts_adj

    def _tabulate_by_read(self):
        """ DataFrame of the bit count for each read and caller. """
        # Initialize an empty DataFrame.
        counts = pd.DataFrame(self.get_null_value(),
                              index=pd.Index(self.get_read_names(),
                                             name=READ_NAME),
                              columns=self.columns)
        # Fill in the DataFrame column by column.
        for col in self.columns:
            order, k, rel = col
            counts[col] = self.bit_counts[rel].nyes_per_read.loc[:, (order, k)]
        return counts

    def tabulate_by_clust(self):
        """ Proportion of each cluster. """
        return self._loader.props.to_frame()


# Helper functions #####################################################

def iter_bit_callers():
    """ Yield a BitCaller for each type of relationship to tabulate. """
    # Call reference matches.
    refc = SemiBitCaller.from_counts(count_ref=True)
    # Call mutations.
    mutc = SemiBitCaller.from_counts(count_sub=True,
                                     count_del=True,
                                     count_ins=True)
    # Create a standard bit caller for all matches and mutations.
    bitc = BitCaller(refc, mutc)
    # Count all base calls (everything but the bytes 0 and 255).
    yield TOTAL_REL, bitc, dict(merge=True)
    # Count matches to the reference sequence.
    yield MATCH_REL, bitc, dict(invert=True)
    # Count all types of mutations, relative to reference matches.
    yield MUTAT_REL, bitc, dict()
    # Count each type of mutation, relative to reference matches.
    yield (SUBST_REL,
           BitCaller(refc, SemiBitCaller.from_counts(count_sub=True)), dict())
    yield (SUB_A_REL,
           BitCaller(refc, SemiBitCaller("ca", "ga", "ta")), dict())
    yield (SUB_C_REL,
           BitCaller(refc, SemiBitCaller("ac", "gc", "tc")), dict())
    yield (SUB_G_REL,
           BitCaller(refc, SemiBitCaller("ag", "cg", "tg")), dict())
    yield (SUB_T_REL,
           BitCaller(refc, SemiBitCaller("at", "ct", "gt")), dict())
    yield (DELET_REL,
           BitCaller(refc, SemiBitCaller.from_counts(count_del=True)), dict())
    yield (INSRT_REL,
           BitCaller(refc, SemiBitCaller.from_counts(count_ins=True)), dict())


def adjust_counts(counts_obs: pd.DataFrame,
                  section: Section,
                  min_mut_gap: int):
    """
    Adjust the given table of masked/clustered bit counts per position
    to correct for observer bias. The table is mutated in-place.

    Parameters
    ----------
    counts_obs: DataFrame
        Counts of the bits for each type of relation (column) at each
        position (index) in the section of interest.
    section: Section
        The section of interest.
    min_mut_gap: int
        Minimum number of non-mutated bases permitted between mutations.
    """
    logger.debug(f"Adjusting mutation counts of table\n{counts_obs}")
    # Compute the observed fraction of mutations at each position.
    nrefs_obs = counts_obs.loc[:, MATCH_REL]
    nmuts_obs = counts_obs.loc[:, MUTAT_REL]
    ninfo_obs = nrefs_obs + nmuts_obs
    with np.errstate(divide="ignore"):
        # Ignore division by zero, which is acceptable here because any
        # NaN values will be zeroed out subsequently.
        fmuts_obs = nmuts_obs / ninfo_obs
    # Fill missing NaN values with zero.
    fmuts_obs.fillna(0., inplace=True)
    # Adjust the fraction of mutations to correct the observer bias.
    mu_adj = calc_mu_adj_df(fmuts_obs.to_frame(),
                            section, min_mut_gap).squeeze(axis=1)
    # Compute the fraction of reads that would be observed.
    f_obs = float(calc_f_obs_df(mu_adj.to_frame(), section, min_mut_gap)[0])
    # Initialize new data for the adjusted counts.
    counts_adj = dict()
    # Assume that the total number of base calls including unobserved
    # reads is the number observed divided by the fraction of the total
    # reads that were observed.
    counts_adj[TOTAL_REL] = counts_obs.loc[:, TOTAL_REL] / f_obs
    # Two conditions must be met:
    # : The number of informative bases (matches + mutations) after
    #   adjustment must equal the number before adjustment divided
    #   by the fraction of reads that were observed:
    #   (nrefs_adj + nmuts_adj) = (nrefs_obs + nmuts_obs) / f_obs
    # : The fraction of mutations at each position after adjustment
    #   must equal the adjusted number of mutations divided by the
    #   adjusted number of informative bases:
    #   nmuts_adj / (nrefs_adj + nmuts_adj) = mu_adj
    # The two unknown variables (nrefs_adj and nmuts_adj) can be
    # solved using the above system of two equations. By solving for
    # (nrefs_adj + nmuts_adj) in both equations, we get:
    # : (nrefs_adj + nmuts_adj) = (nrefs_obs + nmuts_obs) / f_obs
    # : (nrefs_adj + nmuts_adj) = nmuts_adj / mu_adj, if mu_adj > 0
    # Setting both right hand sides equal, nmuts_adj is solved:
    # : (nrefs_obs + nmuts_obs) / f_obs = nmuts_adj / mu_adj
    # : nmuts_adj = (nrefs_obs + nmuts_obs) * (mu_adj / f_obs)
    nmuts_adj = ninfo_obs * (mu_adj / f_obs)
    # Then, plugging this solution into the second equation:
    # : nrefs_adj = nmuts_adj / mu_adj - nmuts_adj
    #             = nmuts_adj * (1 / mu_adj - 1), if mu_adj > 0
    # If, on the other hand, mu_adj = 0, then
    # : nmuts_adj = mu_adj * (nrefs_adj + nmuts_adj) = 0
    # Thus, we solve for nrefs_adj using the first equation instead:
    # : (nrefs_adj + nmuts_adj) = (nrefs_obs + nmuts_obs) / f_obs
    # : nrefs_adj = nrefs_obs / f_obs
    with np.errstate(divide="ignore"):
        # Ignore division-by-zero warnings in np.reciprocal because
        # positions where mu_adj = 0 are just skipped by np.where.
        nrefs_adj = np.where(mu_adj > 0.,
                             nmuts_adj * (np.reciprocal(mu_adj) - 1.),
                             nrefs_obs / f_obs)
    counts_adj[MATCH_REL] = nrefs_adj
    counts_adj[MUTAT_REL] = nmuts_adj
    # Compute the factor by which nmuts was adjusted:
    nmuts_fac = nmuts_adj / nmuts_obs
    # Adjust every type of mutation by this factor.
    for mut in (SUB_A_REL, SUB_C_REL, SUB_G_REL, SUB_T_REL,
                SUBST_REL, DELET_REL, INSRT_REL):
        counts_adj[mut] = counts_obs.loc[:, mut] * nmuts_fac
    # Assemble the adjusted counts into a new DataFrame.
    return pd.DataFrame.from_dict(counts_adj)


def tabulate_loader(report_loader: RelateLoader | MaskLoader | ClustLoader):
    """ Return a new DataLoader, choosing the subclass based on the type
    of the argument `report_loader`. """
    if isinstance(report_loader, RelateLoader):
        return RelTabulator(report_loader)
    if isinstance(report_loader, MaskLoader):
        return MaskTabulator(report_loader)
    if isinstance(report_loader, ClustLoader):
        return ClustTabulator(report_loader)
    raise TypeError(f"Invalid loader type: {type(report_loader).__name__}")
