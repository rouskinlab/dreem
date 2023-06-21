from abc import ABC, abstractmethod
from functools import cached_property
from logging import getLogger

import numpy as np
import pandas as pd

from .base import (DELET_FIELD, INSRT_FIELD, MATCH_FIELD, MUTAT_FIELD,
                   SUB_A_FIELD, SUB_C_FIELD, SUB_G_FIELD, SUB_T_FIELD,
                   SUBST_FIELD, TOTAL_FIELD, POPAVG_TITLE)
from ..cluster.load import ClustLoader
from ..core.bitc import BitCaller, SemiBitCaller
from ..core.bitv import BitCounter, CloseEmptyBitAccumError
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
    @abstractmethod
    def eligible_bits(self):
        """ Return the bits that are eligible to be counted. """
        return BitCaller(SemiBitCaller(), SemiBitCaller())

    def iter_bit_callers(self):
        """ Yield every BitCaller, one for each type of information to
        include in the table of counts. """
        # Initialize two semi-bit-callers for convenience: one to call
        # reference matches, the other to call mutations.
        inter = SemiBitCaller.inter
        union = SemiBitCaller.union
        nosc = self.eligible_bits.nos_call
        yesc = self.eligible_bits.yes_call
        refc = inter(SemiBitCaller.from_counts(count_ref=True), nosc,
                     cache_all=True)
        mutc = inter(SemiBitCaller.from_counts(count_sub=True,
                                               count_del=True,
                                               count_ins=True), yesc,
                     cache_all=True)
        # Count all base calls (everything but the bytes 0 and 255, and
        # any types of base calls excluded by the mask).
        yield (TOTAL_FIELD,
               BitCaller(SemiBitCaller(),
                         inter(SemiBitCaller.from_counts(count_ref=True,
                                                         count_sub=True,
                                                         count_del=True,
                                                         count_ins=True),
                               union(nosc, yesc))))
        # Count matches to the reference sequence.
        yield MATCH_FIELD, BitCaller(mutc, refc)
        # Count all types of mutations, relative to reference matches.
        yield MUTAT_FIELD, BitCaller(refc, mutc)
        # Count each type of mutation, relative to reference matches.
        yield (SUBST_FIELD,
               BitCaller(refc, inter(SemiBitCaller.from_counts(count_sub=True),
                                     yesc)))
        yield (SUB_A_FIELD,
               BitCaller(refc, inter(SemiBitCaller("ca", "ga", "ta"), yesc)))
        yield (SUB_C_FIELD,
               BitCaller(refc, inter(SemiBitCaller("ac", "gc", "tc"), yesc)))
        yield (SUB_G_FIELD,
               BitCaller(refc, inter(SemiBitCaller("ag", "cg", "tg"), yesc)))
        yield (SUB_T_FIELD,
               BitCaller(refc, inter(SemiBitCaller("at", "ct", "gt"), yesc)))
        yield (DELET_FIELD,
               BitCaller(refc, inter(SemiBitCaller.from_counts(count_del=True),
                                     yesc)))
        yield (INSRT_FIELD,
               BitCaller(refc, inter(SemiBitCaller.from_counts(count_ins=True),
                                     yesc)))

    @cached_property
    def fields(self):
        """ Fields of the count tables. """
        return [field for field, _ in self.iter_bit_callers()]

    @classmethod
    def get_bit_counter_type(cls):
        """ Type of the BitCounter. """
        return BitCounter

    @abstractmethod
    def _accumulate_batches(self, callers: dict[str, BitCaller],
                            counters: dict[str, BitCounter]):
        """ Accumulate the counts from the batches. """

    @cached_property
    def bit_counters(self) -> dict[str, BitCounter]:
        """ Return BitCounters for the batches of relation vectors. """
        # Collect all the bit callers.
        callers = dict(self.iter_bit_callers())
        # Create a bit counter for each bit caller.
        counters = {field: self.get_bit_counter_type()() for field in callers}
        # Accumulate the bit counts over all batches.
        self._accumulate_batches(callers, counters)
        # Prevent each bit counter from accepting more data.
        for counter in counters.values():
            counter.close()
        return counters

    @classmethod
    @abstractmethod
    def get_null_value(cls) -> int | float:
        """ The null value for a count: either 0 or NaN. """

    def get_null_counts(self):
        """ Default null data. """
        return pd.DataFrame(self.get_null_value(),
                            index=self._loader.index,
                            columns=self.columns)

    def tabulate_by_pos(self):
        """ DataFrame of the bit count for each position and caller. """
        try:
            return pd.DataFrame.from_dict({field: counter.nyes_per_pos
                                           for field, counter
                                           in self.bit_counters.items()})
        except CloseEmptyBitAccumError:
            return self.get_null_counts()

    def tabulate_by_read(self):
        """ DataFrame of the bit count for each read and caller. """
        try:
            return pd.DataFrame.from_dict({field: counter.nyes_per_read
                                           for field, counter
                                           in self.bit_counters.items()})
        except CloseEmptyBitAccumError:
            # No reads were given.
            return pd.DataFrame(0, index=[], columns=self.fields)


class RelTabulator(Tabulator):

    @classmethod
    def get_null_value(cls):
        return 0

    @property
    def columns(self):
        return self.fields

    @cached_property
    def eligible_bits(self):
        # All bits are eligible.
        all_bits = SemiBitCaller.from_counts(count_ref=True, count_sub=True,
                                             count_del=True, count_ins=True)
        return BitCaller(all_bits, all_bits)

    def _accumulate_batches(self, callers: dict[str, BitCaller],
                            counters: dict[str, BitCounter]):
        # Iterate over all batches.
        for batch in self._loader.iter_batches_public():
            # Use each bit caller to create bit vectors from this batch.
            for field, counter in counters.items():
                # Add the bit vectors to the appropriate counter.
                counter.add_batch(callers[field].call(batch))


class MaskTabulator(Tabulator):

    @classmethod
    def get_null_value(cls):
        return np.nan

    @property
    def columns(self):
        return self.fields

    @cached_property
    def eligible_bits(self):
        # Only the bits specified by the mask are eligible.
        return self._loader.bit_caller

    def _accumulate_batches(self, callers: dict[str, BitCaller],
                            counters: dict[str, BitCounter]):
        # Iterate over all batches.
        for batch_mods in self._loader.iter_batches_modified(callers.values()):
            # Use each bit caller to create bit vectors from this batch.
            for batch_mod, counter in zip(batch_mods, counters.values(),
                                          strict=True):
                # Add the bit vectors to the appropriate counter.
                counter.add_batch(batch_mod)

    def tabulate_by_pos(self):
        # Count every type of relationship at each position, in the same
        # way as for the superclass.
        nrels_per_pos = super().tabulate_by_pos()
        # Adjust the counts to correct for observer bias.
        adjust_counts(nrels_per_pos,
                      self._loader.section,
                      self._loader.min_mut_gap)
        return nrels_per_pos


class ClustTabulator(Tabulator):

    @classmethod
    def get_null_value(cls):
        return np.nan

    @property
    def clusters(self):
        """ Order and number of each cluster. """
        return self._loader.clusters

    @cached_property
    def columns(self):
        return pd.MultiIndex((order, cluster, field)
                             for order, cluster in self.clusters
                             for field in self.fields)

    def tabulate_by_clust(self):
        """ Proportion of each cluster. """
        return self._loader.props.to_frame()


# Helper functions #####################################################

def adjust_counts(nrels_per_pos: pd.DataFrame,
                  section: Section,
                  min_mut_gap: int):
    """
    Adjust the given table of masked/clustered bit counts per position
    to correct for observer bias. The table is mutated in-place.

    Parameters
    ----------
    nrels_per_pos: DataFrame
        Counts of the bits for each type of relation (column) at each
        position (index) in the section of interest. Mutated in-place.
    section: Section
        The section of interest.
    min_mut_gap: int
        Minimum number of non-mutated bases permitted between mutations.
    """
    logger.debug(f"Adjusting mutation counts of table\n{nrels_per_pos}")
    # Compute the observed fraction of mutations at each position.
    nrefs_obs = nrels_per_pos[MATCH_FIELD]
    nmuts_obs = nrels_per_pos[MUTAT_FIELD]
    ninfo_obs = nrefs_obs + nmuts_obs
    fmuts_obs = nmuts_obs / ninfo_obs
    # Adjust the fraction of mutations to correct the observer bias.
    mu_adj = calc_mu_adj_df(fmuts_obs.to_frame(POPAVG_TITLE),
                            section, min_mut_gap)
    # Compute the fraction of reads that would be observed.
    f_obs = calc_f_obs_df(mu_adj, section, min_mut_gap)[POPAVG_TITLE]
    # Assume that the total number of base calls including unobserved
    # reads is the number observed divided by the fraction of the total
    # reads that were observed.
    nrels_per_pos[TOTAL_FIELD] /= f_obs
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
    # : (nrefs_adj + nmuts_adj) = nmuts_adj / mu_adj
    # Setting both right hand sides equal, nmuts_adj is solved:
    # : (nrefs_obs + nmuts_obs) / f_obs = nmuts_adj / mu_adj
    # : nmuts_adj = (nrefs_obs + nmuts_obs) * (mu_adj / f_obs)
    # Then, plugging this solution into the second equation:
    # : nrefs_adj = nmuts_adj / mu_adj - nmuts_adj
    #             = nmuts_adj * (1 / mu_adj - 1)
    nmuts_adj = ninfo_obs * (mu_adj[POPAVG_TITLE] / f_obs)
    nrefs_adj = nmuts_adj * (np.reciprocal(mu_adj[POPAVG_TITLE]) - 1.)
    nrels_per_pos[MATCH_FIELD] = nrefs_adj
    # Compute the factor by which nmuts was adjusted:
    nmuts_fac = nmuts_adj / nmuts_obs
    # Adjust every type of mutation by this factor.
    nrels_per_pos[MUTAT_FIELD] = nmuts_adj
    for mut in (SUB_A_FIELD, SUB_C_FIELD, SUB_G_FIELD, SUB_T_FIELD,
                SUBST_FIELD, DELET_FIELD, INSRT_FIELD):
        nrels_per_pos[mut] *= nmuts_fac


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
