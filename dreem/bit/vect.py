from functools import cache, cached_property
from typing import Any, Iterable

import numpy as np
import pandas as pd

from .call import BitCaller
from .filt import VectorFilter
from ..mut.load import VectorLoader
from ..util.sect import Section, mask_gu, mask_polya, mask_pos


class BitVector(object):
    """ Compute bit vectors from mutation vectors. """

    def __init__(self, /,
                 loader: VectorLoader,
                 section: Section | None = None, *,
                 count_del: bool,
                 count_ins: bool,
                 exclude_polya: int,
                 exclude_gu: bool,
                 exclude_pos: Iterable[int] = (),
                 filter_vec: VectorFilter | None = None):
        """
        Parameters
        ----------
        loader: VectorLoader
            Mutation vector loader
        section: Section | None = None
            Section of the reference sequence to use. If omitted, use
            the entire sequence.
        count_del: bool
            Whether to count deletions as mutations.
        count_ins: bool
            Whether to count insertions as mutations.
        exclude_polya: int
            Exclude stretches of consecutive A bases at least this long.
            If 0, exclude no bases. Must be ≥ 0.
        exclude_gu: bool
            Whether to exclude G and U bases.
        exclude_pos: Iterable[int] = ()
            Additional, arbitrary positions to exclude.
        filter_vec: VectorFilter | None = None
            Filter out low-quality or uninformative reads and positions.
        """
        self.loader = loader
        self.count_del = count_del
        self.count_ins = count_ins
        self.exclude_polya = exclude_polya
        self.exclude_gu = exclude_gu
        self.exclude_pos = list(exclude_pos)
        self.section: Section = loader.section() if section is None else section
        full_mask = np.zeros_like(self.section.positions, dtype=bool)
        # Exclude poly(A) sequences from the section.
        mask = np.logical_and(mask_polya(self.section.seq,
                                         self.exclude_polya),
                              np.logical_not(full_mask))
        self.pos_polya: list[int] = self.section.positions[mask].tolist()
        full_mask = np.logical_or(mask, full_mask)
        # Exclude Gs and Us from the section, if exclude_gu is True.
        mask = np.logical_and(mask_gu(self.section.seq,
                                      self.exclude_gu),
                              np.logical_not(full_mask))
        self.pos_gu: list[int] = self.section.positions[mask].tolist()
        full_mask = np.logical_or(mask, full_mask)
        # Exclude arbitrary, user-specified positions from the section.
        mask = np.logical_and(mask_pos(self.section.positions,
                                       self.exclude_pos),
                              np.logical_not(full_mask))
        self.pos_user: list[int] = self.section.positions[mask].tolist()
        full_mask = np.logical_or(mask, full_mask)
        # Determine which positions remain.
        pos_rem: np.ndarray = self.section.positions[np.logical_not(full_mask)]
        if filter_vec is None:
            # Use all the remaining positions and all the reads.
            self.positions = pos_rem
            self.read_names = None
            # There are no filter parameters.
            self.filter_summary = VectorFilter.null_summary()
        else:
            # Determine the reads and positions to filter.
            for batch in loader.iter_batches(pos_rem):
                # Compute the mutated and matched bits in the batch and
                # feed them through the filter.
                filter_vec.feed(self.mvec_to_muts(batch),
                                self.mvec_to_refs(batch))
            filter_vec.close()
            # Get the positions and names of the reads after filtering.
            self.positions = np.array(filter_vec.pos_kept)
            self.read_names = pd.Index(filter_vec.reads_kept)
            # Copy the filter parameters.
            self.filter_summary = filter_vec.summary

    @property
    def n_pos_init(self):
        return self.section.length

    @property
    def n_pos_polya(self):
        return len(self.pos_polya)

    @property
    def n_pos_gu(self):
        return len(self.pos_gu)

    @property
    def n_pos_user(self):
        return len(self.pos_user)

    @property
    def n_pos_kept(self) -> int:
        return self.positions.size

    @cached_property
    def n_reads_kept(self) -> int:
        if self.read_names is not None:
            return self.read_names.size
        # If the read names are unknown, then need to count the reads.
        return sum(batch.shape[0] for batch in
                   self.loader.iter_batches(self.section.positions[0: 1]))

    @property
    def min_mut_gap(self) -> int:
        return self.filter_summary["min_mut_gap"]

    def mvec_to_muts(self, mvec: pd.DataFrame):
        """ Compute bit vectors of mutations. """
        return BitCaller.mut(self.count_del, self.count_ins).call(mvec)

    def mvec_to_refs(self, mvec: pd.DataFrame):
        """ Compute bit vectors of reference matches. """
        return BitCaller.ref().call(mvec)

    def iter_muts_refs(self, muts: bool = False, refs: bool = True):
        """ For each batch of mutation vectors, select the reads and
        positions that passed the filters and yield a boolean data frame
        indicating the mutations and/or reference matches. """
        # Create a Series of uninitialized values with its index set to
        # the names of the reads. Only the index is needed, for the
        # operation pd.concat(), so the values are never used.
        read_series = (None if self.read_names is None
                       else pd.Series(np.empty_like(self.read_names),
                                      index=self.read_names))
        # Iterate over each batch of mutation vectors.
        for batch in self.loader.iter_batches(self.positions):
            if read_series is None:
                # Yield mutations and/or matches in all reads.
                reads = batch
            else:
                # Find the names of the reads in the current batch
                # that passed the filter. First, create a Series whose
                # index names all reads passing the filter (read_series)
                # and another Series whose index names all reads in the
                # current batch (pd.Series(index=batch.index)). Using
                # pd.concat(join="inner"), compute a new DataFrame whose
                # index is the intersection of the indexes of those two
                # Series. Finally, the index of that DataFrame names all
                # reads in the current batch that passed the filter.
                reads_passing = pd.concat([pd.Series(index=batch.index),
                                           read_series],
                                          axis=1, join="inner").index
                # Yield mutations and/or matches in only passing reads.
                reads = batch.loc[reads_passing]
            # Yield the mutations and/or matches as specified.
            if muts and refs:
                yield self.mvec_to_muts(reads), self.mvec_to_refs(reads)
            elif muts:
                yield self.mvec_to_muts(reads)
            elif refs:
                yield self.mvec_to_refs(reads)

    def iter_muts(self):
        yield from self.iter_muts_refs(muts=True, refs=False)

    def iter_refs(self):
        yield from self.iter_muts_refs(muts=False, refs=True)

    def all_muts(self):
        """ Return a boolean data frame indicating the mutated positions
        of every read that passed the filters. """
        return pd.concat(list(self.iter_muts()), axis=0)

    def all_refs(self):
        """ Return a boolean data frame indicating the matched positions
        of every read that passed the filters. """
        return pd.concat(list(self.iter_refs()), axis=0)

    @cache
    def _get_unique(self):
        """ Find the unique bit vectors and return the indexes """
        # Count each unique bit vector.
        unique, inverse, counts = np.unique(self.all_muts().values, axis=0,
                                            return_inverse=True,
                                            return_counts=True)
        # For each position, find the indexes of the unique bit vectors
        # with a mutation. Storing only the mutations requires much less
        # memory than storing the entire sparse matrix (unique) because
        # mutations are relatively rare.
        mut_idxs = tuple(map(np.flatnonzero, unique.T))
        return mut_idxs, counts, inverse

    @property
    def uniq_muts(self) -> tuple[np.ndarray, ...]:
        return self._get_unique()[0]

    @property
    def uniq_counts(self) -> np.ndarray:
        return self._get_unique()[1]

    @property
    def uniq_inverse(self) -> np.ndarray:
        return self._get_unique()[2]

    @property
    def n_uniq(self):
        return self.uniq_counts.size

    @staticmethod
    def summary_keys():
        return ["count_del", "count_ins",
                "exclude_polya", "exclude_gu", "exclude_pos",
                "n_pos_init", "n_pos_polya", "n_pos_gu",
                "n_pos_user", "n_pos_kept", "n_reads_kept"]

    @property
    def summary(self):
        return {key: self.__getattribute__(key) for key in self.summary_keys()}

    def to_dict(self):
        """ Return a dictionary containing bit vector report data. """
        data: dict[str, Any] = dict()
        # Sample attributes.
        data["sample"] = self.loader.sample
        # Section attributes.
        data.update(self.section.to_dict())
        # BitVector attributes.
        data.update(self.summary)
        # VectorFilter attributes.
        data.update(self.filter_summary)
        return data

    def __str__(self):
        return f"Bit Vectors of '{self.loader.sample}' over {self.section}"
