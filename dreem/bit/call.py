from functools import cache
from logging import getLogger
from sys import byteorder
from typing import Any, Generator, Iterable

import numpy as np
import pandas as pd

from .filt import VectorFilter
from ..mut.load import VectorLoader
from ..util.sect import filter_gu, filter_polya, filter_pos, Section
from ..util.seq import MATCH, DELET, INS_5, INS_3, SUB_N

logger = getLogger(__name__)


QEQ = "="
QSB = "<"
QSP = ">"
QEB = "{"
QEP = "}"
QMETHOD = QSB, QEB, QEQ, QEP, QSP


def muts_to_bits(mut_vectors: pd.DataFrame,
                 query: int,
                 rel: str = QEQ) -> pd.DataFrame:
    """
    Return a boolean array of the same shape as vectors where element
    i,j is True if and only if the byte at element i,j of vectors meets
    the requested relationship between it and the query byte.

    Parameters
    ----------
    mut_vectors: DataFrame
        Mutation vectors
    query: int
        Byte to query; must be in range 0 - 255
    rel: str = "equals"
        Method to decide whether a byte counts, as follows:
        - "=": count only bytes equal to query
        - "<": count strict bitwise subsets of query
        - ">": count strict bitwise supersets of query
        - "{": count bitwise subsets of and bytes equal to query
        - "}": count bitwise supersets of and bytes equal to query

    Returns
    -------
    DataFrame
        Boolean type DataFrame of the same shape as vectors where each
        element is True if and only if the element at the same position
        in vectors fulfilled the relationship with query
    """
    # Validate the query byte.
    if not isinstance(query, int):
        raise TypeError(
            f"Expected query of type int, but got type {type(query).__name__}")
    try:
        query.to_bytes(1, byteorder)
    except OverflowError:
        raise ValueError(f"Expected query in range 0 - 255, but got {query}")
    if rel == QEQ:
        return np.equal(mut_vectors, query)
    if rel == QEB:
        return np.equal(np.bitwise_or(mut_vectors, query), query)
    if rel == QEP:
        return np.equal(np.bitwise_and(mut_vectors, query), query)
    if rel == QSB:
        return np.logical_and(muts_to_bits(mut_vectors, query, QEB),
                              np.not_equal(mut_vectors, query))
    if rel == QSP:
        return np.logical_and(muts_to_bits(mut_vectors, query, QEP),
                              np.not_equal(mut_vectors, query))
    raise ValueError(f"Parameter 'rel' must be in {QMETHOD}, but got '{rel}'")


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
        self.n_pos_init = self.section.length
        # Exclude poly(A) sequences from the section.
        filt_seq_polya, filt_pos_polya = filter_polya(self.exclude_polya,
                                                      self.section.seq,
                                                      self.section.positions)
        self.n_pos_polya = self.section.positions.size - len(filt_pos_polya)
        # Exclude Gs and Us from the section, if exclude_gu is True.
        filt_seq_gu, filt_pos_gu = (filter_gu(filt_seq_polya, filt_pos_polya)
                                    if self.exclude_gu
                                    else (filt_seq_polya, filt_pos_polya))
        self.n_pos_gu = len(filt_pos_polya) - len(filt_pos_gu)
        # Exclude arbitrary, user-specified positions from the section.
        filt_seq_user, filt_pos_user = filter_pos(self.exclude_pos,
                                                  filt_seq_gu,
                                                  filt_pos_gu)
        self.n_pos_user = len(filt_pos_gu) - len(filt_pos_user)
        if filter_vec is None:
            # Use all the filtered positions and all the reads.
            self.positions = np.array(filt_pos_user, dtype=int)
            self.read_names = None
            # There are no filter parameters.
            self.filter_summary = VectorFilter.null_summary()
        else:
            # Determine the reads and postions to filter.
            for batch in loader.iter_batches(filt_pos_user, numeric=True):
                # Compute the mutated and matched bits in the batch and
                # feed them through the filter.
                filter_vec.feed(self.mvec_to_muts(batch),
                                self.mvec_to_refs(batch))
            filter_vec.close()
            # Get the positions and names of the reads after filtering.
            self.positions = np.array(filter_vec.pos_kept, dtype=int)
            self.read_names = pd.Index(filter_vec.reads_kept.keys())
            # Copy the filter parameters.
            self.filter_summary = filter_vec.summary

    @property
    def n_pos_kept(self) -> int:
        return self.positions.size

    @property
    def n_reads_kept(self) -> int:
        return self.read_names.size

    @property
    def min_mut_gap(self) -> int:
        return self.filter_summary["min_mut_gap"]

    @property
    def _query_mut(self):
        query = SUB_N
        if self.count_del:
            query |= DELET
        if self.count_ins:
            query |= INS_3
        return query, QEB

    @property
    def _query_ref(self):
        return MATCH | INS_5, QEB

    def mvec_to_muts(self, mvec: pd.DataFrame):
        """ Compute bit vectors of mutations. """
        return muts_to_bits(mvec, *self._query_mut)

    def mvec_to_refs(self, mvec: pd.DataFrame):
        """ Compute bit vectors of reference matches. """
        return muts_to_bits(mvec, *self._query_ref)

    def iter_muts_refs(self, muts: bool = True, refs: bool = True):
        """ For each batch of mutation vectors, select the reads and
        positions that passed the filters and yield a boolean data frame
        indicating the mutations and/or reference matches. """
        # Create a Series of uninitialized floats with its index set to
        # the names of the reads. Only the index is needed, for the
        # operation pd.concat().
        read_series = pd.Series(np.empty_like(self.read_names),
                                index=self.read_names)
        # Iterate over each batch of mutation vectors.
        for batch in self.loader.iter_batches(self.positions, numeric=True):
            if self.read_names is None:
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

    def iter_muts(self) -> Generator[pd.DataFrame, Any, None]:
        yield from self.iter_muts_refs(muts=True, refs=False)

    def iter_refs(self) -> Generator[pd.DataFrame, Any, None]:
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
