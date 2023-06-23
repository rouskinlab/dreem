"""
DREEM -- Bit Vector Module
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections import Counter
from functools import cached_property
from logging import getLogger
from typing import Callable, Iterable

import numpy as np
import pandas as pd

from .sect import index_to_pos

logger = getLogger(__name__)


class CloseEmptyBitAccumError(Exception):
    """ Close an empty BitAccum. """


class FeedClosedBitAccumError(Exception):
    """ Feed another batch to a closed BitAccum. """


class InconsistentIndexesError(Exception):
    """ Indexes of relation vectors do not match. """


class DuplicateIndexError(Exception):
    """ A value in an index is duplicated. """


class UniqMutBits(object):
    """ Collection of unique bit vectors indicating only mutations. """

    def __init__(self, muts: np.ndarray):
        """ Find the unique bit vectors. """
        if muts.ndim != 2:
            raise ValueError(f"muts needs 2 dimensions, but got {muts.ndim}")
        # Count each unique bit vector.
        uniq, self._inv, self._count = np.unique(muts, axis=0,
                                                 return_inverse=True,
                                                 return_counts=True)
        # For each position, find the indexes of the unique bit vectors
        # with a mutation. Storing only the mutations requires much less
        # memory than storing the entire sparse matrix (uniq) because
        # mutations are relatively rare.
        self._uniq_idxs = tuple(map(np.flatnonzero, uniq.T))

    @property
    def indexes(self):
        """ For each position, the indexes of all unique bit vectors
        that have a mutation at that position. """
        return self._uniq_idxs

    @property
    def counts(self) -> np.ndarray:
        """ Number of times each unique bit vector occurs in the given
        set of possibly redundant bit vectors. """
        return self._count

    @property
    def inverse(self) -> np.ndarray:
        """ Indexes to map the unique bit vectors back to the given set
        of possibly redundant bit vectors. """
        return self._inv

    @property
    def n_uniq(self):
        """ Number of unique bit vectors. """
        return self.counts.size

    @property
    def n_pos(self):
        """ Number of positions in each bit vector. """
        return len(self.indexes)

    def get_full(self):
        """ Full boolean matrix of the unique bit vectors. """
        # Initialize an all-False matrix with one row for each unique
        # bit vector and one column for each position.
        full = np.zeros((self.n_uniq, self.n_pos), dtype=bool)
        # For each position (j), set the mutated elements to True.
        for j, indexes in enumerate(self.indexes):
            full[indexes, j] = True
        return full

    def get_uniq_names(self):
        """ Return the unique bit vectors as byte strings. """
        # Get the full boolean matrix of the unique bit vectors and cast
        # the data from boolean to unsigned 8-bit integer type.
        chars = self.get_full().astype(np.uint8, copy=False)
        # Add ord('0') to transform every 0 into b'0' and every 1 into
        # b'1', and convert each row (bit vector) into a bytes object of
        # b'0' and b'1' characters.
        return pd.Index(np.apply_along_axis(np.ndarray.tobytes, 1,
                                            chars + ord('0')),
                        name="Bit Vector")

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return NotImplemented
        return (self.n_pos == other.n_pos
                and np.array_equal(self.counts, other.counts)
                and np.array_equal(self.inverse, other.inverse)
                and all(np.array_equal(self_idxs, other_idxs)
                        for self_idxs, other_idxs in zip(self.indexes,
                                                         other.indexes,
                                                         strict=True)))


class BitVectorBase(ABC):
    """ Base class for bit vectors. """

    @property
    @abstractmethod
    def nmasked(self) -> Counter[str]:
        """ Number of reads / bit vectors masked for each mask. """
        return Counter()

    @property
    @abstractmethod
    def reads(self):
        """ Names of the reads / bit vectors (after masking). """
        return pd.Index()

    @property
    def nreads(self):
        """ Number of reads / bit vectors (after masking). """
        return self.reads.size

    @property
    def nreads_given(self):
        """ Number of reads / bit vectors given (before masking). """
        return self.nreads + self.nmasked.total()

    @property
    @abstractmethod
    def _main_index(self):
        """ The main index of the bit vectors: index for Series, columns
        for DataFrames. """
        return pd.Index()

    @property
    def seqpos(self):
        """ Bases and positions in the bit vectors. """
        return self._main_index

    @property
    def pos(self):
        """ Numeric positions in the bit vectors. """
        return index_to_pos(self.seqpos)

    @property
    def npos(self):
        """ Number of positions (length of each bit vector). """
        return self.seqpos.size

    @property
    @abstractmethod
    def nall_per_pos(self) -> pd.Series:
        """ Number of informative bits at each position. """
        return pd.Series()

    @property
    @abstractmethod
    def nall_per_read(self) -> pd.Series:
        """ Number of informative bits in each read. """
        return pd.Series()

    @property
    @abstractmethod
    def nyes_per_pos(self) -> pd.Series:
        """ Number of affirmative bits at each position. """
        return pd.Series()

    @property
    @abstractmethod
    def nyes_per_read(self) -> pd.Series:
        """ Number of affirmative bits in each read. """
        return pd.Series()

    @property
    def fall_per_pos(self):
        """ Fraction of informative bits at each position. """
        return self.nall_per_pos / self.reads.size

    @property
    def fall_per_read(self):
        """ Fraction of informative bits in each read. """
        return self.nall_per_read / self.pos.size

    @property
    def fyes_per_pos(self):
        """ Fraction of affirmative bits at each position. """
        return self.nyes_per_pos / self.nall_per_pos

    @property
    def fyes_per_read(self):
        """ Fraction of affirmative bits for each read. """
        return self.nyes_per_read / self.nall_per_read

    def _check_duplicate_reads(self):
        """ Verify that no read name occurs more than once. """
        if self.reads.has_duplicates:
            dups = self.reads[self.reads.duplicated(keep="first")]
            raise DuplicateIndexError(f"Duplicate read names: {dups}")


class BitMatrix(BitVectorBase, ABC):
    """ Bit vectors represented with two explicit boolean matrices of
    informative and affirmative bits. """

    @property
    @abstractmethod
    def all(self):
        """ Boolean DataFrame indicating informative bits. """
        return pd.DataFrame()

    @property
    @abstractmethod
    def yes(self):
        """ Boolean DataFrame indicating affirmative bits. """
        return pd.DataFrame()

    @property
    def _main_index(self):
        return self.all.columns

    @property
    def reads(self):
        return self.all.index

    @property
    def nall_per_pos(self):
        """ Number of informative bits for each position. """
        return pd.Series(np.count_nonzero(self.all, axis=0), index=self.seqpos)

    @property
    def nall_per_read(self):
        """ Number of informative bits for each read. """
        return pd.Series(np.count_nonzero(self.all, axis=1), index=self.reads)

    @property
    def nyes_per_pos(self):
        """ Number of affirmative bits for each position. """
        return pd.Series(np.count_nonzero(self.yes, axis=0), index=self.seqpos)

    @property
    def nyes_per_read(self):
        """ Number of affirmative bits for each read. """
        return pd.Series(np.count_nonzero(self.yes, axis=1), index=self.reads)


class BitBatch(BitMatrix):
    """ One batch of bit vectors. """

    def __init__(self, info: pd.DataFrame, muts: pd.DataFrame,
                 mask: dict[str, Callable[[BitBatch], pd.Index]] | None = None):
        logger.debug(f"Initializing {self}")
        # Initialize the reads.
        self._all = info
        self._yes = muts
        # Validate the reads.
        self._check_indexes()
        self._check_duplicate_reads()
        # Mask the reads.
        self._masked = Counter(self._mask(mask or dict()))

    def _check_indexes(self):
        """ Verify that the read names and positions in info and muts
        are consistent with each other. """
        logger.debug(f"Checking indexes of {self}")
        if not self.all.columns.equals(self.yes.columns):
            raise InconsistentIndexesError(
                f"Got different columns for info ({self.all.columns}) "
                f"and muts ({self.yes.columns})")
        if not self.all.index.equals(self.yes.index):
            raise InconsistentIndexesError(
                f"Got different indexes for info ({self.all.index}) "
                f"and muts ({self.yes.index})")

    @property
    def all(self):
        return self._all

    @property
    def yes(self):
        return self._yes

    @property
    def nmasked(self):
        return self._masked

    def _mask(self, masks: dict[str, Callable[[BitBatch], np.ndarray]]):
        """ Drop reads selected by any of the masks, which should be
        boolean NumPy arrays. """
        logger.debug(f"Masking {self} with {masks}")
        return {name: self._drop(self.reads[mask(self)])
                for name, mask in masks.items()}

    def _drop(self, drop: pd.Index):
        """ Drop the reads in `drop`; return the number dropped. """
        logger.debug(f"Dropping {drop} from {self}")
        self._all.drop(index=drop, inplace=True)
        self._yes.drop(index=drop, inplace=True)
        return drop.size


class ClusterBitBatch(BitBatch):
    """ One batch of bit vectors with cluster membership weights. """

    def __init__(self, info: pd.DataFrame, muts: pd.DataFrame,
                 resps: pd.DataFrame):
        super().__init__(info, muts)
        self._resps = resps

    @property
    def clusters(self):
        """ Index of the clusters. """
        return self._resps.columns

    @property
    def nall_per_pos(self):
        # Count all for each position (row) and cluster (column).
        return self.all.T @ self._resps

    @property
    def nyes_per_pos(self):
        # Count yes for each position (row) and cluster (column).
        return self.yes.T @ self._resps

    @property
    def nall_per_read(self):
        # Count all for each read (row) and cluster (column).
        return self._resps.mul(super().nall_per_read, axis=0)

    @property
    def nyes_per_read(self):
        # Count yes for each read (row) and cluster (column).
        return self._resps.mul(super().nyes_per_read, axis=0)


class BitAccum(BitVectorBase, ABC):
    """ Accumulates batches of bit vectors. """

    def __init__(self, batches: Iterable[BitBatch] = ()):
        logger.debug(f"Initializing {self}")
        # Accept new batches.
        self._closed = False
        # Initialize the number of batches given.
        self._nbatches = 0
        # Initialize the numbers of reads masked.
        self._masked: Counter[str] = Counter()
        # If any batches were given to __init__, then add them.
        self.add_batches(batches)

    def close(self):
        """ Prevent the accumulator from accepting more batches, and
        finalize any attributes that need to be finalized. Return True
        if this is the first time `close()` is called, else False. """
        if self._closed:
            # Already closed: nothing to do.
            return False
        if self.nbatches == 0:
            raise CloseEmptyBitAccumError(
                f"Attempted to close {self} with no batches")
        self._closed = True
        return True

    @abstractmethod
    def _add_all_yes(self, batch: BitBatch):
        """ Count the informative and affirmative bits and add them to
        the accumulator. """

    def add_batch(self, batch: BitBatch):
        """ Add one batch to the accumulator. """
        logger.debug(f"Adding batch {batch} to {self}")
        if self._closed:
            raise FeedClosedBitAccumError(
                f"Attempted to add batch to closed {self}")
        # Update the counts of the numbers of reads masked.
        self._masked += batch.nmasked
        logger.debug(f"Current masked counts for {self}: {self.nmasked}")
        # Add the counts of informative and affirmative bits.
        self._add_all_yes(batch)
        # Increment the number of batches given.
        self._nbatches += 1

    def add_batches(self, batches: Iterable[BitBatch]):
        """ Add batches to the accumulator and, if at least one batch
        was given, close the accumulator to additional batches. """
        any_batches = False
        for batch in batches:
            any_batches = True
            self.add_batch(batch)
        if any_batches:
            self.close()

    @property
    def nmasked(self):
        return self._masked

    @property
    def nbatches(self):
        """ Number of batches given to the accumulator. """
        return self._nbatches


class BitMonolith(BitAccum, BitMatrix):
    """ Accumulates batches of bit vectors into one monolithic unit. """

    def __init__(self, batches: Iterable[BitBatch] = ()):
        # Initialize lists of each batch's total and affirmative bits.
        self._all: list[pd.DataFrame] | pd.DataFrame = list()
        self._yes: list[pd.DataFrame] | pd.DataFrame = list()
        super().__init__(batches)

    def close(self):
        if not super().close():
            return False
        # Merge the total and affirmative bits into data frames.
        self._all = pd.concat(self._all, axis=0)
        self._yes = pd.concat(self._yes, axis=0)
        # Confirm there are no duplicate reads among all batches.
        self._check_duplicate_reads()
        return True

    def _add_all_yes(self, batch: BitBatch):
        # Confirm that the positions in the current batch match those in
        # any previous batches.
        if self._all:
            # Compare this batch to the first batch.
            if not batch.seqpos.equals(ic := self._all[0].columns):
                # The positions disagree.
                raise InconsistentIndexesError(
                    f"Positions in batch {len(self._all) + 1} "
                    f"({batch.seqpos}) conflict with previous batches ({ic})")
        # Add the informative and affirmative bits from this batch to
        # the totals among all batches.
        self._all.append(batch.all)
        self._yes.append(batch.yes)

    @property
    def all(self):
        self.close()
        return self._all

    @property
    def yes(self):
        self.close()
        return self._yes

    def get_unique_muts(self):
        return UniqMutBits(self.yes.values)


class BitCounter(BitAccum):
    """ Accumulates batches of bit vectors into counts of informative
    and affirmative bits per position and per read. """

    def __init__(self, batches: Iterable[BitBatch] = ()):
        # Initialize the counts of informative and affirmative bits.
        self._nall_per_pos: pd.Series | pd.DataFrame | None = None
        self._nyes_per_pos: pd.Series | pd.DataFrame | None = None
        self._nall_per_read = list()
        self._nyes_per_read = list()
        super().__init__(batches)

    def close(self):
        if not super().close():
            return False
        # Confirm there are no duplicate reads among all batches.
        self._check_duplicate_reads()
        return True

    @classmethod
    def _zero_per_pos(cls, batch: BitBatch):
        """ Initialize a count of 0 for each position. """
        return pd.Series(0, index=batch.seqpos)

    def _add_all_yes(self, batch: BitBatch):
        logger.debug(f"Adding batch {len(self._nall_per_read) + 1} to {self}")
        # Confirm that the positions in the current batch match those in
        # previous batches.
        if self._nall_per_pos is None:
            # This is the first batch. Initialize counts to zero.
            self._nall_per_pos = self._zero_per_pos(batch)
            self._nyes_per_pos = self._zero_per_pos(batch)
        else:
            # Compare this batch to the first batch.
            if not batch.seqpos.equals(mi := self._main_index):
                # The positions disagree.
                raise InconsistentIndexesError(
                    f"Positions in batch {len(self._nall_per_read) + 1} "
                    f"({batch.seqpos}) conflict with previous batches ({mi})")
        # Add the counts for this batch to the totals.
        self._nall_per_pos += batch.nall_per_pos
        self._nyes_per_pos += batch.nyes_per_pos
        self._nall_per_read.append(batch.nall_per_read)
        self._nyes_per_read.append(batch.nyes_per_read)
        logger.debug(f"Added batch {len(self._nall_per_read)} to {self}")
        logger.debug(f"Counts:\n{self._nall_per_pos}\n{self._nyes_per_pos}")

    @property
    def _main_index(self):
        return self.nall_per_pos.index

    @property
    def nall_per_pos(self):
        return self._nall_per_pos

    @cached_property
    def nall_per_read(self):
        return pd.concat(self._nall_per_read, axis=0)

    @property
    def nyes_per_pos(self):
        return self._nyes_per_pos

    @cached_property
    def nyes_per_read(self) -> pd.Series:
        return pd.concat(self._nyes_per_read, axis=0)

    @property
    def reads(self):
        """ Read names. """
        return self.nall_per_read.index

    @property
    def nreads(self):
        return sum(batch.size for batch in self._nall_per_read)

    @property
    def read_batches(self):
        """ Read names in each batch. """
        return (batch.index for batch in self._nall_per_read)

    @property
    def seqpos(self):
        return self.nall_per_pos.index


class ClustBitCounter(BitCounter):
    """ Accumulates batches of bit vectors and cluster memberships into
    counts of informative and affirmative bits per position and per read
    for each cluster. """

    def add_batch(self, batch: ClusterBitBatch):
        if not isinstance(batch, ClusterBitBatch):
            raise TypeError(f"{self.__class__.__name__} expected batch of type "
                            f"'{ClusterBitBatch.__name__}', but got "
                            f"'{type(batch).__name__}'")
        return super().add_batch(batch)

    @classmethod
    def _zero_per_pos(cls, batch: ClusterBitBatch):
        """ Initialize a count of 0 for each position and cluster. """
        return pd.DataFrame(0, index=batch.seqpos, columns=batch.clusters)
