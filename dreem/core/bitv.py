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
    def ninfo_per_pos(self) -> pd.Series:
        """ Number of informative bits at each position. """
        return pd.Series()

    @property
    @abstractmethod
    def ninfo_per_read(self) -> pd.Series:
        """ Number of informative bits in each read. """
        return pd.Series()

    @property
    @abstractmethod
    def nmuts_per_pos(self) -> pd.Series:
        """ Number of mutated bits at each position. """
        return pd.Series()

    @property
    @abstractmethod
    def nmuts_per_read(self) -> pd.Series:
        """ Number of mutated bits in each read. """
        return pd.Series()

    @property
    def finfo_per_pos(self):
        """ Fraction of informative bits at each position (index) for
        each bit caller (column). """
        return self.ninfo_per_pos / self.reads.size

    @property
    def finfo_per_read(self):
        """ Fraction of informative bits in each read for each bit
        caller. """
        return self.ninfo_per_read / self.pos.size

    @property
    def fmuts_per_pos(self) -> pd.DataFrame:
        """ Fraction of mutated bits at each position (index) for each
        bit caller (column). """
        return self.nmuts_per_pos / self.ninfo_per_pos

    @property
    def fmuts_per_read(self):
        """ Fraction of mutated bits for each read. """
        return self.nmuts_per_read / self.ninfo_per_read

    def _check_duplicate_reads(self):
        """ Verify that no read name occurs more than once. """
        if self.reads.has_duplicates:
            dups = self.reads[self.reads.duplicated(keep="first")]
            raise ValueError(f"Duplicate read names: {dups}")


class BitMatrix(BitVectorBase, ABC):
    """ Bit vectors represented with two explicit boolean matrices of
    informative and mutated bits. """

    @property
    @abstractmethod
    def info(self):
        """ Boolean DataFrame indicating informative bits. """
        return pd.DataFrame()

    @property
    @abstractmethod
    def muts(self):
        """ Boolean DataFrame indicating mutated bits. """
        return pd.DataFrame()

    @property
    def _main_index(self):
        return self.info.columns

    @property
    def reads(self):
        return self.info.index

    @property
    def ninfo_per_pos(self):
        """ Number of informative bits for each position. """
        return pd.Series(np.count_nonzero(self.info, axis=0), index=self.seqpos)

    @property
    def ninfo_per_read(self):
        """ Number of informative bits for each read. """
        return pd.Series(np.count_nonzero(self.info, axis=1), index=self.reads)

    @property
    def nmuts_per_pos(self):
        """ Number of mutated bits for each position. """
        return pd.Series(np.count_nonzero(self.muts, axis=0), index=self.seqpos)

    @property
    def nmuts_per_read(self):
        """ Number of mutated bits for each read. """
        return pd.Series(np.count_nonzero(self.muts, axis=1), index=self.reads)


class BitBatch(BitMatrix):
    """ One batch of bit vectors. """

    def __init__(self, info: pd.DataFrame, muts: pd.DataFrame,
                 mask: dict[str, Callable[[BitBatch], pd.Index]] | None = None):
        logger.debug(f"Initializing {self}")
        # Initialize the reads.
        self._info = info
        self._muts = muts
        # Validate the reads.
        self._check_indexes()
        self._check_duplicate_reads()
        # Mask the reads.
        self._masked = Counter(self._mask(mask or dict()))

    def _check_indexes(self):
        """ Verify that the read names and positions in info and muts
        are consistent with each other. """
        logger.debug(f"Checking indexes of {self}")
        if not self.info.columns.equals(self.muts.columns):
            raise ValueError("info and muts must have the same columns, but "
                             f"got {self.info.columns} and {self.muts.columns}")
        if not self.info.index.equals(self.muts.index):
            raise ValueError("info and muts must have the same index, but got "
                             f"{self.info.index} and {self.muts.index}")

    @property
    def info(self):
        return self._info

    @property
    def muts(self):
        return self._muts

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
        self._info.drop(index=drop, inplace=True)
        self._muts.drop(index=drop, inplace=True)
        return drop.size


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
            raise ValueError(f"Attempted to close {self} with no batches")
        self._closed = True
        return True

    @abstractmethod
    def _add_info_muts(self, batch: BitBatch):
        """ Count the informative and mutated bits and add them to the
        accumulator. """
        return

    def add_batch(self, batch: BitBatch):
        """ Add one batch to the accumulator. """
        logger.debug(f"Adding batch {batch} to {self}")
        if self._closed:
            raise ValueError(f"Attempted to add batch to closed {self}")
        if self.nbatches > 0:
            # Confirm that the names of the masks in this batch match
            # the names of the masks in previous batches.
            if (bnm := sorted(batch.nmasked)) != (snm := sorted(self.nmasked)):
                raise ValueError(f"Names of masks in current ({bnm}) and "
                                 f"previous ({snm}) batches disagree")
        # Update the counts of the numbers of reads masked.
        self._masked += batch.nmasked
        logger.debug(f"Current masked counts for {self}: {self.nmasked}")
        # Add the counts of informative and mutated bits.
        self._add_info_muts(batch)
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
        # Initialize lists of each batch's informative and mutated bits.
        self._info: list[pd.DataFrame] | pd.DataFrame = list()
        self._muts: list[pd.DataFrame] | pd.DataFrame = list()
        super().__init__(batches)

    def close(self):
        if not super().close():
            return False
        # Merge the informative and mutated bits into data frames.
        self._info = pd.concat(self._info, axis=0)
        self._muts = pd.concat(self._muts, axis=0)
        # Confirm there are no duplicate reads among all batches.
        self._check_duplicate_reads()
        return True

    def _add_info_muts(self, batch: BitBatch):
        # Confirm that the positions in the current batch match those in
        # any previous batches.
        if self._info:
            # Compare this batch to the first batch.
            if not batch.seqpos.equals(ic := self._info[0].columns):
                # The positions disagree.
                raise ValueError(f"Positions in batch {len(self._info) + 1} "
                                 f"({batch.seqpos}) conflict with previous "
                                 f"batches ({ic})")
        # Add the informative and mutated bits from this batch to the
        # totals among all batches.
        self._info.append(batch.info)
        self._muts.append(batch.muts)

    @property
    def info(self):
        self.close()
        return self._info

    @property
    def muts(self):
        self.close()
        return self._muts

    def get_unique_muts(self):
        return UniqMutBits(self.muts.values)


class BitCounter(BitAccum):
    """ Accumulates batches of bit vectors into counts of informative
    and mutated bits per position and per read. """

    def __init__(self, batches: Iterable[BitBatch] = ()):
        # Initialize the total counts of informative and mutated bits.
        self._ninfo_per_pos: pd.Series | None = None
        self._nmuts_per_pos: pd.Series | None = None
        self._ninfo_per_read = list()
        self._nmuts_per_read = list()
        super().__init__(batches)

    def close(self):
        if not super().close():
            return False
        # Confirm there are no duplicate reads among all batches.
        self._check_duplicate_reads()
        return True

    def _add_info_muts(self, batch: BitBatch):
        logger.debug(f"Adding batch {len(self._ninfo_per_read) + 1} to {self}")
        # Confirm that the positions in the current batch match those in
        # previous batches.
        if self._ninfo_per_pos is None:
            # This is the first batch. Initialize counts to zero.
            self._ninfo_per_pos = pd.Series(0, index=batch.info.columns)
            self._nmuts_per_pos = pd.Series(0, index=batch.muts.columns)
        else:
            # Compare this batch to the first batch.
            if not batch.seqpos.equals(mi := self._main_index):
                # The positions disagree.
                raise ValueError(
                    f"Positions in batch {len(self._ninfo_per_read) + 1} "
                    f"({batch.seqpos}) conflict with previous batches ({mi})")
        # Add the counts for this batch to the totals.
        self._ninfo_per_pos += batch.ninfo_per_pos
        self._nmuts_per_pos += batch.nmuts_per_pos
        self._ninfo_per_read.append(batch.ninfo_per_read)
        self._nmuts_per_read.append(batch.nmuts_per_read)
        logger.debug(f"Added batch {len(self._ninfo_per_read)} to {self}")
        logger.debug(f"Counts:\n{self._ninfo_per_pos}\n{self._nmuts_per_pos}")

    @property
    def _main_index(self):
        return self._ninfo_per_pos.index

    @property
    def ninfo_per_pos(self) -> pd.Series:
        return self._ninfo_per_pos

    @cached_property
    def ninfo_per_read(self) -> pd.Series:
        return pd.concat(self._ninfo_per_read, axis=0)

    @property
    def nmuts_per_pos(self) -> pd.Series:
        return self._nmuts_per_pos

    @cached_property
    def nmuts_per_read(self) -> pd.Series:
        return pd.concat(self._nmuts_per_read, axis=0)

    @property
    def reads(self):
        """ Read names. """
        return self.ninfo_per_read.index

    @property
    def nreads(self):
        return sum(batch.size for batch in self._ninfo_per_read)

    @property
    def read_batches(self):
        """ Read names in each batch. """
        return (batch.index for batch in self._ninfo_per_read)

    @property
    def seqpos(self):
        return self._ninfo_per_pos.index
