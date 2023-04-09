from functools import cache, cached_property
from itertools import chain
from logging import getLogger
from pathlib import Path
from sys import byteorder

import numpy as np
import pandas as pd

from .batch import iter_batch_paths
from .report import VectorReport
from ..util import path
from ..util.sect import Section, cols_to_seq_pos
from ..util.seq import (DNA, BLANK, EVERY,
                        MATCH, DELET, INS_5, INS_3,
                        SUB_A, SUB_C, SUB_G, SUB_T)

logger = getLogger(__name__)


class VectorReader(object):
    INDEX_COL = "__index_level_0__"

    def __init__(self, /, *,
                 sample: str,
                 ref: str,
                 seq: DNA,
                 out_dir: Path,
                 n_vectors: int,
                 checksums: list[str]):
        self.sample = sample
        self.ref = ref
        self.seq = seq
        self.out_dir = out_dir
        self.n_vectors = n_vectors
        self.checksums = checksums

    def get_section(self, end5: int, end3: int):
        return Section(ref=self.ref, end5=end5, end3=end3, ref_seq=self.seq)

    @cached_property
    def span(self):
        return self.get_section(1, len(self.seq))

    def get_batch(self, batch_file: Path, *,
                  section: Section | None = None,
                  numeric: bool = False):
        """
        Return the mutation vectors from one batch. Optionally, select
        a subset of the columns of the mutation vectors.

        Parameters
        ----------
        batch_file: Path
            Path to the batch of mutation vectors
        section: Section | None = None
            Section of the batch to return, or None for the entire batch
        numeric: bool = False
            Whether to convert the columns from {base}{position} strings
            to numeric (integer) values of just the positions, e.g.
            ['A14', 'G15', ...] if False, [14, 15, ...] if True

        Return
        ------
        DataFrame
            Mutation vectors; each row is a vector indexed by its name,
            each column a position indexed by its positional number
        """
        # Determine which columns to read from the file.
        if section is None:
            section = self.span
        columns = [self.INDEX_COL] + section.columns
        # Read the vectors from the ORC file using PyArrow as backend.
        vectors = pd.read_orc(batch_file, columns=columns)
        # Remove the column of read names and set it as the index.
        vectors.set_index(self.INDEX_COL, drop=True, inplace=True)
        # Convert the index from bytes to str and give it a name.
        vectors.set_index(pd.Index(vectors.index.map(bytes.decode),
                                   name="Read Name"),
                          inplace=True)
        if numeric:
            # Convert the remaining columns to their integer positions.
            vectors.columns = section.positions
        # The vectors are stored as signed 8-bit integers (np.int8) and
        # must be cast to unsigned 8-bit integers (np.uint8) so that the
        # bitwise operations work. This step must be doneafter removing
        # the column of read names (which cannot be cast to np.uint8).
        return vectors.astype(np.uint8, copy=False)

    def iter_batches(self, *,
                     section: Section | None = None,
                     numeric: bool = True):
        """
        Yield every batch of mutation vectors.

        Parameters
        ----------
        section: Section | None = None
            Section of the batch to return, or None for the entire batch
        numeric: bool (default: False)
            Whether to convert the columns from {base}{position} strings
            to numeric (integer) values of just the positions, e.g.
            ['A14', 'G15', ...] if False, [14, 15, ...] if True

        Yield
        -----
        DataFrame
            Mutation vectors; each row is a vector indexed by its name,
            each column a position indexed by its base and number
        """
        for batch, file in iter_batch_paths(self.out_dir, self.sample,
                                            self.ref, len(self.checksums)):
            yield self.get_batch(file, section=section, numeric=numeric)

    def get_all_vectors(self, *,
                        section: Section | None = None,
                        numeric: bool = True):
        """
        Return all mutation vectors for this vector reader. Note that
        reading all vectors could take more than the available memory
        and cause the program to crash. Thus, use this method only if
        all vectors will fit into memory. Otherwise, use the method
        ```get_all_batches``` to process the vectors in small batches.

        Parameters
        ----------
        section: Section | None = None
            Section of the batch to return, or None for the entire batch
        numeric: bool (default: False)
            Whether to convert the columns from {base}{position} strings
            to numeric (integer) values of just the positions, e.g.
            ['A14', 'G15', ...] if False, [14, 15, ...] if True

        Returns
        -------
        DataFrame
            Mutation vectors; each row is a vector indexed by its name,
            each column a position indexed by its base and number
        """
        # If there are no vectors, then return an empty DataFrame.
        if self.n_vectors == 0:
            if section is None:
                section = self.span
            columns = section.positions if numeric else section.columns
            return pd.DataFrame(columns=columns, dtype=int)
        # Load and concatenate every vector batch into one DataFrame.
        return pd.concat(self.iter_batches(section=section, numeric=numeric),
                         axis=0)

    @classmethod
    def query_vectors(cls, /,
                      vectors: pd.DataFrame,
                      query: int, *,
                      subsets: bool = False,
                      supersets: bool = False) -> pd.DataFrame:
        """
        Return a boolean array of the same shape as vectors where
        element i,j is True if and only if the byte at element i,j of
        vectors matches the given query byte. By default, a byte in
        vectors matches if it equals the query byte and is not blank
        (i.e. 00000000). Matches can be extended to include bitwise
        subsets and supersets of query by setting the corresponding
        parameters to True. Blank bytes in vectors never match.

        Parameters
        ----------
        vectors: DataFrame
            Mutation vectors
        query: int
            Byte to query; must be in range [0, 255]
        subsets: bool
            Whether to count non-blank bitwise subsets of the query
        supersets: bool
            Whether to count non-blank bitwise supersets of the query

        Return
        ------
        DataFrame
            Boolean type DataFrame of the same shape as vectors where
            each element is True if the element at the same position in
            vectors matched the query and False otherwise
        """
        if not isinstance(query, int):
            raise TypeError(
                f"Expected query of type int, but got {type(query).__name__}")
        if not BLANK <= query <= EVERY:
            raise ValueError(
                f"Expected query in range {BLANK} - {EVERY}, but got {query}")
        if supersets and subsets:
            # Count both supersets and subsets.
            return (cls.query_vectors(vectors, query, supersets=True)
                    | cls.query_vectors(vectors, query, subsets=True))
        if supersets:
            # Non-blank vector bytes that are matches and supersets of
            # the query byte count.
            if query == BLANK:
                # If query is BLANK (00000000), then every non-blank
                # byte is a superset.
                return vectors.astype(bool, copy=True)
            if query == EVERY:
                # If the query is EVERY (11111111), then no bitwise
                # supersets exist. Since subsets do not count, only
                # exact matches count.
                return cls.query_vectors(vectors, query)
            # No shortcut method can be used, so the supersets must be
            # computed explicitly. A vector byte is a match or superset
            # of the query byte iff both of the following are true:
            # - The bitwise intersection of the vector and query bytes
            #   equals the query byte, meaning that every bit set to 1
            #   in the query byte is also set to 1 in the vector byte,
            #   and thus the vector byte is a superset of the query.
            #   Equivalently, the union equals the vector byte.
            # - The vector byte is not blank. But since the query byte
            #   is not blank, if a vector byte satisfies the first
            #   condition, then it must also satisfy this one, so this
            #   latter condition need not be checked.
            return np.equal(np.bitwise_and(vectors, query), query)
        if query == BLANK:
            # If supersets do not count, then only matches and subsets
            # count. But if the query is BLANK (00000000), then there
            # are no subsets, and matches do not count because blank
            # vector bytes never count. Thus, no byte counts.
            return pd.DataFrame(False, dtype=bool,
                                index=vectors.index,
                                columns=vectors.columns)
        if subsets:
            # Non-blank vector bytes that are matches and subsets of the
            # query byte count.
            if query == EVERY:
                # If query is EVERY (11111111), then every non-blank
                # byte is a subset.
                return vectors.astype(bool, copy=True)
            if (query & (query - 1)) == 0:
                # If query is a power of 2, then it has exactly one bit
                # set to 1. Thus, the only possible subset of the query
                # is the blank byte, which never counts. Since supersets
                # do not count either, only exact matches count.
                return cls.query_vectors(vectors, query)
            # No shortcut method can be used, so the subsets must be
            # computed explicitly. A vector byte is a match or subset of
            # the query byte iff both of the following are true:
            # - The bitwise union of the vector and query bytes equals
            #   the query byte, meaning that there are no bits set to 1
            #   in the vector byte that are not 1 in the query byte,
            #   and thus the vector byte is a subset of the query.
            #   Equivalently, the intersection equals the vector byte.
            # - The vector byte is not blank.
            return (vectors.astype(bool, copy=False)
                    & np.equal(np.bitwise_or(vectors, query), query))
        # If neither subsets nor supersets count and query is not BLANK,
        # then count vector bytes that match the query byte exactly.
        return np.equal(vectors, query)

    @cache
    def count_muts_by_vec(self, /,
                          query: int, *,
                          subsets: bool = False,
                          supersets: bool = False,
                          section: Section | None = None) -> pd.Series:
        """
        Return the number of mutations that match the query for each
        vector in the mutational profile.

        Parameters
        ----------
        query: int
            Byte to query; must be in range [0, 255]
        subsets: bool
            Whether to count non-blank bitwise subsets of the query
        supersets: bool
            Whether to count non-blank bitwise supersets of the query
        section: Section | None = None
            Section of the batch to return, or None for the entire batch

        Return
        ------
        Series
            Number of mutations matching the query in each vector
        """
        # If there are no vectors, then return an empty Series.
        if not self.checksums:
            return pd.Series([], dtype=int)
        # Initialize empty list to count the mutations in each vector.
        counts = list()
        # Iterate over all batches of vectors.
        for vectors in self.iter_batches(section=section):
            # Count the number of mutations in each vector in the batch
            # and append them to the list of counts.
            counts.append(self.query_vectors(vectors,
                                             query,
                                             subsets=subsets,
                                             supersets=supersets).sum(axis=1))
        # Concatenate and return the number of mutations in each vector
        # among all batches.
        return pd.concat(counts, axis=0)

    @cache
    def count_muts_by_pos(self, /,
                          query: int, *,
                          subsets: bool = False,
                          supersets: bool = False,
                          section: Section | None = None,
                          numeric: bool = True) -> pd.Series:
        """
        Return the number of mutations that match the query at each
        position in the mutational profile.

        Parameters
        ----------
        query: int
            Byte to query; must be in range [0, 255]
        subsets: bool
            Whether to count non-blank bitwise subsets of the query
        supersets: bool
            Whether to count non-blank bitwise supersets of the query
        section: Section | None = None
            Section of the batch to return, or None for the entire batch
        numeric: bool (default: False)
            Whether to convert the columns from {base}{position} strings
            to numeric (integer) values of just the positions, e.g.
            ['A14', 'G15', ...] if False, [14, 15, ...] if True

        Return
        ------
        Series
            Number of mutations matching the query at each position
        """
        if section is None:
            section = self.span
        # Initialize empty Series to count mutations at each position.
        counts = pd.Series(np.zeros(len(section.seq), dtype=int),
                           index=(section.positions if numeric
                                  else section.columns))
        # Iterate over all batches of vectors.
        for vectors in self.iter_batches(section=section, numeric=numeric):
            # Add the number of mutations at each position in the batch
            # to the cumulative count of mutations at each position.
            counts += self.query_vectors(vectors,
                                         query,
                                         subsets=subsets,
                                         supersets=supersets).sum(axis=0)
        return counts

    def get_cluster_mus(self, /,
                        membership: pd.DataFrame,
                        query: int, *,
                        subsets: bool = False,
                        supersets: bool = False,
                        section: Section | None = None,
                        numeric: bool = False) -> pd.DataFrame:
        """
        Calculate the Mutation fraction at each position in a mutational
        profile for one or more clusters.

        Parameters
        ----------
        membership: DataFrame
            Cluster membership: each index (i) is the name of a read,
            each column (k) the name of a cluster, and each value (i, k)
            the likelihood that read (i) came from cluster (k).
        query: int
            Byte to query; must be in range [0, 255]
        subsets: bool
            Whether to count non-blank bitwise subsets of the query
        supersets: bool
            Whether to count non-blank bitwise supersets of the query
        section: Section | None = None
            Section of the batch to return, or None for the entire batch
        numeric: bool (default: False)
            Whether to convert the columns from {base}{position} strings
            to numeric (integer) values of just the positions, e.g.
            ['A14', 'G15', ...] if False, [14, 15, ...] if True

        Return
        ------
        DataFrame
            Mutation rates: each index (j) is a position in the profile,
            each column (k) the name of a cluster, and each value (j, k)
            the Mutation fraction at posision (j) in cluster (k).

        Explanation
        -----------
        To limit memory usage, this function is implemented on a single
        line. This is what each step does:

        1.  Read all mutation vectors into a DataFrame. Assume that if
            the reads were clustered, then they can all fit into memory
            at once. It is much easier to compute products of indexed
            DataFrames all at once than by summing over chunks:

            mutvectors = self.get_all_vectors(section=section,
                                              numeric=numeric)

        2.  Compute the bit vectors: a matrix of boolean values that
            indicate whether each read (row) is mutated at each position
            (column). The matrix has the same shape as mutvectors. The
            query determines which mutations count as True or False:

            bitvectors = self._query_vectors(mutvectors, query,
                                             subsets=subsets,
                                             supersets=supersets)

        3.  Compute the weighted number of mutations at each position in
            each cluster, weighed by the likelihood that each read came
            from each cluster. The resulting matrix has a row for each
            position and a column for each cluster:

            mutsums = bitvectors.T.dot(membership)

        4.  Compute the weighted number of reads in each cluster by just
            summing the likelihood that each read came from the cluster:

            readsums = membership.sum(axis=0)

        5.  Compute the mutation rates for each cluster by dividing the
            cluster-weighted number of mutations at each position by the
            weighted number of reads in the cluster:

            return mutsums / readsums
        """
        return (self.query_vectors(self.get_all_vectors(section=section,
                                                        numeric=numeric),
                                   query,
                                   subsets=subsets,
                                   supersets=supersets).T.dot(membership)
                / membership.sum(axis=0))

    @classmethod
    def load(cls, report_file: Path, validate_checksums: bool = True):
        rep = VectorReport.load(report_file, validate_checksums)
        out_dir = path.parse(report_file, path.ModSeg, path.SampSeg,
                             path.RefSeg, path.VecRepSeg)[path.TOP]
        return cls(out_dir=out_dir,
                   sample=rep[rep.SampleField],
                   ref=rep[rep.RefField],
                   seq=rep[rep.SeqField],
                   n_vectors=rep[rep.NumVectorsField],
                   checksums=rep[rep.ChecksumsField])


vector_trans_table = bytes.maketrans(*map(b"".join, zip(*[(
    i.to_bytes(length=1, byteorder=byteorder),
    (
        b"." if i == BLANK
        else b"~" if i == MATCH
        else b"/" if i == DELET
        else b"{" if i == (INS_5 | MATCH)
        else b"}" if i == (INS_3 | MATCH)
        else b"A" if i == SUB_A
        else b"C" if i == SUB_C
        else b"G" if i == SUB_G
        else b"T" if i == SUB_T
        else b"?"
    )
) for i in range(256)])))


def trans_vectors_iter(vectors: pd.DataFrame):
    for index, row in zip(vectors.index, vectors.values, strict=True):
        vector = row.tobytes(order='C').translate(vector_trans_table).decode()
        yield f"{index}\t{vector}\n"


def trans_vectors_block(vectors: pd.DataFrame, reference: bool = False):
    lines = trans_vectors_iter(vectors)
    if reference:
        # Display the reference sequence above the vectors.
        try:
            # Get the reference sequence from the column names.
            seq, _ = cols_to_seq_pos(vectors.columns.tolist())
            # Prepend the reference sequence to the lines of vectors.
            lines = chain([f"Reference\t{seq.decode()}\n"], lines)
        except Exception as error:
            logger.error(f"Could not determine sequence from columns of the "
                         f"vectors (perhaps you used numeric=True): {error} ")
    return "".join(lines)
