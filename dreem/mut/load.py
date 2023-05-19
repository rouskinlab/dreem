from logging import getLogger
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .write import iter_batch_paths, VectorReport
from ..util import path
from ..util.sect import Section, RefSections, seq_pos_to_cols
from ..util.seq import DNA

logger = getLogger(__name__)


POS = "by_position"
VEC = "by_vector"
READ = "Read"


class VectorLoader(object):
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

    @property
    def n_batches(self):
        return len(self.checksums)

    def section(self, end5: int | None = None, end3: int | None = None):
        return Section(ref=self.ref, ref_seq=self.seq, end5=end5, end3=end3)

    def get_batch(self,
                  batch_file: Path,
                  positions: Sequence[int] | None = None):
        """
        Return the mutation vectors from one batch. Optionally, select
        a subset of the columns of the mutation vectors.

        Parameters
        ----------
        batch_file: Path
            Path to the batch of mutation vectors
        positions: Sequence[int] | None = None
            Select 1-indexed positions to return. If None, return all.

        Return
        ------
        DataFrame
            Mutation vectors; each row is a vector indexed by its name,
            each column a position indexed by its positional number
        """
        # Read the vectors from the ORC file using PyArrow as backend.
        cols = (None if positions is None
                else [self.INDEX_COL] + seq_pos_to_cols(self.seq, positions))
        vectors = pd.read_orc(batch_file, columns=cols)
        # Remove the column of read names and set it as the index.
        vectors.set_index(self.INDEX_COL, drop=True, inplace=True)
        # Convert the index from bytes to str and give it a name.
        vectors.set_index(pd.Index(vectors.index.map(bytes.decode),
                                   name=READ),
                          inplace=True)
        # The vectors are stored as signed 8-bit integers (np.int8) and
        # must be cast to unsigned 8-bit integers (np.uint8) so that the
        # bitwise operations work. This step must be doneafter removing
        # the column of read names (which cannot be cast to np.uint8).
        return vectors.astype(np.uint8, copy=False)

    def iter_batches(self, positions: Sequence[int] | None = None):
        """
        Yield every batch of mutation vectors as a DataFrame.

        Parameters
        ----------
        positions: Sequence[int] | None = None
            Select 1-indexed positions to return. If None, return all.
        """
        for batch, file in iter_batch_paths(self.out_dir, self.sample,
                                            self.ref, self.n_batches):
            yield self.get_batch(file, positions)

    def all_vectors(self, positions: Sequence[int] | None = None):
        """
        Return all mutation vectors for this vector reader. Note that
        reading all vectors could take more than the available memory
        and cause the program to crash. Thus, use this method only if
        all vectors will fit into memory. Otherwise, use the method
        ```get_all_batches``` to process the vectors in small batches.

        Parameters
        ----------
        positions: Sequence[int] | None = None
            Select 1-indexed positions to return. If None, return all.

        Returns
        -------
        DataFrame
            Mutation vectors; each row is a vector indexed by its name,
            each column a position indexed by its base and/or number
        """
        if not self.n_batches:
            # If there are no batches, then return an empty DataFrame.
            if positions is None:
                positions = self.section().positions
            return pd.DataFrame(columns=seq_pos_to_cols(self.seq, positions))
        # Load and concatenate every vector batch into one DataFrame.
        return pd.concat(self.iter_batches(positions), axis=0)

    @classmethod
    def open(cls, report_file: Path):
        """ Create a VectorLoader from a vectoring report file. """
        rep = VectorReport.open(report_file)
        return cls(out_dir=path.parse(report_file,
                                      *VectorReport.path_segs())[path.TOP],
                   sample=rep.sample, ref=rep.ref, seq=rep.seq,
                   n_vectors=rep.n_vectors, checksums=rep.checksums)

    def __str__(self):
        return f"Mutation Vectors from '{self.sample}' aligned to '{self.ref}'"


def open_reports(report_files: Iterable[Path]):
    """ Load an arbitrary number of vector reports. """
    reports = dict()
    for report_file in report_files:
        try:
            # Load the report and collect basic information.
            report = VectorLoader.open(report_file)
            key = report.sample, report.ref
            if key in reports:
                logger.warning(f"Got multiple reports for {key}")
            else:
                reports[key] = report
        except Exception as error:
            logger.error(f"Failed to open {report_file}: {error}")
    return list(reports.values())


def open_sections(report_paths: Iterable[Path],
                  coords: Iterable[tuple[str, int, int]],
                  primers: Iterable[tuple[str, DNA, DNA]],
                  primer_gap: int,
                  library: Path | None = None):
    report_files = path.find_files_multi(report_paths, [path.VecRepSeg])
    reports = open_reports(report_files)
    sections = RefSections({(rep.ref, rep.seq) for rep in reports},
                           coords=coords, primers=primers, primer_gap=primer_gap,
                           library=library)
    return reports, sections
