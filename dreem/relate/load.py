from functools import cache
from logging import getLogger
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .write import RelateReport
from ..core import path
from ..core.sect import Section, RefSections, seq_pos_to_cols
from ..core.seq import DNA

logger = getLogger(__name__)

POS = "by_position"
VEC = "by_vector"
READ = "Read"


class RelaVecLoader(object):
    """ Load batches of relation vectors. Wrapper around RelateReport
    that exposes only the attributes of the report that are required for
    loading batches of relation vectors. """

    INDEX_COL = "__index_level_0__"

    def __init__(self, report: RelateReport):
        self._rep = report

    # Select the necessary attributes of the Vector Report.

    @property
    def out_dir(self):
        return self._rep.out_dir

    @property
    def sample(self):
        return self._rep.sample

    @property
    def ref(self):
        return self._rep.ref

    @property
    def seq(self):
        return self._rep.seq

    @property
    def n_batches(self):
        return self._rep.n_batches

    # Define new methods.

    @cache
    def section(self, end5: int | None = None, end3: int | None = None):
        return Section(ref=self.ref, ref_seq=self.seq, end5=end5, end3=end3)

    def load_batch(self, batch: int, positions: Sequence[int] | None = None):
        """
        Return the mutation vectors from one batch. Optionally, select
        a subset of the columns of the mutation vectors.

        Parameters
        ----------
        batch: int
            Number of the batch
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
        vectors = pd.read_orc(self._rep.get_batch_path(batch), columns=cols)
        # Remove the column of read names and set it as the index.
        vectors.set_index(self.INDEX_COL, drop=True, inplace=True)
        # Convert the index from bytes to str and give it a name.
        vectors.set_index(pd.Index(vectors.index.map(bytes.decode), name=READ),
                          inplace=True)
        # The vectors are stored as signed 8-bit integers (np.int8) and
        # must be cast to unsigned 8-bit integers (np.uint8) so that the
        # bitwise operations work. This step must be done after removing
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
        for batch in range(self.n_batches):
            yield self.load_batch(batch, positions)

    @classmethod
    def open(cls, report_file: Path):
        """ Create a VectorLoader from a vectoring report file. """
        return cls(RelateReport.open(report_file))

    def __str__(self):
        return f"Mutation Vectors from '{self.sample}' aligned to '{self.ref}'"


def open_reports(report_files: Iterable[Path]):
    """ Load an arbitrary number of vector reports. """
    reports = dict()
    for report_file in report_files:
        try:
            # Load the report and collect basic information.
            report = RelaVecLoader.open(report_file)
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
    report_files = path.find_files_multi(report_paths, [path.RelateRepSeg])
    reports = open_reports(report_files)
    sections = RefSections({(rep.ref, rep.seq) for rep in reports},
                           coords=coords, primers=primers, primer_gap=primer_gap,
                           library=library)
    return reports, sections
