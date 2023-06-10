from functools import cache
from logging import getLogger
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .write import RelateReport
from ..core import path
from ..core.load import DataLoader
from ..core.report import NumBatchF, SeqF
from ..core.sect import Section, RefSections, seq_pos_to_index
from ..core.seq import DNA

logger = getLogger(__name__)

POS = "by_position"
VEC = "by_vector"
READ = "Read Name"


class RelateLoader(DataLoader):
    """ Load batches of relation vectors. Wrapper around RelateReport
    that exposes only the attributes of the report that are required for
    loading batches of relation vectors. """

    IDX_COL = "__index_level_0__"

    @classmethod
    def report_type(cls):
        return RelateReport

    @property
    def seq(self):
        return self._rep.get_field(SeqF)

    @property
    def n_batches(self):
        return self._rep.get_field(NumBatchF)

    @property
    def sect(self):
        return None

    @property
    def section(self):
        return self.get_section()

    @cache
    def get_section(self, end5: int | None = None, end3: int | None = None):
        return Section(ref=self.ref, refseq=self.seq, end5=end5, end3=end3)

    def load_rel_batch(self, batch: int, positions: np.ndarray | None = None):
        """
        Return the relation vectors from one batch. Optionally, select
        a subset of the columns of the mutation vectors.

        Parameters
        ----------
        batch: int
            Number of the batch
        positions: np.ndarray | None = None
            Select 1-indexed positions to return. If None, return all.

        Return
        ------
        DataFrame
            Mutation vectors; each row is a vector indexed by its name,
            each column a position indexed by its positional number
        """
        # Read the vectors from the ORC file using PyArrow as backend.
        columns = (None if positions is None
                   else [self.IDX_COL] + seq_pos_to_index(self.seq,
                                                          positions,
                                                          start=1).to_list())
        vectors = pd.read_orc(self._rep.get_batch_path(batch), columns=columns)
        # Remove the column of read names and set it as the index.
        vectors.set_index(self.IDX_COL, drop=True, inplace=True)
        # Convert the index from bytes to str and give it a name.
        vectors.set_index(pd.Index(vectors.index.map(bytes.decode), name=READ),
                          inplace=True)
        # The vectors are stored as signed 8-bit integers (np.int8) and
        # must be cast to unsigned 8-bit integers (np.uint8) so that the
        # bitwise operations work. This step must be done after removing
        # the column of read names (which cannot be cast to np.uint8).
        return vectors.astype(np.uint8, copy=False)

    def iter_rel_batches(self, positions: Sequence[int] | None = None):
        """
        Yield every batch of relation vectors as a DataFrame.

        Parameters
        ----------
        positions: Sequence[int] | None = None
            Select 1-indexed positions to return. If None, return all.
        """
        for batch in range(self.n_batches):
            yield self.load_rel_batch(batch, positions)

    def __str__(self):
        return f"Mutation Vectors from '{self.sample}' aligned to '{self.ref}'"


def open_reports(report_files: Iterable[Path]):
    """ Load an arbitrary number of vector reports. """
    reports = dict()
    for report_file in report_files:
        try:
            # Load the report and collect basic information.
            report = RelateLoader.open(report_file)
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
