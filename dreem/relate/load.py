from logging import getLogger
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .report import RelateReport, BATCH_INDEX_COL
from ..core import path
from ..core.load import BatchLoader
from ..core.report import SeqF
from ..core.sect import seq_pos_to_index

logger = getLogger(__name__)

POS = "by_position"
VEC = "by_vector"
READ = "Read Name"


class RelateLoader(BatchLoader):
    """ Load batches of relation vectors. """

    @classmethod
    def get_report_type(cls):
        return RelateReport

    @property
    def seq(self):
        return self._rep.get_field(SeqF)

    def load_data(self, batch_file: Path, positions: np.ndarray | None = None):
        """
        Return the relation vectors from one batch. Optionally, return
        a subset of the positions (columns) in the relation vectors.

        Parameters
        ----------
        batch_file: Path
            File of the batch
        positions: np.ndarray | None = None
            Select 1-indexed positions to return. If None, return all.

        Return
        ------
        DataFrame
            Relation vectors; each row is a vector indexed by its name,
            each column a position indexed by its base and number.
        """
        # Select an engine to read the batch file using its extension.
        if batch_file.suffix in path.PARQ_EXTS:
            read_engine = pd.read_parquet
            set_index = False
        elif batch_file.suffix in path.ORC_EXTS:
            read_engine = pd.read_orc
            set_index = True
        else:
            raise ValueError(f"Relate batch file format '{batch_file.suffix}' "
                             f"is not supported for {batch_file}")
        # Determine the columns to read from the file.
        if positions is None:
            # Read all columns.
            columns = None
        else:
            # Read the columns indicated by the integer positions.
            columns = seq_pos_to_index(self.seq, positions, start=1).to_list()
            if set_index:
                # If the index is a regular column, then read it too.
                columns.insert(0, BATCH_INDEX_COL)
        # Read the batch file using the selected engine and columns.
        vectors = read_engine(batch_file, columns=columns)
        # If a column contains the read names, then set it as the index.
        try:
            vectors.set_index(BATCH_INDEX_COL, drop=True, inplace=True)
        except KeyError:
            # No column contains the read names: nothing must be done.
            pass
        # Convert the index from bytes to str and give it a name.
        vectors.set_index(pd.Index(vectors.index.map(bytes.decode), name=READ),
                          inplace=True)
        # The vectors are stored as signed 8-bit integers (np.int8) and
        # must be cast to unsigned 8-bit integers (np.uint8) so that the
        # bitwise operations work. This step must be done after removing
        # the column of read names (which cannot be cast to np.uint8).
        return vectors.astype(np.uint8, copy=False)

    def iter_batches(self, positions: np.ndarray | None = None):
        yield from super().iter_batches(positions=positions)


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
