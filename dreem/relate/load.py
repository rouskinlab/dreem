from logging import getLogger
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .report import RelateReport
from ..core.bitc import BitCaller
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

    def get_refseq(self):
        return self._report.get_field(SeqF)

    def load_data_personal(self, batch_file: Path, *,
                           positions: Sequence[int] | None = None):
        """
        Return the relation vectors from one batch. Optionally, return
        a subset of the positions (columns) in the relation vectors.

        Parameters
        ----------
        batch_file: Path
            File of the batch

        Return
        ------
        DataFrame
            Relation vectors; each row is a vector indexed by its name,
            each column a position indexed by its base and number.
        """
        # Determine the columns to load.
        if positions is None:
            # Load all columns.
            columns = None
        else:
            # Load the columns corresponding to the given positions.
            columns = seq_pos_to_index(self.seq, positions, start=1)
        # Read the batch file using the selected columns.
        vectors = pd.read_parquet(batch_file, columns=columns)
        # Convert the index from bytes to str and give it a name.
        vectors.set_index(pd.Index(vectors.index.map(bytes.decode), name=READ),
                          inplace=True)
        # The vectors are stored as signed 8-bit integers (np.int8) and
        # must be cast to unsigned 8-bit integers (np.uint8) so that the
        # bitwise operations work. This step must be done after removing
        # the column of read names (which cannot be cast to np.uint8).
        return vectors.astype(np.uint8, copy=False)

    def iter_batches_personal(self, *, positions: Sequence[int] | None = None):
        yield from super().iter_batches_personal(positions=positions)

    def iter_batches_processed(self, *, positions: Sequence[int] | None = None):
        yield from super().iter_batches_processed(positions=positions)


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
