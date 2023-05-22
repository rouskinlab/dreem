from collections import Counter
from logging import getLogger
from pathlib import Path
from typing import Iterable

import pandas as pd

from ..util.files import digest_file

logger = getLogger(__name__)


def write_batch(read_names: Iterable[str], batch_file: Path):
    """ Write the names of the reads in one batch to a file. """
    # Write the read names to the batch file.
    read_data = pd.Series(read_names)
    read_data.to_csv(batch_file, index=False)
    logger.debug(f"Wrote {read_data.size} read names to {batch_file}:\n"
                 f"{read_data}")
    # Compute the checksum of the batch file.
    return digest_file(batch_file)
