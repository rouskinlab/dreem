from pathlib import Path
from typing import Iterable

from ..call.load import CallVecLoader


def callvec(report_file: Path):
    # Load the called bit vectors.
    loader = CallVecLoader.open(report_file)
    mus = loader.mus
