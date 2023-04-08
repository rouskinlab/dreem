from logging import getLogger
from pathlib import Path

import numpy as np
import pandas as pd

from ..util import path
from ..util.sect import seq_pos_to_cols


logger = getLogger(__name__)


BATCH_NUM_START = 1


def mib_to_bytes(batch_size: float):
    """
    Return the number of bytes per batch of a given size in mebibytes.

    Parameters
    ----------
    batch_size: float
        Size of the batch in mebibytes (1 MiB = 2^20 bytes)

    Return
    ------
    int
        Number of bytes per batch, to the nearest integer
    """
    return round(batch_size * 1048576)  # 1048576 = 2^20


def _get_batch_dir(out_dir: Path, sample: str, ref: str):
    return path.build(path.ModSeg, path.SampSeg, path.RefSeg,
                      top=out_dir, module=path.MOD_VECT,
                      sample=sample, ref=ref)


def _get_batch_path(out_dir: Path, sample: str, ref: str, batch: int):
    batch_seg = path.VecBatSeg.build({path.BATCH: batch,
                                      path.EXT: path.ORC_EXT})
    return _get_batch_dir(out_dir, sample, ref).joinpath(batch_seg)


def _get_batch_paths(out_dir: Path, sample: str, ref: str, n_batches: int):
    return {batch: _get_batch_path(out_dir, sample, ref, batch)
            for batch in range(BATCH_NUM_START, n_batches + BATCH_NUM_START)}


def _write_batch(batch: int,
                 vectors: tuple[bytearray, ...],
                 read_names: list[bytes], *,
                 sample: str,
                 ref: str,
                 seq: bytes,
                 out_dir: Path):
    """ Write a batch of mutation vectors to an ORC file. """
    logger.info(
        f"Began writing sample '{sample}' reference '{ref}' batch {batch}")
    # Process the mutation vectors into a 2D NumPy array.
    array = np.frombuffer(b"".join(vectors), dtype=np.byte)
    array.shape = (array.size // len(seq), len(seq))
    # Data must be converted to pd.DataFrame for PyArrow to write.
    # Set copy=False to prevent copying the mutation vectors.
    positions = np.arange(1, len(seq) + 1)
    mut_frame = pd.DataFrame(data=array,
                             index=read_names,
                             columns=seq_pos_to_cols(seq, positions),
                             copy=False)
    batch_path = _get_batch_path(out_dir, sample, ref, batch)
    batch_path.parent.mkdir(parents=True, exist_ok=True)
    mut_frame.to_orc(batch_path, index=True, engine="pyarrow")
    logger.info(f"Ended writing sample '{sample}' reference '{ref}' "
                f"batch {batch} to {batch_path}")
    return batch_path
