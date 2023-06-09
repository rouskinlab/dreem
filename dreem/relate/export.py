from logging import getLogger

import pandas as pd

from ..core.sect import index_to_seq_pos
from ..core.rel import translate_relvec


logger = getLogger(__name__)


def as_iter(vectors: pd.DataFrame, reference: bool = False):
    """ For each vector (row) in vectors, yield a string of its name and
    its mutations as human-readable text.

    Parameters
    ----------
    vectors: DataFrame
        The vectors to display as a DataFrame of bytes, typically from
        the method get_batch(), get_all_batches(), or get_all_vectors()
        of a .load.VectorLoader instance.
    reference: bool = False
        Whether to yield the reference sequence as the first item, prior
        to any mutation vectors. The reference sequence will be the same
        length as every mutation vector. This option can be True only if
        the columns of vectors include the nucleotide and the position,
        e.g. "G54"; otherwise, the reference sequence will be skipped
        and an error logged.
    """
    if reference:
        # Display the reference sequence above the vectors.
        try:
            # Get the reference sequence from the column names.
            seq, _ = index_to_seq_pos(vectors.columns.tolist())
            # Prepend the reference sequence to the lines of vectors.
            yield f"Reference\t{seq.decode()}"
        except Exception as error:
            logger.error(f"Could not determine sequence from columns of the "
                         f"vectors (perhaps you used numeric=True): {error} ")
    for index, row in zip(vectors.index, vectors.values, strict=True):
        yield f"{index}\t{translate_relvec(row).decode()}"


def as_block(vectors: pd.DataFrame, reference: bool = False):
    """ Display all mutation vectors as a block of human-readable text,
    with each vector on a new line and all the positions in the vector
    aligned vertically. Parameters are the same as those of as_iter. """
    return "\n".join(as_iter(vectors, reference=reference))
