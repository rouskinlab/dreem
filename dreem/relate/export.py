from logging import getLogger
from sys import byteorder

import pandas as pd

from ..core.sect import cols_to_seq_pos
from ..core.seq import (MATCH, DELET, INS_5, INS_3, SUB_A, SUB_C, SUB_G, SUB_T,
                        NOCOV, IRREC)


logger = getLogger(__name__)


# Number of unique bytes
N_BYTES = 256  # = 2^8
# Map each common byte in the vector encoding (e.g. MATCH) to a byte of
# human-readable text.
BYTE_SYMBOLS = {NOCOV: b"_",
                MATCH: b"=",
                DELET: b".",
                INS_5 | MATCH: b"{",
                INS_5 | MATCH | INS_3: b"+",
                INS_3 | MATCH: b"}",
                SUB_A: b"A",
                SUB_C: b"C",
                SUB_G: b"G",
                SUB_T: b"T",
                SUB_A | SUB_C | SUB_G | MATCH: b"N",
                SUB_T | SUB_A | SUB_C | MATCH: b"N",
                SUB_G | SUB_T | SUB_A | MATCH: b"N",
                SUB_C | SUB_G | SUB_T | MATCH: b"N",
                IRREC: b"!"}
# Default for uncommon bytes in the vector encoding.
OTHER_SYMBOL = b"?"
# Create an array that maps each vector byte to its readable character.
map_array = bytearray(OTHER_SYMBOL) * N_BYTES
for byte, symbol in BYTE_SYMBOLS.items():
    map_array[byte] = int.from_bytes(symbol, byteorder)


# Create a translation table from vector to human-readable encodings.
map_table = bytes.maketrans(bytes(range(N_BYTES)), map_array)


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
            seq, _ = cols_to_seq_pos(vectors.columns.tolist())
            # Prepend the reference sequence to the lines of vectors.
            yield f"Reference\t{seq.decode()}"
        except Exception as error:
            logger.error(f"Could not determine sequence from columns of the "
                         f"vectors (perhaps you used numeric=True): {error} ")
    for index, row in zip(vectors.index, vectors.values, strict=True):
        vector = row.tobytes(order='C').translate(map_table).decode()
        yield f"{index}\t{vector}"


def as_block(vectors: pd.DataFrame, reference: bool = False):
    """ Display all mutation vectors as a block of human-readable text,
    with each vector on a new line and all the positions in the vector
    aligned vertically. Parameters are the same as those of as_iter. """
    return "\n".join(as_iter(vectors, reference=reference))
