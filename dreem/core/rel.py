"""
DREEM Relation Core Module
--------------------------
Auth: Matty
"""

from collections import defaultdict
from itertools import (product, combinations,
                       combinations_with_replacement as cwr)
from sys import byteorder
from typing import Sequence

import numpy as np

from .sect import Section
from .seq import (BASES, A_INT, C_INT, G_INT, T_INT, N_INT, DNA,
                  expand_degenerate_seq)

# Minimum number of non-indel bases between an insertion and a deletion.
MIN_INDEL_GAP = 1

# Integer encodings for relation vectors
IRREC = b"\x00"[0]  # 00000000 (000): irreconcilable paired mates
MATCH = b"\x01"[0]  # 00000001 (001): match with reference
DELET = b"\x02"[0]  # 00000010 (002): deletion from reference
INS_5 = b"\x04"[0]  # 00000100 (004): insertion 5' of base in reference
INS_3 = b"\x08"[0]  # 00001000 (008): insertion 3' of base in reference
SUB_A = b"\x10"[0]  # 00010000 (016): substitution to A
SUB_C = b"\x20"[0]  # 00100000 (032): substitution to C
SUB_G = b"\x40"[0]  # 01000000 (064): substitution to G
SUB_T = b"\x80"[0]  # 10000000 (128): substitution to T
NOCOV = b"\xff"[0]  # 11111111 (255): not covered by read
SUB_N = SUB_A | SUB_C | SUB_G | SUB_T
ANY_N = SUB_N | MATCH
INS_B = INS_5 | INS_3
INDEL = DELET | INS_5 | INS_3
MINS5 = INS_5 | MATCH
MINS3 = INS_3 | MATCH
MINSB = INS_B | MATCH

# Map bases to integer encodings and vice versa.
BASE_ENCODINGS = {A_INT: SUB_A, C_INT: SUB_C, G_INT: SUB_G, T_INT: SUB_T}
BASE_DECODINGS = {code: base for base, code in BASE_ENCODINGS.items()}

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


def validate_relvec(relvec: np.ndarray):
    """ Confirm that a relation vector is valid, then return it. """
    if not isinstance(relvec, np.ndarray):
        raise TypeError(f"Expected {np.ndarray}, but got {type(relvec)}")
    if relvec.dtype.type is not np.uint8:
        raise TypeError(f"Expected an array of type {np.uint8}, but got "
                        f"type {relvec.dtype.type}")
    if relvec.ndim != 1:
        raise ValueError("Expected an array with 1 dimension, but got "
                         f"{relvec.ndim} dimensions")
    called = np.flatnonzero(relvec != NOCOV)
    if not called.any():
        raise ValueError("Relation vector is entirely blank")
    if (nexp := called[-1] - called[0] + 1) != called.size:
        raise ValueError(f"Expected {nexp} base calls "
                         f"({np.arange(called[0], called[-1] + 1)}), but got "
                         f"{called.size} ({called})")
    return relvec


def translate_relvec(relvec: np.ndarray):
    """ Translate a binary relation vector to human-readable text. """
    return validate_relvec(relvec).tobytes().translate(map_table)


def encode_relate(ref_base: int, read_base: int, read_qual: int, min_qual: int):
    """
    Encode the relation between a base in the read and a base in the
    reference sequence. If the read quality is sufficient, then return
    the match encoding if the read and reference bases match, otherwise
    the encoding of the substitution for the base in the read.
    If the read quality is insufficient, then return the fully ambiguous
    base encoding, that is a match or substitution to any base except
    the reference base, since a "substitution to the reference base"
    would be a match, not a substitution.

    Parameters
    ----------
    ref_base: int
        ASCII encoding for the base in the reference sequence.
    read_base: int
        ASCII encoding for the base in the read sequence.
    read_qual: int
        ASCII encoding for the Phred quality score of the read base.
    min_qual: int
        Minimum value of `read_qual` to not call the relation ambiguous.
    """
    return ((MATCH if ref_base == read_base
             else BASE_ENCODINGS[read_base]) if read_qual >= min_qual
            else ANY_N ^ BASE_ENCODINGS[ref_base])


def encode_match(read_base: int, read_qual: int, min_qual: int):
    """
    A more efficient version of `encode_compare` given prior knowledge
    from the CIGAR string that the read and reference match at this
    position. Note that there is no analagous version when there is a
    known substitution because substitutions are relatively infrequent,
    so optimizing their processing would speed the program only slightly
    while making the source code more complex and harder to maintain.
    """
    return (MATCH if read_qual >= min_qual
            else ANY_N ^ BASE_ENCODINGS[read_base])


def iter_relvecs_q53(ref_seq: DNA, low_qual: Sequence[int] = (),
                     end5: int | None = None, end3: int | None = None):
    """
    For a given reference sequence, yield every possible unambiguous
    relation vector between positions end5 and end3 that follows the
    given low-quality positions and has at most two insertions.

    Parameters
    ----------
    ref_seq: DNA
        Sequence of the reference.
    low_qual: Sequence[int]
        List of positions in the read that are low-quality.
    end5: int | None = None
        5' end of the read; 1-indexed with respect to `ref_seq`.
    end3: int | None = None
        3' end of the read; 1-indexed with respect to `ref_seq`.
    """
    low_qual = set(low_qual)
    # Determine the section of the reference sequence that is occupied
    # by the read.
    sect = Section("", ref_seq, end5=end5, end3=end3)
    if low_qual - set(sect.positions):
        raise ValueError(f"Invalid positions in low_qual: "
                         f"{sorted(low_qual - set(sect.positions))}")
    # Find the possible relationships at each position in the section,
    # not including insertions.
    rel_opts: list[tuple[int, ...]] = [(NOCOV,)] * len(ref_seq)
    for pos in sect.positions:
        # Find the base in the reference sequence (pos is 1-indexed).
        ref_base = ref_seq[pos - 1]
        if low_qual:
            # The only option is low-quality.
            opts = encode_relate(ref_base, ref_base, 0, 1),
        else:
            # The options are a match and three substitutions.
            opts = tuple(encode_relate(ref_base, read_base, 1, 1)
                         for read_base in BASES)
        # A deletion is always an option.
        rel_opts[pos - 1] = opts + (DELET,)
    # Iterate through all possible relationships at each position.
    margin = MIN_INDEL_GAP + 1
    for rel in product(*rel_opts):
        # Generate and yield a relation vector from the relationships.
        relvec = np.array(rel, dtype=np.uint8)
        yield relvec
        # Allow up to two insertions in one read.
        for ins1_pos5, ins2_pos5 in cwr(range(sect.end5 - 1, sect.end3 + 1), 2):
            # Skip insertions nearby any deletion.
            mrg15 = max(ins1_pos5 - margin, 0)
            mrg13 = min(ins1_pos5 + margin, relvec.size)
            if np.any(relvec[mrg15: mrg13] == DELET):
                continue
            mrg25 = max(ins2_pos5 - margin, 0)
            mrg23 = min(ins2_pos5 + margin, relvec.size)
            if np.any(relvec[mrg25: mrg23] == DELET):
                continue
            # Yield the relation vector with these insertions.
            relvec_ins = relvec.copy()
            if ins1_pos5 >= 1:
                relvec_ins[ins1_pos5 - 1] |= INS_5
            if ins1_pos5 < relvec.size:
                relvec_ins[ins1_pos5] |= INS_3
            if ins2_pos5 >= 1:
                relvec_ins[ins2_pos5 - 1] |= INS_5
            if ins2_pos5 < relvec.size:
                relvec_ins[ins2_pos5] |= INS_3
            yield relvec_ins


def iter_relvecs_all(ref_seq: DNA):
    """
    For a given reference sequence, yield every possible unambiguous
    relation vector that has at most two insertions.

    Parameters
    ----------
    ref_seq: DNA
        Sequence of the reference.
    """
    # Use every possible pair of 5' and 3' end positions.
    for end5, end3 in cwr(range(1, len(ref_seq) + 1), 2):
        # Use every possible number of low-quality base calls.
        for n_low_qual in range((end3 - end5 + 1) + 1):
            # Use every possible set of low-quality positions.
            for low_qual in combinations(range(end5, end3 + 1), n_low_qual):
                # Yield every possible relation vector.
                yield from iter_relvecs_q53(ref_seq, low_qual, end5, end3)


def relvec_to_read(ref: DNA, relvec: np.ndarray, hi_qual: int, lo_qual: int):
    """ Make a read sequence from a reference and relation vector. """
    # Validate the relation vector and reference sequence.
    relvec = validate_relvec(relvec)
    if len(ref) != len(relvec):
        raise ValueError(f"Reference sequence has {len(ref)} nt, "
                         f"but relation vector has {len(relvec)} nt")
    # Build the read sequence and quality scores one position at a time.
    read: list[int] = list()
    qual: list[int] = list()
    end5 = 0
    end3 = relvec.size
    ins3_pending = False
    ins3_require = False
    for pos, (base, rel) in enumerate(zip(ref, relvec), start=1):
        if rel == NOCOV:
            # Skip positions not covered by the read.
            if end5:
                # The 5' end was already found, so the 3' end of the
                # read has been reached. Set the 3' end to the previous
                # position (i.e. the last one with coverage) and stop.
                end3 = pos - 1
                break
            # Otherwise, skip to the next position.
            continue
        if end5 == 0:
            # Set the 5' end to the first covered position.
            end5 = pos
        if rel & INS_B:
            # Specially handle insertions because they may overlap any
            # other relation except for a deletion.
            if rel & DELET:
                # Being both a deletion and an insertion is forbidden.
                raise ValueError(f"Relation {rel} in {relvec} is del and ins")
            if rel & INS_3:
                # Being 3' of an insertion is allowed iff ins3_require
                # is True or this is the first covered position.
                if not ins3_require and pos != end5:
                    raise ValueError(f"Unexpected 3' ins at {pos} in {relvec}")
                # The current position is 5' of an insertion, so the
                # inserted base must be added before the base at the
                # current position is added. Add the inserted base.
                read.append(N_INT)
                qual.append(hi_qual)
                # Being 3' of an insertion is not allowed until the next
                # position 5' of an insertion is reached.
                ins3_require = False
            if rel & INS_5:
                # The current position is 5' of an insertion, so the
                # inserted base must be added after the base at the
                # current position is added. Unless the 5' insertion
                # lies at the very 3' end of the read, a 3' insertion
                # will follow it, so typically defer the responsibility
                # of adding the inserted base to the 3' insertion.
                # However, in case the 5' insertion is at the 3' end of
                # the read, flag that a base must be inserted after the
                # current position.
                ins3_pending = True
            # Switch off the insertion flag(s) in the relation so that
            # the base at this position can be added below.
            rel = rel & ~INS_B
        # Check if this position should be 3' of an insertion.
        if ins3_require:
            raise ValueError(f"Required 3' ins not found at {pos} in {relvec}")
        # Update whether the next position must be 3' of an insertion.
        if ins3_pending:
            ins3_require = True
            ins3_pending = False
        if rel == MATCH:
            # Match: Add the reference base to the read.
            read.append(base)
            qual.append(hi_qual)
        elif rel ^ ANY_N in (SUB_A, SUB_C, SUB_G, SUB_T):
            # Ambiguous Sub: Add any nucleotide as low quality.
            read.append(N_INT)
            qual.append(lo_qual)
        elif rel == SUB_T:
            # Substitution to T.
            if base == T_INT:
                raise ValueError(f"Cannot substitute T to itself in {relvec}")
            read.append(T_INT)
            qual.append(hi_qual)
        elif rel == SUB_G:
            # Substitution to G.
            if base == G_INT:
                raise ValueError(f"Cannot substitute G to itself in {relvec}")
            read.append(G_INT)
            qual.append(hi_qual)
        elif rel == SUB_C:
            # Substitution to C.
            if base == C_INT:
                raise ValueError(f"Cannot substitute C to itself in {relvec}")
            read.append(C_INT)
            qual.append(hi_qual)
        elif rel == SUB_A:
            # Substitution to A.
            if base == A_INT:
                raise ValueError(f"Cannot substitute A to itself in {relvec}")
            read.append(A_INT)
            qual.append(hi_qual)
        elif rel != DELET:
            raise ValueError(f"Invalid relation {rel} in {relvec}")
    # Check that the 5' end was found.
    if end5 == 0:
        raise ValueError(f"Relation vector had no 5' end: {relvec}")
    if ins3_require:
        # One base must be inserted at the 3' end because a previous
        # insertion was not handled.
        read.append(N_INT)
        qual.append(hi_qual)
    if len(read) != len(qual):
        raise ValueError(
            f"Lengths of read ({len(read)}) and qual ({len(qual)}) differed")
    if not read:
        raise ValueError("Read contained no bases")
    # Assemble and return the read and quality sequences.
    return bytes(read), bytes(qual), end5, end3


def reads_to_relvecs(ref: DNA):
    """ Return a map from each possible read to the (possibly ambiguous)
    relation vector for the read. """
    # Initialize an empty map.
    read_to_relvec = defaultdict(lambda: np.zeros(len(DNA), dtype=np.uint8))
    # Iterate through all possible relation vectors.
    for relvec in iter_relvecs_all(ref):
        # Determine the read(s) corresponding to this relation vector.
        for read in expand_degenerate_seq(relvec_to_read(ref, relvec, 1, 0)[0]):
            # Accumulate the bitwise OR of all relation vectors.
            read_to_relvec[read] |= relvec
    return read_to_relvec
