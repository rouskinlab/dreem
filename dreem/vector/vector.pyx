# cython: profile=False

import re
from typing import Generator

from ..util.seq import (BLANK, MATCH, DELET, INS_5, INS_3,
                        SUB_A, SUB_C, SUB_G, SUB_T, SUB_N,
                        A_INT, C_INT, G_INT, T_INT, ANY_N)


# C types for sequence and mutation vector characters
cdef unsigned char BLANK_CHR = BLANK
cdef unsigned char MATCH_CHR = MATCH
cdef unsigned char DELET_CHR = DELET
cdef unsigned char INS_5_CHR = INS_5
cdef unsigned char INS_3_CHR = INS_3
cdef unsigned char SUB_A_CHR = SUB_A
cdef unsigned char SUB_C_CHR = SUB_C
cdef unsigned char SUB_G_CHR = SUB_G
cdef unsigned char SUB_T_CHR = SUB_T
cdef unsigned char SUB_N_CHR = SUB_N
cdef unsigned char ANY_N_CHR = ANY_N
cdef unsigned char A_CHR = A_INT
cdef unsigned char C_CHR = C_INT
cdef unsigned char G_CHR = G_INT
cdef unsigned char T_CHR = T_INT

# CIGAR string operation codes
CIG_ALIGN = b"M"[0]  # alignment match
CIG_MATCH = b"="[0]  # sequence match
CIG_SUBST = b"X"[0]  # substitution
CIG_DELET = b"D"[0]  # deletion
CIG_INSRT = b"I"[0]  # insertion
CIG_SCLIP = b"S"[0]  # soft clipping

# C types for CIGAR string operation codes
cdef char CIG_A_CHR = CIG_ALIGN
cdef char CIG_M_CHR = CIG_MATCH
cdef char CIG_S_CHR = CIG_SUBST
cdef char CIG_D_CHR = CIG_DELET
cdef char CIG_I_CHR = CIG_INSRT
cdef char CIG_C_CHR = CIG_SCLIP

# Regular expression pattern that matches a single CIGAR operation
# (length ≥ 1 and operation code, defined above)
CIG_PATTERN = re.compile(b"".join([rb"(\d+)([",
                                   bytearray([CIG_ALIGN,
                                              CIG_MATCH,
                                              CIG_SUBST,
                                              CIG_DELET,
                                              CIG_INSRT,
                                              CIG_SCLIP]),
                                   b"])"]))


class VectorError(Exception):
    """ Any error that occurs during vectoring """


class VectorValueError(VectorError, ValueError):
    """ Any ValueError that occurs during vectoring """


class VectorNotImplementedError(VectorError, NotImplementedError):
    """ Any NotImplementedError that occurs during vectoring """


cdef inline unsigned char encode_base(unsigned char base):
    if base == T_CHR:
        return SUB_T_CHR
    if base == G_CHR:
        return SUB_G_CHR
    if base == C_CHR:
        return SUB_C_CHR
    if base == A_CHR:
        return SUB_A_CHR
    return SUB_N_CHR


cdef inline unsigned char encode_compare(unsigned char ref_base,
                                         unsigned char read_base,
                                         unsigned char read_qual,
                                         unsigned char min_qual):
    return ((MATCH_CHR if ref_base == read_base
             else encode_base(read_base)) if read_qual >= min_qual
            else ANY_N_CHR ^ encode_base(ref_base))


cdef inline unsigned char encode_match(unsigned char read_base,
                                       unsigned char read_qual,
                                       unsigned char min_qual):
    """
    A more efficient version of encode_compare given the prior knowledge
    from the CIGAR string that the read and reference match at this
    position. Note that there is no analagous version when there is a
    known substitution because substitutions are relatively infrequent,
    so optimizing their processing would speed the program only sligtly
    while making the source code more complex and harder to maintain.
    :param read_base:
    :param read_qual:
    :param min_qual:
    :return:
    """
    return (MATCH_CHR if read_qual >= min_qual
            else ANY_N_CHR ^ encode_base(read_base))


class Indel(object):
    """
    Base class for an Insertion or Deletion (collectively, "indel")

    It is used to find alternative positions for indels by keeping track
    of an indel's current coordinates (as it is moved) and determining
    whether a specific move is valid.

    Arguments
    ---------
    rel_ins_idx: int
        The 0-indexed position of the indel with respect to the sequence
        (ref or read) with the relative insertion. This position points
        to one specific base. If the mutation is labeled an insertion,
        then the read is the sequence with the relative insertion (since
        it has a base that is not in the reference), and rel_ins_idx is
        the 0-based index of the inserted base in the coordinates of the
        read sequence. If the mutation is labeled a deletion, then the
        reference is the sequence with the relative insertion (since it
        has a base that is not in the read), and rel_ins_idx is the
        0-based index of the deleted base in the coordinates of the
        reference sequence.
    rel_del_idx (int):  The opposite of rel_ins_idx: the 0-indexed position
                        of the indel with respect to the sequence with a
                        relative deletion (that is, the read if the mutation
                        is denoted a deletion, and the ref if an insertion).
                        Because the deleted base does not actually exist
                        within the sequence whose coordinates it is based on,
                        rel_del_idx does not refer to a specific position in
                        the sequence, but instead to the two extant positions
                        in the sequence that flank the deleted position.
                        It is most convenient for the algorithm to have this
                        argument refer to the position 3' of the deleted base
                        and define the 5' position as a property.
    """

    # Define __slots__ to improve speed and memory performance.
    __slots__ = ["_ins_idx", "_ins_init", "_del_idx", "_del_init", "_tunneled"]

    # Minimum distance between an insertion and a deletion
    MIN_INDEL_DIST = 2

    def __init__(self, rel_ins_idx: int, rel_del_idx: int) -> None:
        self._ins_idx = rel_ins_idx
        self._ins_init = rel_ins_idx
        self._del_idx = rel_del_idx
        self._del_init = rel_del_idx
        self._tunneled = False

    @property
    def ins_idx(self):
        return self._ins_idx

    @property
    def del_idx5(self):
        return self._del_idx - 1

    @property
    def del_idx3(self):
        return self._del_idx

    @property
    def tunneled(self):
        return self._tunneled

    @property
    def rank(self) -> int:
        raise VectorNotImplementedError

    def reset(self):
        """ Reset the indel to its initial position, and erase its
        history of tunneling. """
        self._ins_idx = self._ins_init
        self._del_idx = self._del_init
        self._tunneled = False

    @staticmethod
    def _get_indel_by_idx(indels: list, idx: int):
        for indel in indels:
            if indel.ins_idx == idx:
                return indel

    def _peek_out_of_indel(self, indels: list, from3to5: bool):
        inc = -1 if from3to5 else 1
        idx = self.ins_idx + inc
        tunneled_indels: list = list()
        while True:
            indel = self._get_indel_by_idx(indels, idx)
            if not indel:
                break
            idx += inc
            tunneled_indels.append(indel)
        self._tunneled = bool(tunneled_indels)
        return idx, tunneled_indels

    def _collision(self, other: "Indel", swap_idx: int):
        return self.MIN_INDEL_DIST > (min(abs(swap_idx - other.del_idx5),
                                          abs(swap_idx - other.del_idx3)))

    def _collisions(self, indels: list, swap_idx: int):
        return any(self._collision(indel, swap_idx) for indel in indels)

    def step_del_idx(self, swap_idx: int):
        # Move the indel's position (self._ins_idx) to swap_idx.
        # Move self._del_idx one step in the same direction.
        if swap_idx == self.ins_idx:
            raise VectorValueError(f"swap ({swap_idx}) = ins ({self.ins_idx})")
        self._del_idx += 1 if swap_idx > self.ins_idx else -1

    def _step(self, swap_idx: int):
        self.step_del_idx(swap_idx)
        self._ins_idx = swap_idx

    @staticmethod
    def _consistent_rels(curr_rel: int, swap_rel: int):
        if curr_rel & swap_rel or (curr_rel & SUB_N_CHR
                                   and swap_rel & SUB_N_CHR):
            # Relationship between reference and read base (read_code) and
            # relationship between reference and swap base (swap_code)
            # are consistent, meaning either
            # - both match the reference
            # - one matches and the other potentially matches (i.e. low qual)
            # - one is a substitution and the other could be a substitution
            # - both are substitutions (for each code, code & SUB_N_INT == code)
            return curr_rel
        # Otherwise, i.e.g if one base matches and the other is a substitution,
        # then the relationships are not consistent.
        return 0

    def _encode_swap(self, *args, **kwargs) -> bool:
        raise VectorNotImplementedError

    def _try_swap(self, *args, **kwargs) -> bool:
        raise VectorNotImplementedError

    def sweep(self, muts: bytearray, ref: bytes, read: bytes, qual: bytes,
              min_qual: int, dels: list, inns: list,
              from3to5: bool, tunnel: bool):
        # Move the indel as far as possible in either the 5' or 3' direction.
        while self._try_swap(muts, ref, read, qual, min_qual, dels, inns,
                             from3to5, tunnel):
            # All actions happen in _try_swap, so loop body is empty.
            pass


class Deletion(Indel):
    @property
    def rank(self):
        return self._ins_idx

    @classmethod
    def _encode_swap(cls, ref_base: int, swap_base: int, read_base: int,
                     read_qual: int, min_qual: int):
        curr_rel = encode_compare(ref_base, read_base, read_qual, min_qual)
        swap_rel = encode_compare(swap_base, read_base, read_qual, min_qual)
        return cls._consistent_rels(curr_rel, swap_rel)

    def _swap(self, muts: bytearray, swap_idx: int, relation: int):
        """
        Arguments
        ---------
        muts: bytearray
            Mutation vector
        swap_idx: int
            Index in the region to which the deletion moves during this
            swap
        relation: int
            Relationship (match, sub, etc.) between the base located at
            swap_idx and the base in the read
        """
        # The base at swap_idx moves to self.ref_idx, so after the swap, the
        # relationship between self.ref_idx and the read base will be swap_code.
        muts[self.ins_idx] = muts[self.ins_idx] | relation
        # The base at self.ref_idx is marked as a deletion (by definition), so
        # mark the position it moves to (swap_idx) as a deletion too.
        muts[swap_idx] = muts[swap_idx] | DELET_CHR
        self._step(swap_idx)

    def _try_swap(self, muts: bytearray, ref: bytes, read: bytes, qual: bytes,
                  min_qual: int, dels: list, inns: list,
                  from3to5: bool, tunnel: bool) -> bool:
        swap_idx, tunneled_indels = self._peek_out_of_indel(dels, from3to5)
        read_idx = self.del_idx5 if from3to5 else self.del_idx3
        if (0 <= swap_idx < len(ref) and 0 <= read_idx < len(read)
                and (tunnel or not self.tunneled)
                and not self._collisions(inns, swap_idx)):
            relation = self._encode_swap(ref[self.ins_idx], ref[swap_idx],
                                         read[read_idx], qual[read_idx],
                                         min_qual)
            if relation:
                self._swap(muts, swap_idx, relation)
                for indel in tunneled_indels:
                    indel.step_del_idx(swap_idx)
                return True
        return False


class Insertion(Indel):
    @property
    def rank(self):
        return self._del_idx

    def stamp(self, muts: bytearray):
        if 0 <= self.del_idx5 < len(muts):
            muts[self.del_idx5] = muts[self.del_idx5] | INS_5_CHR
        if 0 <= self.del_idx3 < len(muts):
            muts[self.del_idx3] = muts[self.del_idx3] | INS_3_CHR

    @classmethod
    def _encode_swap(cls, ref_base: int, read_base: int, read_qual: int,
                     swap_base: int, swap_qual: int, min_qual: int):
        curr_rel = encode_compare(ref_base, read_base, read_qual, min_qual)
        swap_rel = encode_compare(ref_base, swap_base, swap_qual, min_qual)
        return cls._consistent_rels(curr_rel, swap_rel)

    def _swap(self, muts: bytearray, ref_idx: int,
              swap_idx: int, relation: int):
        """
        Arguments
        mut_vectors (bytearray): mutation vector
        swap_idx (int): the index in the read to which the deletion moves
                        during this swap
        swap_code (int): the relationship (match, sub, etc.) between the
                         base located at swap_idx and the base in the ref
        """
        # The base at ins_idx moves to swap_idx, so after the swap, the
        # relationship between ref_idx and the read base will be relation.
        muts[ref_idx] = muts[ref_idx] | relation
        self._step(swap_idx)
        # Mark the new positions of the insertion.
        self.stamp(muts)

    def _try_swap(self, muts: bytearray, ref: bytes, read: bytes, qual: bytes,
                  min_qual: int, dels: list, inns: list,
                  from3to5: bool, tunnel: bool) -> bool:
        swap_idx, tunneled_indels = self._peek_out_of_indel(inns, from3to5)
        ref_idx = self.del_idx5 if from3to5 else self.del_idx3
        if (0 <= swap_idx < len(read) and 0 <= ref_idx < len(ref)
                and (tunnel or not self.tunneled)
                and not self._collisions(dels, swap_idx)):
            relation = self._encode_swap(ref[ref_idx], read[self.ins_idx],
                                         qual[self.ins_idx], read[swap_idx],
                                         qual[swap_idx], min_qual)
            if relation:
                self._swap(muts, ref_idx, swap_idx, relation)
                for indel in tunneled_indels:
                    indel.step_del_idx(swap_idx)
                return True
        return False


def sweep_indels(muts: bytearray, ref: bytes, read: bytes, qual: bytes,
                 min_qual: int, dels: list, inns: list,
                 from3to5: bool, tunnel: bool):
    """
    For every insertion and deletion,

    Parameters
    ----------
    muts: bytearray
        Mutation vector
    ref: bytes
        Sequence of the region of the reference that the mutation vector
        covers
    read: bytes
        Sequence of the read
    qual: bytes
        Phred quality scores of the read, encoded as ASCII characters
    min_qual: int
        The minimum Phred quality score needed to consider a base call
        informative: integer value of the ASCII character
    dels: list
        Deletions identified by ```vectorize_read```
    inns: list
        Insertions identified by ```vectorize_read```
    from3to5: bool
        Whether to move indels in the 3' -> 5' direction (True) or the
        5' -> 3' direction (False)
    tunnel: bool
        Whether to allow tunneling
    """
    # Collect all indels into one list.
    indels: list = list()
    indels.extend(dels)
    indels.extend(inns)
    # Reset each indel to its initial state. This operation does nothing
    # the first time sweep_indels is called because all indels start in
    # their initial state (by definition). But the indels may move when
    # this function runs, so resetting is necessary at the beginning of
    # the second and subsequent calls to sweep_indels to ensure that the
    # algorithm starts from the initial state every time.
    for indel in indels:
        indel.reset()
    # Sort the indels by their rank, which is
    sort_rev = from3to5 != tunnel
    indels.sort(key=lambda idl: idl.rank, reverse=sort_rev)
    while indels:
        indel = indels.pop()
        indel.sweep(muts, ref, read, qual, min_qual,
                    dels, inns, from3to5, tunnel)
        i = len(indels)
        if sort_rev:
            while i > 0 and indel.rank > indels[i - 1].rank:
                i -= 1
        else:
            while i > 0 and indel.rank < indels[i - 1].rank:
                i -= 1
        if i < len(indels):
            indels.insert(i, indel)


def get_ambids(muts: bytearray, ref: bytes, read: bytes, qual: bytes,
               min_qual: int, dels: list, inns: list):
    """
    Find and label all positions in the vector that are ambiguous due to
    insertions and deletions.

    Parameters
    ----------
    muts: bytearray
        Mutation vector
    ref: bytes
        Sequence of the region of the reference that the mutation vector
        covers
    read: bytes
        Sequence of the read
    qual: bytes
        Phred quality scores of the read, encoded as ASCII characters
    min_qual: int
        The minimum Phred quality score needed to consider a base call
        informative: integer value of the ASCII character
    dels: list
        Deletions identified by ```vectorize_read```
    inns: list
        Insertions identified by ```vectorize_read```
    """
    # Each indel might be able to be moved in the 5' -> 3' direction
    # (from3to5 is False) or 3' -> 5' direction (from3to5 is True).
    # Test both directions.
    for from3to5 in (False, True):
        # For each indel, try to move it as far as it can go in the
        # direction indicated by from3to5. Allow tunneling so that any
        # runs of consecutive insertions or consecutive deletions can
        # effectively move together.
        sweep_indels(muts, ref, read, qual, min_qual,
                     dels, inns, from3to5, tunnel=True)
        if any(d.tunneled for d in dels) or any(i.tunneled for i in inns):
            # If any indel tunneled,
            sweep_indels(muts, ref, read, qual, min_qual,
                         dels, inns, from3to5, tunnel=False)


def parse_cigar(cigar_string: bytes) -> Generator[tuple[int, int], None, None]:
    """
    Yield the fields of a CIGAR string as pairs of (operation, length),
    where operation is 1 byte indicating the CIGAR operation and length
    is a positive integer indicating the number of bases from the read
    that the operation consumes. Note that in the CIGAR string itself,
    each length precedes its corresponding operation.

    Parameters
    ----------
    cigar_string: bytes
        CIGAR string from a SAM file. For full documentation, refer to
        https://samtools.github.io/hts-specs/

    Yield
    -----
    bytes (length = 1)
        Current CIGAR operation
    int (≥ 1)
        Length of current CIGAR operation

    Examples
    --------
    >>> list(parse_cigar(b"17=1X43=1D26="))
    [(b'=', 17), (b'X', 1), (b'=', 43), (b'D', 1), (b'=', 26)]
    """
    # Length-0 CIGAR strings are forbidden.
    if not cigar_string:
        raise VectorValueError("CIGAR string is empty")
    # If the CIGAR string has any invalid bytes (e.g. an unrecognized
    # operation byte, an operation longer than 1 byte, a length that is
    # not a positive integer, or any extraneous characters), then the
    # regular expression parser will simply skip these invalid bytes.
    # In order to catch such problems, keep track of the number of
    # bytes matched from the CIGAR string. After reading the CIGAR, if
    # the number of bytes matched is smaller than the length of the
    # CIGAR string, then some bytes must have been skipped, indicating
    # that the CIGAR string contained at least one invalid byte.
    num_bytes_matched = 0
    # Find every operation in the CIGAR string that matches the regular
    # expression.
    for match in CIG_PATTERN.finditer(cigar_string):
        length_bytes, operation = match.groups()
        # Convert the length field from bytes to int.
        length_int = int(length_bytes)
        # Add the total number of bytes in the current operation to the
        # count of the number of bytes matched from the CIGAR string.
        num_bytes_matched += len(length_bytes) + len(operation)
        # Note that the fields are yielded as (operation, length), but
        # in the CIGAR string itself, the order is (length, operation).
        yield operation[0], length_int
    # Confirm that all bytes in the CIGAR string were matched by the
    # regular expression. Note: This check will only be performed if
    # the entire CIGAR string is read. Thus, it is essential to read
    # the entire CIGAR string, even if the read extends beyond the
    # region for which the mutation vector is being computed.
    if num_bytes_matched != len(cigar_string):
        raise VectorValueError(f"Invalid CIGAR: '{cigar_string.decode()}'")


cdef inline bint op_consumes_ref(char op):
    """ Return whether the CIGAR operation consumes the reference. """
    return op != CIG_I_CHR and op != CIG_C_CHR


cdef inline bint op_consumes_read(char op):
    """ Return whether the CIGAR operation consumes the read. """
    return op != CIG_D_CHR


class SamFlag(object):
    """ Represents the set of 12 boolean flags for a SAM record. """

    # Define __slots__ to improve speed and memory performance.
    __slots__ = ["paired", "proper", "unmap", "munmap", "rev", "mrev",
                 "first", "second", "secondary", "qcfail", "dup", "supp"]

    # Maximum value of a valid SAM flag representation, corresponding
    # to all 12 flags set to 1: 111111111111 (binary) = 4095 (decimal)
    MAX_FLAG: int = 2 ** len(__slots__) - 1
    # Pattern for padding the left of the binary string with 0s.
    PATTERN = "".join(["{:0>", str(len(__slots__)), "}"])

    def __init__(self, flag: int):
        """
        Validate the integer value of the SAM flag, then set the 12
        individual flag values. To maximize speed, all flags are set in
        a single one-line operation, each step of which is explained:

        1.  Convert the flag (int) to a binary representation (str) that
            starts with '0b':
            >>> flag_int = 83
            >>> flag_bin = bin(flag_int)
            >>> flag_bin
            '0b1010011'

        2.  Remove the prefix '0b':
            >>> flag_bits = flag_bin[2:]
            >>> flag_bits
            '1010011'

        3.  Pad the left side of the string with 0 up to a length of 12:
            >>> PATTERN = "".join(["{:0>12}"])
            >>> PATTERN
            '{:0>12}'
            >>> all_bits = PATTERN.format(flag_bits)
            '000001010011'

        4.  Convert '1' to True and '0' to False, and assign to the 12
            flag bits in order from greatest (supp: 2048) to least
            (paired: 1) numerical value of the flag bit:
            >>> (supp, dup, qcfail, secondary, second, first, mrev, rev,
            ...  munmap, unmap, proper, paired) = map(bool, map(int, all_bits))
            >>> paired, rev, second, qcfail
            (True, True, False, False)

        Parameters
        ----------
        flag: int
            The integer value of the SAM flag. For documentation, see
            https://samtools.github.io/hts-specs/

        Examples
        --------
        >>> flag0099 = SamFlag(99)
        >>> flag0099.paired, flag0099.rev
        (True, False)
        """
        if not 0 <= flag <= self.MAX_FLAG:
            raise VectorValueError(f"Invalid flag: '{flag}'")
        (self.supp, self.dup, self.qcfail, self.secondary,
         self.second, self.first, self.mrev, self.rev,
         self.munmap, self.unmap, self.proper, self.paired) = (
            x == '1' for x in self.PATTERN.format(bin(flag)[2:])
        )


class SamRead(object):
    # Define __slots__ to improve speed and memory performance.
    __slots__ = ["qname", "flag", "rname", "pos", "mapq", "cigar",
                 "tlen", "seq", "qual", "min_qual"]

    # Minimum number of fields in a valid SAM record
    MIN_FIELDS = 11

    def __init__(self, line: bytes):
        fields = line.rstrip().split(b"\t")
        if len(fields) < self.MIN_FIELDS:
            raise VectorValueError(f"Invalid SAM line:\n{line}")
        self.qname = fields[0]
        self.flag = SamFlag(int(fields[1]))
        self.rname = fields[2]
        self.pos = int(fields[3])
        self.mapq = int(fields[4])
        self.cigar = fields[5]
        # RNEXT, PNEXT, and TLEN are not used during vectoring,
        # so they are commented out to reduce processing time
        # but still shown here in case they are ever needed.
        # self.rnext = fields[6]
        # self.pnext = int(fields[7])
        # self.tlen = int(fields[8])
        self.seq = fields[9]
        self.qual = fields[10]
        if len(self.seq) != len(self.qual):
            raise VectorValueError(f"Lengths of seq ({len(self.seq)}) and qual "
                                   f"string {len(self.qual)} did not match.")


cdef inline void validate_op_end3_read(int op_end3_read, int read_len):
    """ Validate that the position of the 3' end of the read has not
    exceeded the length of the read. """
    if op_end3_read > read_len:
        raise VectorValueError(f"Operation ending at {op_end3_read} "
                               f"overshot the read (length = {read_len})")
    


cdef inline void skip_extra5(int* op_end5_region,
                             int* op_end5_read,
                             int* op_len):
    """ Skip extra positions (if any) at the 5' end of the operation,
    before the 5' end of the region. """
    if op_end5_region[0] < 0:
        # If the operation starts 5' of the region (i.e. the 5' position
        # of the operation with respect to the region is negative), then
        # decrease the length of the operation,
        op_len[0] += op_end5_region[0]
        # advance the 5' end of the operation with respect to the read,
        op_end5_read[0] -= op_end5_region[0]
        # and advance the 5' end of the operation with respect to the
        # region to exactly 5' end of the region (i.e. position 0).
        op_end5_region[0] = 0


cdef _vectorize_read(int region_end5_ref,
                     int region_len,
                     unsigned char* region_seq,
                     bytes cigar,
                     int read_end5_ref,
                     int read_len,
                     unsigned char* read_seq,
                     unsigned char* read_qual,
                     unsigned char min_qual,
                     bint ambid):
    """
    Generate and return a mutation vector of an aligned read over a
    given region of the reference sequence.

    Parameters
    ----------
    read: SamRead
        Read from SAM file to be vectorized
    region_seq: bytes
        Sequence of the region for which to compute the mutation vector
        (only the region, not any other part of the reference sequence)
    region_end5: int (≥ 1)
        5'-most coordinate of the region with respect to the entire
        reference sequence (1-indexed, includes coordinate)
    region_end3: int (≥ region_end5)
        3'-most coordinate of the region with respect to the entire
        reference sequence (1-indexed, includes coordinate)
    min_qual: int
        ASCII encoding of the minimum Phred score to accept a base call
    ambid: bool
        Whether to find and label all ambiguous insertions and deletions

    Return
    ------
    bytearray
        Mutation vector, whose length either equals that of the region
        or is zero to indicate an error occurred during vectorization.
    """
    # Indexes of the 5' and 3' ends of the current CIGAR operation with
    # respect to the read; 0-indexed, uses Python half-open intervals
    cdef int op_end5_read = 0
    cdef int op_end3_read = op_end5_read
    # Indexes of the 5' and 3' ends of the current CIGAR operation with
    # respect to the region; 0-indexed, uses Python half-open intervals,
    # and is negative if the operation is upstream of the region 5' end
    cdef int op_end5_region = read_end5_ref - region_end5_ref
    cdef int op_end3_region = op_end5_region
    # Initialize a blank mutation vector covering the entire region.
    muts = bytearray([BLANK] * region_len)
    # Initialize lists to record all deletions and insertions.
    dels = list()
    inns = list()
    # Read each operation in the CIGAR string.
    cdef unsigned char cigar_op
    cdef int op_len
    for cigar_op, op_len in parse_cigar(cigar):
        # Act based on the CIGAR operation and its length.
        if cigar_op == CIG_M_CHR:  # match
            # Update the position of the current operation's 3' end.
            op_end3_region += op_len
            op_end3_read += op_len
            validate_op_end3_read(op_end3_read, read_len)
            # If at least one position of the operation overlaps the
            # region at all, then compute bytes for this operation.
            if op_end3_region > 0 and op_end5_region < region_len:
                # Skip any extra positions after the 5' end of the
                # operation and before the 5' end of the region.
                skip_extra5(&op_end5_region, &op_end5_read, &op_len)
                # Find the position (with respect to the region) of the
                # furthest 3' byte to compute in the current operation.
                # Usually this position is the 3' end of the operation,
                # but also it cannot exceed the length of the region.
                op_end5_region_limit3 = min(op_end3_region, region_len)
                # Visit every position in the operation, stopping at the
                # limit on the 3' side.
                while op_end5_region < op_end5_region_limit3:
                    # Compute the byte of the mutation vector at the
                    # current position in the region.
                    muts[op_end5_region] = encode_match(read_seq[op_end5_read],
                                                        read_qual[op_end5_read],
                                                        min_qual)
                    # Advance one position in the region and read.
                    op_end5_region += 1
                    op_end5_read += 1
        elif cigar_op == CIG_A_CHR or cigar_op == CIG_S_CHR:  # no indel
            # Update the position of the current operation's 3' end.
            op_end3_region += op_len
            op_end3_read += op_len
            validate_op_end3_read(op_end3_read, read_len)
            # If at least one position of the operation overlaps the
            # region at all, then compute bytes for this operation.
            if op_end3_region > 0 and op_end5_region < region_len:
                # Skip any extra positions after the 5' end of the
                # operation and before the 5' end of the region.
                skip_extra5(&op_end5_region, &op_end5_read, &op_len)
                # Find the position (with respect to the region) of the
                # furthest 3' byte to compute in the current operation.
                # Usually this position is the 3' end of the operation,
                # but also it cannot exceed the length of the region.
                op_end5_region_limit3 = min(op_end3_region, region_len)
                # Visit every position in the operation, stopping at the
                # limit on the 3' side.
                while op_end5_region < op_end5_region_limit3:
                    # Compute the byte of the mutation vector at the
                    # current position in the region.
                    muts[op_end5_region] = encode_compare(region_seq[
                                                              op_end5_region],
                                                          read_seq[
                                                              op_end5_read],
                                                          read_qual[
                                                              op_end5_read],
                                                          min_qual)
                    # Advance one position in the region and read.
                    op_end5_region += 1
                    op_end5_read += 1
        elif cigar_op == CIG_D_CHR:  # deletion from read
            # Update the position of the current operation's 3' end.
            op_end3_region += op_len
            # If at least one position of the operation overlaps the
            # region at all, then compute bytes for this operation.
            if op_end3_region > 0 and op_end5_region < region_len:
                # Skip any extra positions after the 5' end of the
                # operation and before the 5' end of the region.
                skip_extra5(&op_end5_region, &op_end5_read, &op_len)
                # Find the position (with respect to the region) of the
                # furthest 3' byte to compute in the current operation.
                # Usually this position is the 3' end of the operation,
                # but also it cannot exceed the length of the region.
                op_end5_region_limit3 = min(op_end3_region, region_len)
                # Visit every position in the operation, stopping at the
                # limit on the 3' side.
                while op_end5_region < op_end5_region_limit3:
                    # Put a deletion byte into the mutation vector.
                    muts[op_end5_region] = DELET_CHR
                    # Create a deletion object (used by get_ambids).
                    dels.append(Deletion(op_end5_region, op_end5_read))
                    # Advance one position in the region.
                    op_end5_region += 1
        elif cigar_op == CIG_I_CHR:  # insertion into read
            # Update the position of the current operation's 3' end.
            op_end3_read += op_len
            validate_op_end3_read(op_end3_read, read_len)
            # If at least one position of the operation overlaps the
            # region at all, then compute bytes for this operation.
            if op_end3_region > 0 and op_end5_region < region_len:
                # Visit every position in the operation, stopping at the
                # limit on the 3' side.
                while op_end5_read < op_end3_read:
                    # Create an insertion object (used by get_ambids).
                    # Every mutation needs to be assigned a coordinate
                    # in the region, which is the coordinate at which it
                    # appears in the mutation vector. But each inserted
                    # base, being absent from the reference, does not
                    # correspond to a single coordinate in the region.
                    # Instead, each inserted base lies between two: the
                    # coordinates of the region immediately 5' and 3' of
                    # the insertion. Either could be designated as the
                    # coordinate of the insertion within the region.
                    # This code uses the 3' coordinate: for example, a
                    # base inserted between coordinates 45 and 46 of the
                    # region would be assigned coordinate 46. Doing so
                    # simplifies the math, as the 3' coordinate already
                    # equals both op_end5_region and op_end3_region.
                    inns.append(Insertion(op_end5_read, op_end5_region))
                    # Advance one position in the read.
                    op_end5_read += 1
                    # Insertions do not consume the reference, so add
                    # no information to muts yet; it is added by stamp.
        elif cigar_op == CIG_C_CHR:  # soft clipping from read
            # Like insertions, soft clippings consume the read but not
            # the reference. Unlike insertions, they are not mutations,
            # so they require no further processing.
            op_end3_read += op_len
            validate_op_end3_read(op_end3_read, read_len)
        else:
            raise VectorValueError(
                f"Invalid CIGAR operation: '{cigar_op.decode()}'")
        # Advance the 5' positions in the region and read to the current
        # 3' positions so that the next CIGAR operation lies immediately
        # 3' of the current operation. The 3' positions will be advanced
        # at the beginning of the next iteration (if any) of the loop.
        op_end5_region = op_end3_region
        op_end5_read = op_end3_read
    # Verify that the sum of all CIGAR operations that consumed the read
    # equals the length of the read.
    if op_end5_read != read_len:
        raise VectorValueError(
            f"CIGAR string '{cigar.decode()}' consumed {op_end5_read} bases "
            f"from read of length {len(read_len)}.")
    # Add insertions to mut_vectors.
    for ins in inns:
        ins.stamp(muts)
    # Label all positions that are ambiguous due to indels.
    if ambid and (dels or inns):
        get_ambids(muts, region_seq, read_seq, read_qual, min_qual, dels, inns)
    return muts


def vectorize_read(read: SamRead, *,
                   region_seq: bytes,
                   region_end5: int,
                   min_qual: int,
                   ambid: bool):
    return _vectorize_read(region_end5, len(region_seq), region_seq, read.cigar,
                           read.pos, len(read.seq), read.seq, read.qual,
                           min_qual, ambid)


cdef unsigned char get_consensus_mut(unsigned char byte1,
                                     unsigned char byte2):
    if byte1:
        if byte2:
            intersect = byte1 & byte2
            if intersect:
                return intersect
            return byte1 | byte2
        return byte1
    return byte2


def vectorize_pair(read1: SamRead, read2: SamRead, **kwargs):
    return bytearray(map(get_consensus_mut,
                         vectorize_read(read1, **kwargs),
                         vectorize_read(read2, **kwargs)))


class SamRecord(object):
    __slots__ = ["read1", "read2", "strict"]

    def __init__(self,
                 read1: SamRead,
                 read2: SamRead | None = None,
                 strict: bool = True):
        if read2 is None:
            if read1.flag.paired and strict:
                # If the region does not span the reference sequence,
                # then it is possible for read 1 but not read 2 to
                # overlap the region, or vice versa. If this happens,
                # then strict mode is turned off to allow processing
                # the one read that overlaps the region even though its
                # mate (which does not overlap the region) is absent.
                raise VectorValueError(f"Read 1 '{self.read1.qname.decode()}' "
                                       "was paired, but no read 2 was given")
        else:
            if read1.flag.paired:
                if read1.qname != read2.qname:
                    raise VectorValueError(f"Mates 1 '{read1.qname.decode()}' "
                                           f"and 2 '{read2.qname.decode()}') "
                                           "had different read names")
                if read1.rname != read2.rname:
                    raise VectorValueError(f"Read '{read1.qname.decode()}' had "
                                           "different references for mates 1 "
                                           f"('{read1.rname.decode()}') and 2 "
                                           f"('{read2.rname.decode()}')")
                if not read2.flag.paired:
                    raise VectorValueError(f"Read '{read1.qname.decode()}' had "
                                           "paired mate 1 not unpaired mate 2")
                if not (read1.flag.first and read2.flag.second):
                    raise VectorValueError(f"Read '{read1.qname.decode()}' had "
                                           f"mate 1 = {2 - read1.flag.first}, "
                                           f"mate 2 = {1 + read2.flag.second}")
                if read1.flag.rev == read2.flag.rev:
                    raise VectorValueError(f"Read '{read1.qname.decode()}' had "
                                           "mates 1 and 2 facing the same way")
            else:
                raise VectorValueError(f"Mate 1 ('{read1.qname.decode()}') "
                                       "was not paired, but mate 2 "
                                       f"('{read2.qname.decode()}') was given")
        self.read1 = read1
        self.read2 = read2
        self.strict = strict

    def vectorize(self, **kwargs):
        if self.read2 is None:
            return vectorize_read(self.read1, **kwargs)
        else:
            return vectorize_pair(self.read1, self.read2, **kwargs)
