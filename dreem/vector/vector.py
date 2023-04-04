from __future__ import annotations
import re

from ..util.seq import (MATCH, DELET, INS_5, INS_3,
                        SUB_A, SUB_C, SUB_G, SUB_T, SUB_N,
                        A_INT, C_INT, G_INT, T_INT, ANY_N)


class VectorError(Exception):
    """ Any error that occurs during vectoring """


class VectorValueError(VectorError, ValueError):
    """ Any ValueError that occurs during vectoring """


class VectorNotImplementedError(VectorError, NotImplementedError):
    """ Any NotImplementedError that occurs during vectoring """


# CIGAR string operation codes
CIG_ALIGN = b"M"  # alignment match
CIG_MATCH = b"="  # sequence match
CIG_SUBST = b"X"  # substitution
CIG_DELET = b"D"  # deletion
CIG_INSRT = b"I"  # insertion
CIG_SCLIP = b"S"  # soft clipping

# Regular expression pattern that matches a single CIGAR operation
# (length ≥ 1 and operation code, defined above)
CIG_PATTERN = re.compile(b"".join([rb"(\d+)([",
                                   CIG_ALIGN,
                                   CIG_MATCH,
                                   CIG_SUBST,
                                   CIG_DELET,
                                   CIG_INSRT,
                                   CIG_SCLIP,
                                   b"])"]))


def encode_base(base: int):
    if base == T_INT:
        return SUB_T
    if base == G_INT:
        return SUB_G
    if base == C_INT:
        return SUB_C
    if base == A_INT:
        return SUB_A
    raise VectorValueError(f"Invalid base: '{chr(base)}'")


def encode_compare(ref_base: int, read_base: int, read_qual: int, min_qual: int):
    return ((MATCH if ref_base == read_base
             else encode_base(read_base)) if read_qual >= min_qual
            else ANY_N ^ encode_base(ref_base))


def encode_match(read_base: int, read_qual: int, min_qual: int):
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
    return (MATCH if read_qual >= min_qual
            else ANY_N ^ encode_base(read_base))


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
    def _get_indel_by_idx(indels: list[Indel], idx: int):
        for indel in indels:
            if indel.ins_idx == idx:
                return indel

    def _peek_out_of_indel(self, indels: list[Indel], from3to5: bool):
        inc = -1 if from3to5 else 1
        idx = self.ins_idx + inc
        tunneled_indels: list[Indel] = list()
        while indel := (self._get_indel_by_idx(indels, idx)):
            idx += inc
            tunneled_indels.append(indel)
        self._tunneled = bool(tunneled_indels)
        return idx, tunneled_indels

    def _collision(self, other: Indel, swap_idx: int):
        return self.MIN_INDEL_DIST > (min(abs(swap_idx - other.del_idx5),
                                          abs(swap_idx - other.del_idx3)))

    def _collisions(self, indels: list[Indel], swap_idx: int):
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
        if curr_rel & swap_rel or (curr_rel & SUB_N
                                   and swap_rel & SUB_N):
            # Relationship between reference and read base (read_code) and
            # relationship between reference and swap base (swap_code)
            # are consistent, meaning either
            # - both match the reference
            # - one matches and the other potentially matches (i.e. low qual)
            # - one is a substitution and the other could be a substitution
            # - both are substitutions (for each code, code & SUB_N == code)
            return curr_rel
        # Otherwise, i.e.g if one base matches and the other is a substitution,
        # then the relationships are not consistent.
        return 0

    def _encode_swap(self, *args, **kwargs) -> bool:
        raise VectorNotImplementedError

    def _try_swap(self, *args, **kwargs) -> bool:
        raise VectorNotImplementedError

    def sweep(self, muts: bytearray, ref: bytes, read: bytes, qual: bytes,
              min_qual: int, dels: list[Deletion], inns: list[Insertion],
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
        muts[swap_idx] = muts[swap_idx] | DELET
        self._step(swap_idx)

    def _try_swap(self, muts: bytearray, ref: bytes, read: bytes, qual: bytes,
                  min_qual: int, dels: list[Deletion], inns: list[Insertion],
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
            muts[self.del_idx5] = muts[self.del_idx5] | INS_5
        if 0 <= self.del_idx3 < len(muts):
            muts[self.del_idx3] = muts[self.del_idx3] | INS_3

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
                  min_qual: int, dels: list[Deletion], inns: list[Insertion],
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
                 min_qual: int, dels: list[Deletion], inns: list[Insertion],
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
    dels: list[Deletion]
        List of deletions identified by ```vectorize_read```
    inns: list[Insertion]
        List of insertions identified by ```vectorize_read```
    from3to5: bool
        Whether to move indels in the 3' -> 5' direction (True) or the
        5' -> 3' direction (False)
    tunnel: bool
        Whether to allow tunneling
    """
    # Collect all indels into one list.
    indels: list[Indel] = list()
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
               min_qual: int, dels: list[Deletion], inns: list[Insertion]):
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
    dels: list[Deletion]
        List of deletions identified by ```vectorize_read```
    inns: list[Insertion]
        List of insertions identified by ```vectorize_read```
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


def parse_cigar(cigar_string: bytes):
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
        # Convert the length field from bytes to int and verify that it
        # is a positive integer.
        if (length_int := int(length_bytes)) < 1:
            raise VectorValueError("length of CIGAR operation must be ≥ 1")
        # Add the total number of bytes in the current operation to the
        # count of the number of bytes matched from the CIGAR string.
        num_bytes_matched += len(length_bytes) + len(operation)
        # Note that the fields are yielded as (operation, length), but
        # in the CIGAR string itself, the order is (length, operation).
        yield operation, length_int
    # Confirm that all bytes in the CIGAR string were matched by the
    # regular expression. Note: This check will only be performed if
    # the entire CIGAR string is read. Thus, it is essential to read
    # the entire CIGAR string, even if the read extends beyond the
    # region for which the mutation vector is being computed.
    if num_bytes_matched != len(cigar_string):
        raise VectorValueError(f"Invalid CIGAR: '{cigar_string.decode()}'")


def op_consumes_ref(op: bytes) -> bool:
    """ Return whether the CIGAR operation consumes the reference. """
    return op != CIG_INSRT and op != CIG_SCLIP


def op_consumes_read(op: bytes) -> bool:
    """ Return whether the CIGAR operation consumes the read. """
    return op != CIG_DELET


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
            >>> (flag_bin := bin(flag_int := 83))
            '0b1010011'
        2.  Remove the prefix '0b':
            >>> (flag_bits := flag_bin[2:])
            '1010011'
        3.  Pad the left side of the string with 0 up to a length of 12:
            >>> (PATTERN := "".join(["{:0>", str(num_flags := 12), "}"]))
            '{:0>12}'
            >>> (all_bits := PATTERN.format(flag_bits))
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
    __slots__ = ["qname", "flag", "paired", "rname", "pos", "mapq", "cigar",
                 "tlen", "seq", "qual", "min_qual"]

    # Minimum number of fields in a valid SAM record
    MIN_FIELDS = 11

    def __init__(self, line: bytes):
        fields = line.rstrip().split(b"\t")
        if len(fields) < self.MIN_FIELDS:
            raise VectorValueError(f"Invalid SAM line:\n{line}")
        self.qname = fields[0]
        self.flag = SamFlag(int(fields[1]))
        self.paired = self.flag.paired
        self.rname = fields[2].decode()
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


def vectorize_read(read: SamRead,
                   muts: bytearray,
                   region_seq: bytes,
                   region_length: int,
                   region_end5: int,
                   min_qual: int,
                   ambid: bool):
    """
    Generate and return a mutation vector of an aligned read over a
    given region of the reference sequence.
    Parameters
    ----------
    read: SamRead
        Read from SAM file to be vectorized
    muts: bytearray
        Mutation vector (initially blank) into which to write bytes
    region_seq: bytes
        Sequence of the region for which to compute the mutation vector
        (only the region, not any other part of the reference sequence)
    region_length: int (≥ region_end5)
        3'-most coordinate of the region with respect to the entire
        reference sequence (1-indexed, includes coordinate)
    region_end5: int (≥ 1)
        5'-most coordinate of the region with respect to the entire
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
    if region_length != len(muts):
        raise VectorValueError(f"Region is {region_length} nt, but mutation "
                               f"vector is {len(muts)} nt.")
    if region_length != len(region_seq):
        raise VectorValueError(f"Region is {region_length} nt, but sequence "
                               f"is {len(region_seq)} nt.")
    # Indexes of the 5' and 3' ends of the current CIGAR operation with
    # respect to the read; 0-indexed, using Python half-open intervals
    read_idx5 = 0
    read_idx3 = 0
    # Indexes of the 5' and 3' ends of the current CIGAR operation with
    # respect to the region; 0-indexed, using Python half-open intervals
    region_idx5 = read.pos - region_end5
    region_idx3 = region_idx5
    # Number of bases truncated from the end of the operation
    truncated = 0
    # Record all deletions and insertions.
    dels: list[Deletion] = list()
    inns: list[Insertion] = list()
    # Read the CIGAR string one operation at a time.
    for cigar_op, op_length in parse_cigar(read.cigar):
        # Update the coordinates, with respect to the region and read,
        # that correspond to the 3' end of the current CIGAR operation.
        if op_consumes_ref(cigar_op):
            region_idx3 += op_length
        if op_consumes_read(cigar_op):
            read_idx3 += op_length
        # Check if the part of the read that corresponds to the current
        # CIGAR operation at all overlaps the region of interest.
        # Note: This loop does not terminate when the CIGAR operation
        # exits the region of interest, even though no more matches or
        # mutations will be appended to the end of mut_vectors afterwards.
        # This behavior forces the entire CIGAR string to be read, which
        # is necessary for parse_cigar to validate the CIGAR string.
        if region_idx3 > 0 and region_idx5 < region_length:
            if region_idx5 < 0:
                # If the 3' end of the CIGAR operation overlaps the
                # region but the 5' end does not, then truncate the
                # 5' end of the current CIGAR operation so that the
                # operation starts at the 5' end of the region.
                # Decrease the length of the CIGAR operation because it
                # is being truncated. Note: region_idx5 < 0.
                op_length += region_idx5
                if op_consumes_read(cigar_op):
                    # If the CIGAR operation consumes the read, advance
                    # the index of the 5' end of the CIGAR operation
                    # with respect to the read. Note: region_idx5 < 0.
                    read_idx5 -= region_idx5
                # Advance the index of the 5' end of the CIGAR operation
                # with respect to the region to 0; that is, the CIGAR
                # operation now starts at the 5' end of the region.
                region_idx5 = 0
            if region_idx3 > region_length:
                # If the 5' end of the CIGAR operation overlaps the
                # region but the 3' end does not, then truncate the
                # 3' end of the current CIGAR operation so that the
                # operation ends at the 3' end of the region.
                # First, find the number of positions to truncate from
                # the 3' end; truncated is guaranteed to be > 0.
                truncated = region_idx3 - region_length
                # Decrease the length of the CIGAR operation because it
                # is being truncated.
                op_length -= truncated
                if op_consumes_read(cigar_op):
                    # If the CIGAR operation consumes the read, reduce
                    # the index of the 3' end of the CIGAR operation
                    # with respect to the read.
                    read_idx3 -= truncated
                # Reduce the index of the 3' end of the CIGAR operation
                # with respect to the region to equal the length of the
                # region; that is, the CIGAR operation now ends at the
                # 3' end of the region.
                region_idx3 = region_length
            # Act based on the CIGAR operation and its length.
            if cigar_op == CIG_MATCH:
                # The read and reference sequences match over the entire
                # CIGAR operation.
                for base, qual in zip(read.seq[read_idx5: read_idx3],
                                      read.qual[read_idx5: read_idx3]):
                    muts[region_idx5] = encode_match(base, qual, min_qual)
                    region_idx5 += 1
            elif cigar_op == CIG_ALIGN or cigar_op == CIG_SUBST:
                # The read contains only matches or substitutions (no
                # indels) relative to the reference over the entire
                # CIGAR operation.
                for ref_base, read_base, read_qual in zip(
                        region_seq[region_idx5: region_idx3],
                        read.seq[read_idx5: read_idx3],
                        read.qual[read_idx5: read_idx3],
                        strict=True):
                    muts[region_idx5] = encode_compare(ref_base, read_base,
                                                       read_qual, min_qual)
                    region_idx5 += 1
            elif cigar_op == CIG_DELET:
                # The portion of the reference sequence corresponding
                # to the CIGAR operation is deleted from the read.
                # Create one Deletion object for each base in the
                # reference sequence that is missing from the read.
                while region_idx5 < region_idx3:
                    dels.append(Deletion(region_idx5, read_idx5))
                    muts[region_idx5] = DELET
                    region_idx5 += 1
            elif cigar_op == CIG_INSRT:
                # The read contains an insertion of one or more bases
                # that are not present in the reference sequence.
                # Create one Insertion object for each base in the read
                # sequence that is not present in the reference. Every
                # mutation needs to be assigned a coordinate in the
                # region in order to appear at that coordinate in the
                # mutation vector. But each inserted base, being absent
                # from the reference, does not correspond to a single
                # coordinate in the region; instead, each inserted base
                # lies between two coordinates in the region. Either of
                # these coordinates could be chosen; this code assigns
                # the 3' coordinate to the insertion. For example, if
                # two bases are inserted between coordinates 45 and 46
                # of the region, then both will be given coordinate 46.
                # The reason for this convention is that the math is
                # simpler than it would be if using the 5' coordinate.
                # Because region_idx5 is, by definition, immediately 3'
                # of the previous CIGAR operation; and the bases that
                # are inserted lie between the previous and subsequent
                # CIGAR operations; region_idx5 is the coordinate
                # immediately 3' of the inserted bases. In this special
                # case, region_idx5 also equals region_idx3 (because
                # the insertion does not consume the reference, so
                # region_idx3 += op_length was not run at the beginning
                # of this loop iteration), as well as the length of mut_vectors
                # (because Python is 0-indexed, so the length of a range
                # of indexes such as [0, 1, ... , 45] equals the value
                # of the next index in the range, 46). Thus, there are
                # three variables that already equal the 3' coordinate
                # and none that equal the 5' coordinate.
                while read_idx5 < read_idx3:
                    inns.append(Insertion(read_idx5, region_idx5))
                    read_idx5 += 1
                # Insertions do not consume the reference, so do not add
                # any information to mut_vectors yet; it will be added later.
            elif cigar_op == CIG_SCLIP:
                # Bases were soft-clipped from the 5' or 3' end of the
                # read during alignment. Like insertions, they consume
                # the read but not the reference. Unlike insertions,
                # they are not mutations, so they do not require any
                # additional processing.
                pass
            else:
                raise VectorValueError(
                    f"Invalid CIGAR operation: '{cigar_op.decode()}'")
        if truncated:
            # If the current operation was truncated because it extended
            # past the 3' end of the region, then the 3' ends of the
            # region and read need to be reset to their values before
            # truncation. Otherwise, all subsequent CIGAR operations
            # will start and end 5' of where they should.
            region_idx3 += truncated
            if op_consumes_read(cigar_op):
                read_idx3 += truncated
            truncated = 0
        # Advance the 5' positions in the region and read to the current
        # 3' positions so that the next CIGAR operation lies immediately
        # 3' of the current operation. The 3' positions will be advanced
        # at the beginning of the next iteration (if any) of the loop.
        region_idx5 = region_idx3
        read_idx5 = read_idx3
    if len(muts) != region_length:
        raise VectorValueError(f"Mutation vector is {len(muts)} nt, "
                               f"but region is {region_length} nt.")
    # Verify that the sum of all CIGAR operations that consumed the read
    # equals the length of the read. The former equals read_idx5 because
    # for each operation that consumed the read, the length of the
    if read_idx5 != len(read.seq):
        raise VectorValueError(
            f"CIGAR string '{read.cigar.decode()}' consumed {read_idx5} "
            f"bases from read, but read is {len(read.seq)} bases long.")
    # Add insertions to mut_vectors.
    for ins in inns:
        ins.stamp(muts)
    # Label all positions that are ambiguous due to indels.
    if ambid and (dels or inns):
        get_ambids(muts, region_seq, read.seq, read.qual, min_qual, dels, inns)


def vectorize_line(line: bytes, muts: bytearray, seq: bytes,
                   length: int, end5: int, ref: str, qmin: int, ambid: bool):
    read = SamRead(line)
    if read.rname != ref:
        raise VectorValueError(f"Read '{read.qname.decode()}' had reference "
                               f"'{read.rname}' (≠ '{ref}')")
    vectorize_read(read, muts, seq, length, end5, qmin, ambid)


def vectorize_pair(line1: bytes, line2: bytes, muts: bytearray, seq: bytes,
                   length: int, end5: int, ref: str, qmin: int, ambid: bool):
    # Parse lines 1 and 2 into SAM reads.
    read1 = SamRead(line1)
    read2 = SamRead(line2)
    # Ensure that reads 1 and 2 are compatible mates.
    if not read1.flag.paired:
        raise VectorValueError(f"Read 1 ({read1.qname.decode()}) was not "
                               f"paired, but read 2 ('{read2.qname.decode()}') "
                               f"was given")
    if not read2.flag.paired:
        raise VectorValueError(f"Read 2 ({read2.qname.decode()}) was not "
                               f"paired, but read 1 ({read1.qname.decode()}) "
                               f"was given")
    if read1.qname != read2.qname:
        raise VectorValueError(f"Reads 1 ({read1.qname.decode()}) and 2 "
                               f"({read2.qname.decode()}) had different names")
    if read1.rname != read2.rname:
        raise VectorValueError(f"Read '{read1.qname.decode()}' had "
                               "different references for mates 1 "
                               f"('{read1.rname}') and 2 "
                               f"('{read2.rname}')")
    if read1.rname != ref:
        raise VectorValueError(f"Read '{read1.qname.decode()}' had reference "
                               f"'{read1.rname}' (≠ '{ref}')")
    if not (read1.flag.first and read2.flag.second):
        raise VectorValueError(f"Read '{read1.qname.decode()}' had mate 1 "
                               f"labeled {2 - read1.flag.first} and mate 2 "
                               f"labeled {1 + read2.flag.second}")
    if read1.flag.rev == read2.flag.rev:
        raise VectorValueError(f"Read '{read1.qname.decode()}' had "
                               "mates 1 and 2 facing the same way")
    # Vectorize read 1.
    muts1 = muts.copy()
    vectorize_read(read1, muts1, seq, length, end5, qmin, ambid)
    # Vectorize read 2.
    muts2 = muts.copy()
    vectorize_read(read2, muts2, seq, length, end5, qmin, ambid)
    # Compute the consensus of reads 1 and 2.
    for i, (m1, m2) in enumerate(zip(muts1, muts2, strict=True)):
        muts[i] = inter if (inter := m1 & m2) else m1 | m2
