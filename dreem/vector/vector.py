from __future__ import annotations
import re
from typing import List, Optional


from dreem.util.util import BASES, SUB_A, SUB_C, SUB_G, SUB_T, SUB_N, MATCH, DELET, ANY_N, INS_3, INS_5, BLANK, DNA, PHRED_ENCODING



CIG_ALN = b"M"  # alignment match
CIG_MAT = b"="  # sequence match
CIG_SUB = b"X"  # substitution
CIG_DEL = b"D"  # deletion
CIG_INS = b"I"  # insertion
CIG_SCL = b"S"  # soft clipping
CIG_PATTERN = re.compile(b"".join(
    [rb"(\d+)([", CIG_ALN, CIG_MAT, CIG_SUB, CIG_DEL, CIG_INS, CIG_SCL, b"])"]
))
SAM_HEADER = b"@"
NOLIM = -1
A_INT, C_INT, G_INT, T_INT = BASES
ANY_N_INT = ANY_N[0]
MATCH_INT = MATCH[0]
SUB_A_INT = SUB_A[0]
SUB_C_INT = SUB_C[0]
SUB_G_INT = SUB_G[0]
SUB_T_INT = SUB_T[0]
MIN_QUAL = 25
MIN_QUAL_CODE = MIN_QUAL + PHRED_ENCODING


def encode_base(base: int):
    if base == T_INT:
        return SUB_T_INT
    if base == G_INT:
        return SUB_G_INT
    if base == C_INT:
        return SUB_C_INT
    if base == A_INT:
        return SUB_A_INT
    raise ValueError(f"Invalid base: {base.to_bytes().decode()}")


def encode_compare(ref_base: int, read_base: int, read_qual: int):
    return ((MATCH_INT if ref_base == read_base
             else encode_base(read_base)) if read_qual >= MIN_QUAL_CODE
             else ANY_N_INT ^ encode_base(ref_base))


def encode_match(ref_base: int, read_qual: int):
    # A more efficient version of encode_compare given the prior knowledge from
    # the CIGAR string that the read and reference match at this position.
    # NOTE: there is no analagous version when there is a known substitution
    # because substitutions are infrequent, so optimizing their processing
    # would speed the program insignificantly while making the source code
    # more complex and harder to maintain.
    return (MATCH_INT if read_qual >= MIN_QUAL_CODE
            else ANY_N_INT ^ encode_base(ref_base))


class Indel(object):
    """
    Base class for an Insertion or Deletion.

    It is used to find alternative positions for insertions and deletions,
    by keeping track of an indel's current coordinates (as they are moved)
    and determining whether a specific move is valid.

    Arguments
    rel_ins_idx (int):  The 0-indexed position of the indel with respect to
                        the sequence (ref or read) with a relative insertion
                        (that is, the read if the mutation is denoted an
                        insertion, and the ref if a deletion).
                        This index points to one specific base, either
                        (for insertions) a base inserted into the read,
                        in the coordinates of the read or (for deletions)
                        a base present in the ref and absent from the read,
                        in the coordinates of the ref.
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

    __slots__ = ["_ins_idx", "_ins_init", "_del_idx", "_del_init"]

    MIN_INDEL_DIST = 2

    def __init__(self, rel_ins_idx: int, rel_del_idx: int) -> None:
        self._ins_idx = rel_ins_idx
        self._ins_init = rel_ins_idx
        self._del_idx = rel_del_idx
        self._del_init = rel_del_idx

    @property
    def ins_idx(self):
        return self._ins_idx
    
    @property
    def del_idx5(self):
        return self._del_idx - 1

    @property
    def del_idx3(self):
        return self._del_idx
    
    def ordering(self) -> int:
        raise NotImplementedError

    def reset_idx(self):
        self._ins_idx = self._ins_init
        self._del_idx = self._del_init
    
    @staticmethod
    def _get_indel_by_idx(indels: List[Indel], idx: int):
        for indel in indels:
            if indel.ins_idx == idx:
                return indel

    def _peek_out_of_indel(self, indels: List[Indel], from5_to3: bool):
        idx = self.ins_idx + (inc := 1 if from5_to3 else -1)
        while self._get_indel_by_idx(indels, idx):
            idx += inc
        return idx
    
    def _collision(self, other: Indel, swap_idx: int):
        return self.MIN_INDEL_DIST > (min(abs(swap_idx - other.del_idx5),
                                          abs(swap_idx - other.del_idx3)))
    
    def _collisions(self, indels: List[Indel], swap_idx: int):
        return any(self._collision(indel, swap_idx) for indel in indels)
    
    def _try_swap(self, *args, **kwargs) -> bool:
        raise NotImplementedError
    
    def sweep(self, muts: bytearray, ref: bytes, read: SamRead,
              dels: List[Deletion], inns: List[Insertion], from5_to3: bool):
        while self._try_swap(muts, ref, read, dels, inns, from5_to3):
            pass


class Deletion(Indel):
    def ordering(self):
        return self._ins_idx
    
    @staticmethod
    def _encode_swap(ref_base: int, swap_base: int,
                     read_base: int, read_qual: int):
        swap_code = encode_compare(swap_base, read_base, read_qual)
        if ref_base == swap_base:
            # Swap occurs between two identical positions in reference:
            # always valid
            return swap_code
        if swap_code == MATCH_INT:
            # The destination in the reference matches the read and differs
            # from the starting position in the reference:
            # never valid 
            return 0
        if swap_code & MATCH_INT:
            # The destination in the read is low quality, so it could match the
            # reference or be a substitution to anything but the reference.
            return ANY_N_INT ^ encode_base(ref_base)
        # The destination in the read must be a substitution.
        if ref_base == read_base:
            # The substitution cannot be to the same base as the reference.
            return 0
        # Otherwise, return the substitution.
        return encode_base(read_base)

    def _swap(self, muts: bytearray, swap_idx: int, swap_code: int):
        """
        Arguments
        muts (bytearray): mutation vector
        swap_idx (int): the index in the region to which the deletion moves
                        during this swap
        swap_code (int): the relationship (match, sub, etc.) between the
                         base located at swap_idx and the base in the read
        """
        # The base at swap_idx moves to self.ref_idx, so after the swap, the
        # relationship between self.ref_idx and the read base will be swap_code.
        muts[self.ins_idx] = muts[self.ins_idx] | swap_code
        # The base at self.ref_idx is marked as a deletion (by definition), so
        # mark the position it moves to (swap_idx) as a deletion too.
        muts[swap_idx] = muts[swap_idx] | DELET
        # Move the indel's position (self.ins_idx) to swap_idx.
        # Move self.del_idx by the same amount.
        self._del_idx += swap_idx - self.ins_idx
        self.ins_idx = swap_idx
    
    def _try_swap(self, muts: bytearray, ref: bytes, read: SamRead,
                  dels: List[Deletion], inns: List[Insertion],
                  from5_to3: bool) -> bool:
        swap_idx = self._peek_out_of_indel(dels, from5_to3)
        if 0 <= swap_idx < len(ref) and not self._collisions(inns, swap_idx):
            read_idx = self.del_idx3 if from5_to3 else self.del_idx5
            swap_code = self._encode_swap(ref[self.ins_idx],
                                          ref[swap_idx],
                                          read.seq[read_idx],
                                          read.qual[read_idx])
            if swap_code:
                self._swap(muts, swap_idx, swap_code)
                return True
        return False


class Insertion(Indel):
    def ordering(self):
        return self._del_idx


def sweep_indels(muts: bytearray, ref: bytes, read: SamRead,
                 dels: List[Deletion], inns: List[Insertion],
                 from5_to3: bool):
    for indel in (indels := dels + inns):
        indel.reset_idx()
    indels.sort(key=Indel.ordering, reverse=from5_to3)
    while indels:
        indel = indels.pop()
        indel.sweep(muts, ref, read, dels, inns, from5_to3)
        idx = len(indels)
        while idx > 0 and (from5_to3 == indels[idx - 1].ordering()
                           < indel.ordering()):
            idx -= 1
        if idx < len(indels):
            indels.insert(idx, indel)


def all_indels(muts: bytearray, ref: bytes, read: SamRead,
               dels: List[Deletion], inns: List[Insertion]):
    if dels or inns:
        for from5_to3 in (True, False):
            sweep_indels(muts, ref, read, dels, inns, from5_to3)


def parse_cigar(cigar_string: bytes):
    if not cigar_string:
        raise ValueError("CIGAR string was empty.")
    length_matched = 0
    for match in CIG_PATTERN.finditer(cigar_string):
        length_bytes, operation = match.groups()
        length_matched += len(length_bytes) + len(operation)
        if (length_int := int(length_bytes)) < 1:
            raise ValueError("length of CIGAR operation must be >= 1")
        yield operation, length_int
    if length_matched != len(cigar_string):
        raise ValueError(f"Invalid CIGAR string: '{cigar_string.decode()}'")


def op_consumes_ref(op: bytes):
    return op == CIG_ALN or op == CIG_MAT or op == CIG_SUB or op == CIG_DEL


class SamFlag(object):
    __slots__ = ["paired", "proper", "unmap", "munmap", "rev", "mrev",
                 "first", "second", "secondary", "qcfail", "dup", "supp"]

    max_flag: int = 2**len(__slots__) - 1
    pattern = "".join(["{:0<", str(len(__slots__)), "}"])

    def __init__(self, flag: int):
        if not 0 <= flag <= self.max_flag:
            raise ValueError(f"Invalid flag: '{flag}'")
        (self.paired, self.proper, self.unmap, self.munmap, self.rev,
         self.mrev, self.first, self.second, self.secondary, self.qcfail,
         self.dup, self.supp) = (x == "1" for x in
                                 self.pattern.format(bin(flag)[:1:-1]))


class SamRead(object):
    __slots__ = ["qname", "flag", "rname", "pos", "mapq", "cigar",
                 "tlen", "seq", "qual"]
    
    min_fields = 11

    def __init__(self, line: bytes):
        fields = line.rstrip().split(b"\t")
        if len(fields) < self.min_fields:
            raise ValueError(f"Invalid SAM line:\n{line}")
        self.qname = fields[0]
        self.flag = SamFlag(int(fields[1]))
        self.rname = fields[2]
        self.pos = int(fields[3])
        self.mapq = int(fields[4])
        self.cigar = fields[5]
        #self.rnext = fields[6]
        #self.pnext = int(fields[7])
        self.tlen = int(fields[8])
        self.seq = DNA(fields[9])
        self.qual = fields[10]
        if len(self) != len(self.qual):
            raise ValueError(f"Lengths of seq ({len(self)}) and qual "
                             f"string {len(self.qual)} did not match.")
        
    def __len__(self):
        return len(self.seq)


def vectorize_read(region_seq: bytes, region_first: int, region_last: int,
                   read: SamRead):
    """
    :param region_seq: str, reference sequence (must contain T, not U)
    :param region_first: int, the first coordinate (w.r.t. the reference
        sequence) of the region of interest (1-indexed)
    :param region_last: int, the last coordinate (w.r.t. the reference
        sequence) of the region of interest, inclusive (1-indexed)
    :param read: SamRead, the read for which to compute mutations
    :return:
    """
    region_length = region_last - region_first + 1
    assert region_length == len(region_seq)
    # current position in the read
    # 0-indexed from beginning of read
    read_idx = 0
    # position at which the current CIGAR operation starts
    # 0-indexed from beginning of region
    op_start_idx = read.pos - region_first
    # position at which the current CIGAR operation ends (0-indexed)
    # does not include the last position of the operation (like Python slicing)
    # 0-indexed from beginning of region
    op_end_idx = op_start_idx
    # Initialize the mutation vector. Pad the beginning with missing bytes
    # if the read starts after the first position in the region.
    muts = bytearray(BLANK * min(op_start_idx, region_length))
    # Record all deletions and insertions.
    dels: List[Deletion] = list()
    inns: List[Insertion] = list()
    # Read the CIGAR string one operation at a time.
    for cigar_op, op_length in parse_cigar(read.cigar):
        if op_consumes_ref(cigar_op):
            # Advance the end of the operation if it consumes the reference.
            op_end_idx += op_length
        if op_end_idx > 0 and op_start_idx < region_length:
            # Run this block once the operation has entered the region.
            if op_start_idx < 0:
                # If the current operation starts before the region,
                # then advance it to the beginning of the region.
                op_length += op_start_idx
                op_start_idx = 0
            if op_end_idx > region_length:
                # If the current operation ends after the region,
                # then truncate it to the end of the region.
                op_length -= op_end_idx - region_length
                op_end_idx = region_length
            # Perform an action based on the CIGAR operation and its length.
            if cigar_op == CIG_MAT:
                # Condition: read matches reference
                next_read_idx = read_idx + op_length
                muts.extend(map(encode_match,
                                read.seq[read_idx: next_read_idx],
                                read.qual[read_idx: next_read_idx]))
                read_idx = next_read_idx
            elif cigar_op == CIG_ALN or cigar_op == CIG_SUB:
                # Condition: read has a match or substitution relative to ref
                next_read_idx = read_idx + op_length
                muts.extend(map(encode_compare,
                                region_seq[len(muts): len(muts) + op_length],
                                read.seq[read_idx: next_read_idx],
                                read.qual[read_idx: next_read_idx]))
                read_idx = next_read_idx
            elif cigar_op == CIG_DEL:
                # Condition: read contains a deletion w.r.t. the reference.
                dels.extend(Deletion(ref_idx, read_idx) for ref_idx
                            in range(len(muts), len(muts) + op_length))
                muts.extend(DELET * op_length)
            elif cigar_op == CIG_INS:
                # Condition: read contains an insertion w.r.t. the reference.
                # Position added to insertions is of the base 3' of the insert.
                next_read_idx = read_idx + op_length
                inns.extend(Insertion(read_idx, len(muts)) for read_idx
                              in range(read_idx, next_read_idx))
                # Insertions do not consume the reference, so do not add any
                # information to muts yet. That information is added later.
                read_idx = next_read_idx
            elif cigar_op == CIG_SCL:
                # Condition: read contains a soft clipping
                read_idx += op_length
            else:
                raise ValueError(
                    f"Invalid CIGAR operation: '{cigar_op.decode()}'")
        # Advance the start position to the end of the current operation.
        op_start_idx = op_end_idx
    # Pad the end of the mutation vector with any non-covered positions.
    muts.extend(BLANK * (region_length - len(muts)))
    assert len(muts) == region_length
    # Ensure the CIGAR string matched the length of the read.
    if read_idx != len(read):
        raise ValueError(
            f"CIGAR string '{read.cigar.decode()}' consumed {read_idx} "
            f"bases from read, but read is {len(read)} bases long.")
    # Flag all positions that are ambiguous due to indels.
    all_indels(muts, region_seq, read, dels, inns)
    return muts


def get_consensus_mut(byte1: int, byte2: int):
    return intersect if (intersect := byte1 & byte2) else byte1 | byte2


def vectorize_pair(region_seq: bytes, region_first: int, region_last: int,
                   read1: SamRead, read2: SamRead):
    muts1 = vectorize_read(region_seq, region_first, region_last, read1)
    muts2 = vectorize_read(region_seq, region_first, region_last, read2)
    return bytearray(map(get_consensus_mut, muts1, muts2))


class SamRecord(object):
    __slots__ = ["read1", "read2"]

    def __init__(self, read1: SamRead, read2: Optional[SamRead] = None):
        self.read1 = read1
        if read2 is not None:
            if self.paired:
                if not read2.flag.paired:
                    raise ValueError("read1 is paired but read2 is not")
                if read1.qname != read2.qname:
                    raise ValueError(
                        "Paired reads had inconsistent query names: "
                        f"'{read1.qname}' and '{read2.qname}'")
                if read1.rname != read2.rname:
                    raise ValueError(
                        "Paired reads had inconsistent reference names: "
                        f"'{read1.rname}' and '{read2.rname}'")
                if abs(read1.tlen) != abs(read2.tlen):
                    raise ValueError(
                        "Paired reads had inconsistent template lengths: "
                        f"{read1.tlen} and {read2.tlen}")
                if read1.flag.second or not read1.flag.first:
                    raise ValueError("read1 is not flagged as first read")
                if read2.flag.first or not read2.flag.second:
                    raise ValueError("read2 is not flagged as second read")
                if read1.flag.rev == read2.flag.rev:
                    raise ValueError(
                        "read1 and read2 are in the same orientation")
            else:
                raise ValueError("read1 is unpaired, but read2 was given")
        self.read2 = read2

    @property
    def ref_name(self):
        return self.read1.rname

    @property
    def paired(self):
        return self.read1.flag.paired

    def vectorize(self, region_seq, region_start, region_end):
        if self.read2 is None:
            return vectorize_read(region_seq, region_start, region_end,
                                  self.read1)
        else:
            return vectorize_pair(region_seq, region_start, region_end,
                                  self.read1, self.read2)


'''
def query_muts(muts: np.ndarray, bits: int):
    """
    Count the number of times a query mutation occurs in each column
    or one column of a set of mutation vectors.
    The counting operation comprises three steps:
    1. bitwise AND to confirm at least one "1" bit is shared, e.g.
       bits: 11110000 & muts: 00100000 -> 00100000 (True)
       bits: 11110000 & muts: 00100010 -> 00100000 (True)
       bits: 11110000 & muts: 00000000 -> 00000000 (False)
    2. bitwise OR to confirm no "1" bit in muts is not in bits, e.g.
       bits: 11110000 | muts: 00100000 -> 11110000 =? 11110000 (True)
       bits: 11110000 | muts: 00100010 -> 11110010 =? 11110000 (False)
       bits: 11110000 | muts: 00000000 -> 11110000 =? 11110000 (True)
    3. logical AND to confirm that both tests pass, e.g.
       bits: 11110000, muts: 00100000 -> True  AND True  (True)
       bits: 11110000, muts: 00100010 -> True  AND False (False)
       bits: 11110000, muts: 00000000 -> False AND True  (False)

    Arguments
    muts: NDArray of a set of mutation vectors (2-dimensional)
          or one column in a set of mutation vectors (1-dimensional).
          Data type must be uint8.
    bits: One-byte int in the range [0, 256) representing the mutation
          to be queried. The bits in the int encode the mutation as
          defined above, e.g.
          - 00000010 (int 2) is a deletion
          - 11010001 (int 209) is either substitution to A, G, or T
                               or a match to C
    
    Returns
    count: If muts is 1-dimensional, int of the number of times the
           query mutation occurs in muts.
           If muts is 2-dimensional, NDArray with one int for each
           column in muts.
    """
    assert muts.dtype is np.uint8
    assert isinstance(bits, int) and 0 <= bits < 256
    count = np.logical_and(muts & bits, (muts | bits) == bits).sum(axis=0)
    return count

def count_muts(muts: np.ndarray):
    out = dict()
    out["match_bases"] = query_muts(muts, MATCH[0])
    out["mod_bases_A"] = query_muts(muts, SUB_A[0])
    out["mod_bases_C"] = query_muts(muts, SUB_C[0])
    out["mod_bases_G"] = query_muts(muts, SUB_G[0])
    out["mod_bases_T"] = query_muts(muts, SUB_T[0])
    out["mod_bases_N"] = query_muts(muts, SUB_N[0])
    out["del_bases"]   = query_muts(muts, DELET[0])
    out["ins_bases"]   = query_muts(muts, INS_5[0] | INS_3[0])
    # Can have any mutation, but not a match
    out["mut_bases"] = query_muts(muts, SUB_N[0] | DELET[0] | INS_5[0] | INS_3[0])
    out["cov_bases"] = muts.astype(bool).sum(axis=0)  # i.e. not BLANK
    # Unambiguously matching or mutated (informative)
    out["info_bases"] = out["match_bases"] + out["mut_bases"]
    # Mutation rate (fraction mutated among all unambiguously matching/mutated)
    try:
        out["mut_rates"] = out["mut_bases"] / out["info_bases"]
    except ZeroDivisionError:
        out["mut_rates"] = np.nan
    return out
'''
