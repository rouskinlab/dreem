from __future__ import annotations
import re
from typing import List, Optional

from dreem.util.util import BASES, SUB_A, SUB_C, SUB_G, SUB_T, SUB_N, MATCH, DELET, ANY_N, INS_3, INS_5, BLANK, DNA, DEFAULT_PHRED_ENCODING


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
MATCH_INT = MATCH[0]
DELET_INT = DELET[0]
INS_5_INT = INS_5[0]
INS_3_INT = INS_3[0]
SUB_A_INT = SUB_A[0]
SUB_C_INT = SUB_C[0]
SUB_G_INT = SUB_G[0]
SUB_T_INT = SUB_T[0]
SUB_N_INT = SUB_N[0]
ANY_N_INT = ANY_N[0]
MIN_QUAL_PHRED = 20
MIN_QUAL_PCODE = MIN_QUAL_PHRED + DEFAULT_PHRED_ENCODING


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
    return ((MATCH_INT if ref_base == read_base else encode_base(read_base))
            if read_qual >= MIN_QUAL_PCODE
            else ANY_N_INT ^ encode_base(ref_base))


def encode_match(read_base: int, read_qual: int):
    # A more efficient version of encode_compare given the prior knowledge from
    # the CIGAR string that the read and reference match at this position.
    # NOTE: there is no analagous version when there is a known substitution
    # because substitutions are infrequent, so optimizing their processing
    # would speed the program insignificantly while making the source code
    # more complex and harder to maintain.
    return (MATCH_INT
            if read_qual >= MIN_QUAL_PCODE
            else ANY_N_INT ^ encode_base(read_base))


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

    __slots__ = ["_ins_idx", "_ins_init", "_del_idx", "_del_init", "_tunneled"]

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
        raise NotImplementedError

    def reset(self):
        self._ins_idx = self._ins_init
        self._del_idx = self._del_init
        self._tunneled = False
    
    @staticmethod
    def _get_indel_by_idx(indels: List[Indel], idx: int):
        for indel in indels:
            if indel.ins_idx == idx:
                return indel

    def _peek_out_of_indel(self, indels: List[Indel], from3to5: bool):
        inc = -1 if from3to5 else 1
        idx = self.ins_idx + inc
        tunneled_indels: List[Indel] = list()
        while indel := (self._get_indel_by_idx(indels, idx)):
            idx += inc
            tunneled_indels.append(indel)
        self._tunneled = bool(tunneled_indels)
        return idx, tunneled_indels
    
    def _collision(self, other: Indel, swap_idx: int):
        return self.MIN_INDEL_DIST > (min(abs(swap_idx - other.del_idx5),
                                          abs(swap_idx - other.del_idx3)))
    
    def _collisions(self, indels: List[Indel], swap_idx: int):
        return any(self._collision(indel, swap_idx) for indel in indels)
    
    def step_del_idx(self, swap_idx: int):
        # Move the indel's position (self._ins_idx) to swap_idx.
        # Move self._del_idx one step in the same direction.
        assert swap_idx != self.ins_idx
        self._del_idx += 1 if swap_idx > self.ins_idx else -1
    
    def _step(self, swap_idx: int):
        self.step_del_idx(swap_idx)
        self._ins_idx = swap_idx
    
    @staticmethod
    def _consistent_rels(curr_rel: int, swap_rel: int):
        if curr_rel & swap_rel or (curr_rel & SUB_N_INT
                                   and swap_rel & SUB_N_INT):
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
        raise NotImplementedError
    
    def _try_swap(self, *args, **kwargs) -> bool:
        raise NotImplementedError
    
    def sweep(self, muts: bytearray, ref: bytes, read: bytes, qual: bytes,
              dels: List[Deletion], inns: List[Insertion], from3to5: bool,
              tunnel: bool):
        # Move the indel as far as possible in either the 5' or 3' direction.
        while self._try_swap(muts, ref, read, qual, dels, inns, from3to5,
                             tunnel):
            # All actions happen in _try_swap, so loop body is empty.
            pass


class Deletion(Indel):
    @property
    def rank(self):
        return self._ins_idx
    
    @classmethod
    def _encode_swap(cls, ref_base: int, swap_base: int,
                     read_base: int, read_qual: int):
        curr_rel = encode_compare(ref_base, read_base, read_qual)
        swap_rel = encode_compare(swap_base, read_base, read_qual)
        return cls._consistent_rels(curr_rel, swap_rel)

    def _swap(self, muts: bytearray, swap_idx: int, relation: int):
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
        muts[self.ins_idx] = muts[self.ins_idx] | relation
        # The base at self.ref_idx is marked as a deletion (by definition), so
        # mark the position it moves to (swap_idx) as a deletion too.
        muts[swap_idx] = muts[swap_idx] | DELET_INT
        self._step(swap_idx)
    
    def _try_swap(self, muts: bytearray, ref: bytes, read: bytes, qual: bytes,
                  dels: List[Deletion], inns: List[Insertion],
                  from3to5: bool, tunnel: bool) -> bool:
        swap_idx, tunneled_indels = self._peek_out_of_indel(dels, from3to5)
        read_idx = self.del_idx5 if from3to5 else self.del_idx3
        if (0 <= swap_idx < len(ref) and 0 <= read_idx < len(read)
                and (tunnel or not self.tunneled)
                and not self._collisions(inns, swap_idx)):
            relation = self._encode_swap(ref[self.ins_idx], ref[swap_idx],
                                         read[read_idx], qual[read_idx])
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
            muts[self.del_idx5] = muts[self.del_idx5] | INS_5_INT
        if 0 <= self.del_idx3 < len(muts):
            muts[self.del_idx3] = muts[self.del_idx3] | INS_3_INT

    @classmethod
    def _encode_swap(cls, ref_base: int, read_base: int,
                     read_qual: int, swap_base: int, swap_qual: int):
        curr_rel = encode_compare(ref_base, read_base, read_qual)
        swap_rel = encode_compare(ref_base, swap_base, swap_qual)
        return cls._consistent_rels(curr_rel, swap_rel)

    def _swap(self, muts: bytearray, ref_idx: int,
              swap_idx: int, relation: int):
        """
        Arguments
        muts (bytearray): mutation vector
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
                  dels: List[Deletion], inns: List[Insertion],
                  from3to5: bool, tunnel: bool) -> bool:
        swap_idx, tunneled_indels = self._peek_out_of_indel(inns, from3to5)
        ref_idx = self.del_idx5 if from3to5 else self.del_idx3
        if (0 <= swap_idx < len(read) and 0 <= ref_idx < len(ref)
                and (tunnel or not self.tunneled)
                and not self._collisions(dels, swap_idx)):
            relation = self._encode_swap(ref[ref_idx],
                                         read[self.ins_idx],
                                         qual[self.ins_idx],
                                         read[swap_idx],
                                         qual[swap_idx])
            if relation:
                self._swap(muts, ref_idx, swap_idx, relation)
                for indel in tunneled_indels:
                    indel.step_del_idx(swap_idx)
                return True
        return False


def sweep_indels(muts: bytearray, ref: bytes, read: bytes, qual: bytes,
                 dels: List[Deletion], inns: List[Insertion], from3to5: bool,
                 tunnel: bool):
    indels = dels + inns
    for indel in indels:
        indel.reset()
    sort_rev = from3to5 != tunnel
    indels.sort(key=lambda indel: indel.rank, reverse=sort_rev)
    while indels:
        indel = indels.pop()
        indel.sweep(muts, ref, read, qual, dels, inns, from3to5, tunnel)
        i = len(indels)
        if sort_rev:
            while i > 0 and indel.rank > indels[i - 1].rank:
                i -= 1
        else:
            while i > 0 and indel.rank < indels[i - 1].rank:
                i -= 1
        if i < len(indels):
            indels.insert(i, indel)


def allindel(muts: bytearray, ref: bytes, read: bytes, qual: bytes,
             dels: List[Deletion], inns: List[Insertion]):
    for from3to5 in (False, True):
        sweep_indels(muts, ref, read, qual, dels, inns, from3to5, True)
        if any(d.tunneled for d in dels) or any(i.tunneled for i in inns):
            sweep_indels(muts, ref, read, qual, dels, inns, from3to5, False)


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
    return op != CIG_INS and op != CIG_SCL


def op_consumes_read(op: bytes):
    return op != CIG_DEL


class SamFlag(object):
    __slots__ = ["paired", "proper", "unmap", "munmap", "rev", "mrev",
                 "first", "second", "secondary", "qcfail", "dup", "supp"]

    MAX_FLAG: int = 2**len(__slots__) - 1
    PATTERN = "".join(["{:0<", str(len(__slots__)), "}"])

    def __init__(self, flag: int):
        if not 0 <= flag <= self.MAX_FLAG:
            raise ValueError(f"Invalid flag: '{flag}'")
        (self.paired, self.proper, self.unmap, self.munmap, self.rev,
         self.mrev, self.first, self.second, self.secondary, self.qcfail,
         self.dup, self.supp) = (x == "1" for x in
                                 self.PATTERN.format(bin(flag)[:1:-1]))


class SamRead(object):
    __slots__ = ["qname", "flag", "rname", "pos", "mapq", "cigar",
                 "tlen", "seq", "qual"]
    
    MIN_FIELDS = 11

    def __init__(self, line: bytes):
        fields = line.rstrip().split(b"\t")
        if len(fields) < self.MIN_FIELDS:
            raise ValueError(f"Invalid SAM line:\n{line}")
        self.qname = fields[0]
        self.flag = SamFlag(int(fields[1]))
        self.rname = fields[2]
        self.pos = int(fields[3])
        self.mapq = int(fields[4])
        self.cigar = fields[5]
        #self.rnext = fields[6]
        #self.pnext = int(fields[7])
        #self.tlen = int(fields[8])
        self.seq = fields[9]
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
    read_start_idx = 0
    read_end_idx = 0
    # position at which the current CIGAR operation starts
    # 0-indexed from beginning of region
    op_start_idx = read.pos - region_first
    # position at which the current CIGAR operation ends (0-indexed)
    # does not include the last position of the operation (like Python slicing)
    # 0-indexed from beginning of region
    op_end_idx = op_start_idx
    # Number of bases truncated from the end of the operation.
    truncated = 0
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
        if op_consumes_read(cigar_op):
            read_end_idx += op_length
        if op_end_idx > 0 and op_start_idx < region_length:
            # Run this block once the operation has entered the region.
            if op_start_idx < 0:
                # If the current operation starts before the region,
                # then advance it to the beginning of the region.
                if op_consumes_read(cigar_op):
                    read_start_idx -= op_start_idx
                op_length += op_start_idx
                op_start_idx = 0
            if op_end_idx > region_length:
                # If the current operation ends after the region,
                # then truncate it to the end of the region.
                truncated = op_end_idx - region_length
                op_length -= truncated
                if op_consumes_read(cigar_op):
                    read_end_idx -= truncated
                op_end_idx = region_length
            # Perform an action based on the CIGAR operation and its length.
            if cigar_op == CIG_MAT:
                # Condition: read matches reference
                muts.extend(map(encode_match,
                                read.seq[read_start_idx: read_end_idx],
                                read.qual[read_start_idx: read_end_idx]))
            elif cigar_op == CIG_ALN or cigar_op == CIG_SUB:
                # Condition: read has a match or substitution relative to ref
                muts.extend(map(encode_compare,
                                region_seq[len(muts): len(muts) + op_length],
                                read.seq[read_start_idx: read_end_idx],
                                read.qual[read_start_idx: read_end_idx]))
            elif cigar_op == CIG_DEL:
                # Condition: read contains a deletion w.r.t. the reference.
                dels.extend(Deletion(ref_idx, read_start_idx) for ref_idx
                            in range(len(muts), len(muts) + op_length))
                muts.extend(DELET * op_length)
            elif cigar_op == CIG_INS:
                # Condition: read contains an insertion w.r.t. the reference.
                # Position added to insertions is of the base 3' of the insert.
                inns.extend(Insertion(idx, len(muts)) for idx
                            in range(read_start_idx, read_end_idx))
                # Insertions do not consume the reference, so do not add any
                # information to muts yet. That information is added later.
            elif cigar_op == CIG_SCL:
                # Condition: read contains a soft clipping
                pass
            else:
                raise ValueError(
                    f"Invalid CIGAR operation: '{cigar_op.decode()}'")
        # Advance the start positions to the end of the current operation.
        if truncated:
            op_end_idx += truncated
            if op_consumes_read(cigar_op):
                read_end_idx += truncated
            truncated = 0
        op_start_idx = op_end_idx
        read_start_idx = read_end_idx
    # Pad the end of the mutation vector with any non-covered positions.
    muts.extend(BLANK * (region_length - len(muts)))
    assert len(muts) == region_length
    # Ensure the CIGAR string matched the length of the read.
    if read_start_idx != len(read):
        raise ValueError(
            f"CIGAR string '{read.cigar.decode()}' consumed {read_start_idx} "
            f"bases from read, but read is {len(read)} bases long.")
    # Add insertions to muts.
    for ins in inns:
        ins.stamp(muts)
    # Label all positions that are ambiguous due to indels.
    if dels or inns:
        allindel(muts, region_seq, bytes(read.seq), read.qual, dels, inns)
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
                assert read2.flag.paired
                assert read1.qname == read2.qname
                assert read1.rname == read2.rname
                assert read1.flag.first and read2.flag.second
                assert read1.flag.rev != read2.flag.rev
            else:
                raise ValueError("read1 is unpaired, but read2 was given")
        self.read2 = read2
    
    @property
    def read_name(self):
        return self.read1.qname.decode()

    @property
    def ref_name(self):
        return self.read1.rname.decode()

    @property
    def paired(self):
        return self.read1.flag.paired

    def vectorize(self, region_seq: bytes, region_start: int, region_end: int):
        if self.read2 is None:
            return vectorize_read(region_seq, region_start, region_end,
                                  self.read1)
        else:
            return vectorize_pair(region_seq, region_start, region_end,
                                  self.read1, self.read2)
