from functools import cached_property
import logging
from pathlib import Path

from ..util import path


logger = logging.getLogger(__name__)

# Byte encodings for nucleic acid alphabets
BASES = b"ACGT"
COMPS = b"TGCA"
RBASE = b"ACGU"
RCOMP = b"UGCA"
BASEN = b"N"
BASES_SET = set(BASES)
BASEN_SET = set(BASES + BASEN)

# Integer encodings for nucleic acid alphabets
A_INT = BASES[0]
C_INT = BASES[1]
G_INT = BASES[2]
T_INT = BASES[3]
N_INT = BASEN[0]

# Integer encodings for mutation vectors
BLANK = b"\x00"[0]  # 00000000 (000): no coverage at this position
MATCH = b"\x01"[0]  # 00000001 (001): match with reference
DELET = b"\x02"[0]  # 00000010 (002): deletion from reference
INS_5 = b"\x04"[0]  # 00000100 (004): insertion 5' of base in reference
INS_3 = b"\x08"[0]  # 00001000 (008): insertion 3' of base in reference
SUB_A = b"\x10"[0]  # 00010000 (016): substitution to A
SUB_C = b"\x20"[0]  # 00100000 (032): substitution to C
SUB_G = b"\x40"[0]  # 01000000 (064): substitution to G
SUB_T = b"\x80"[0]  # 10000000 (128): substitution to T
SUB_N = SUB_A | SUB_C | SUB_G | SUB_T
ANY_N = SUB_N | MATCH
INDEL = DELET | INS_5 | INS_3
EVERY = ANY_N | INDEL


def get_diffs(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences were different lengths: "
                         f"'{seq1}' and '{seq2}'")
    diffs = [i for i, (x1, x2) in enumerate(zip(seq1, seq2)) if x1 != x2]
    return diffs


class Seq(bytes):
    __slots__ = []

    alph = b""
    comp = b""
    alphaset = set(alph)
    trans = alph.maketrans(alph, comp)

    def __init__(self, seq: bytes):
        self.validate_seq(seq)
        super().__init__()
    
    @classmethod
    def validate_seq(cls, seq):
        if not seq:
            raise ValueError("seq is empty")
        if set(seq) - cls.alphaset:
            raise ValueError(f"Invalid characters in seq: '{seq.decode()}'")

    @property
    def rc(self):
        return self.__class__(self[::-1].translate(self.trans))

    def __getitem__(self, item):
        return self.__class__(super().__getitem__(item))

    def __str__(self):
        return self.decode()


class DNA(Seq):
    alph = BASES
    comp = COMPS
    alphaset = set(alph)
    trans = alph.maketrans(alph, comp)

    def tr(self):
        """ Transcribe DNA to RNA. """
        return RNA(self.replace(b"T", b"U"))


class RNA(Seq):
    alph = RBASE
    comp = RCOMP
    alphaset = set(alph)
    trans = alph.maketrans(alph, comp)

    def rt(self):
        """ Reverse transcribe RNA to DNA. """
        return DNA(self.replace(b"U", b"T"))


class FastaIO(object):
    recsym = b">"
    deftrunc = len(recsym)

    def __init__(self, file: str | Path | None):
        if file is None:
            # If the user forgets to give a FASTA file when one is
            # required, then the program will crash with an ugly error
            # when it tries to open a path that is None. This check
            # raises an error that describes the specific problem.
            raise TypeError("No FASTA file was given.")
        self._path = str(file)

    @cached_property
    def refset(self):
        """ Return the name of the set of references. """
        return path.RefsetSeqInFilePath.parse(self._path).refset


class FastaParser(FastaIO):
    def __init__(self, file: str | Path):
        super().__init__(file)
        self._refs: set[str] = set()

    @classmethod
    def _parse_fasta_record(cls, fasta, line: bytes):
        if not line.startswith(cls.recsym):
            raise ValueError("FASTA definition line does not start with "
                             f"'{cls.recsym.decode()}'")
        name = line.split()[0][cls.deftrunc:].decode()
        seq = bytearray()
        while (line := fasta.readline()) and not line.startswith(cls.recsym):
            seq.extend(line.rstrip())
        seq = DNA(bytes(seq))
        return line, name, seq

    def parse(self):
        has_ref_named_refset = False
        with open(self._path, "rb") as f:
            line = f.readline()
            while line:
                # Read the name from the current line and the sequence
                # from the next line, and then load the next, next line.
                line, name, seq = self._parse_fasta_record(f, line)
                if name in self._refs:
                    # If there are two or more references with the same
                    # name, then the sequence of only the first is used.
                    logging.warning(f"Duplicate ref '{name}' in {self._path}")
                    continue
                # Record the name of every reference in the file.
                self._refs.add(name)
                # If a reference has the same name as the file, then it
                # must be the only reference in the file. Otherwise, if
                # the FASTA is split into individual references during
                # demultiplexing, then it would be overwritten by the
                # FASTA containing the reference with the same name.
                # If that reference is the only one in the file, though,
                # then the file after splitting would be identical to
                # the file before splitting.
                if name == self.refset:
                    # Whether a reference has the same name as the file.
                    has_ref_named_refset = True
                if has_ref_named_refset and len(self._refs) > 1:
                    raise ValueError(
                        f"Because {self._path} includes a reference whose name "
                        f"matches the file name ({name}), that reference must "
                        f"be the only one in the file, but also got {self._refs}")
                yield name, seq


class FastaWriter(FastaIO):
    def __init__(self, file: str | Path, refs: dict[str, DNA]):
        super().__init__(file)
        self._refs = refs
    
    def write(self):
        with open(self._path, "wb") as f:
            for ref, seq in self._refs.items():
                f.write(b"".join(
                    [self.recsym, ref.encode(), b"\n", seq, b"\n"]
                ))
                

def parse_fasta(file: str):
    return FastaParser(file).parse()
