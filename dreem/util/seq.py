from logging import getLogger
from pathlib import Path
from typing import Iterable

from ..util import path


logger = getLogger(__name__)

# FASTA format
FASTA_NAMESYM = b">"

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
    comptrans = alph.maketrans(alph, comp)

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
        return self.__class__(self[::-1].translate(self.comptrans))

    def __getitem__(self, item):
        return self.__class__(super().__getitem__(item))

    def __str__(self):
        return self.decode()


class DNA(Seq):
    alph = BASES
    comp = COMPS
    alphaset = set(alph)
    comptrans = alph.maketrans(alph, comp)

    def tr(self):
        """ Transcribe DNA to RNA. """
        return RNA(self.replace(b"T", b"U"))


class RNA(Seq):
    alph = RBASE
    comp = RCOMP
    alphaset = set(alph)
    comptrans = alph.maketrans(alph, comp)

    def rt(self):
        """ Reverse transcribe RNA to DNA. """
        return DNA(self.replace(b"U", b"T"))


def parse_fasta(fasta: str | Path):
    """ Parse a FASTA file and iterate through the reference names and
    sequences. """
    if not fasta:
        raise TypeError("No FASTA file given")
    logger.info(f"Began parsing FASTA: {fasta}")
    # Get the name of the set of references.
    refset = path.RefsetSeqInFilePath.parse(fasta).refset
    has_ref_named_refset = False
    # Record the names of all the references.
    names = set()
    with open(fasta, "rb") as f:
        line = f.readline()
        while line:
            # Read the name from the current line.
            if not line.startswith(FASTA_NAMESYM):
                logger.error(f"Name line '{line.strip()}' in {fasta} does not "
                             f"start with name symbol '{FASTA_NAMESYM}'")
                continue
            # Get the name of the reference up to the first whitespace.
            name = line.split(maxsplit=1)[0][len(FASTA_NAMESYM):].decode()
            # Read the sequence of the reference up until the next
            # reference or the end of the file, whichever comes first.
            seqarray = bytearray()
            while (line := f.readline()) and not line.startswith(FASTA_NAMESYM):
                seqarray.extend(line.rstrip())
            # Confirm that the sequence is a valid DNA sequence.
            try:
                seq = DNA(seqarray)
            except Exception as error:
                logger.error(f"Failed to parse sequence in {fasta}: {error}")
                continue
            # Confirm that the name is not blank.
            if not name:
                logger.error(f"Blank name line '{line.strip()}' in {fasta}")
                continue
            # If there are two or more references with the same name,
            # then the sequence of only the first is used.
            if name in names:
                logger.warning(f"Duplicate reference '{name}' in {fasta}")
                continue
            # If any reference has the same name as the file, then the
            # file is not allowed to have any additional references
            # because, if it did, then the files of all references and
            # of only the self-named reference would have the same names
            # and thus be indistinguishable by their paths.
            if name == refset:
                has_ref_named_refset = True
                logger.debug(f"Reference '{name}' had same name as {fasta}")
            if has_ref_named_refset and names:
                raise ValueError(f"Because {fasta} had a reference with the "
                                 f"same name as the file ('{name}'), it was "
                                 f"not allowed to have any other references, "
                                 f"but it also had {', '.join(names)}")
            # Yield the validated name and sequence.
            names.add(name)
            logger.debug(f"Read reference '{name}' of length {len(seq)} "
                         f"from {fasta}")
            yield name, seq
    logger.info(f"Ended parsing {len(names)} references from FASTA: {fasta}")


def write_fasta(fasta: str | Path, refs: Iterable[tuple[str, DNA]]):
    """ Write an iterable of reference names and DNA sequences to a
    FASTA file. """
    if not fasta:
        raise TypeError("No FASTA file given")
    logger.info(f"Began writing FASTA file: {fasta}")
    # Get the name of the set of references.
    refset = path.RefsetSeqInFilePath.parse(fasta).refset
    has_ref_named_refset = False
    # Record the names of all the references.
    names = set()
    with open(fasta, "xb") as f:
        for name, seq in refs:
            # Confirm that the name is not blank.
            if not name:
                logger.error(f"Blank reference name")
                continue
            # If there are two or more references with the same name,
            # then the sequence of only the first is used.
            if name in names:
                logger.warning(f"Duplicate reference '{name}'")
                continue
            # If any reference has the same name as the file, then the
            # file is not allowed to have any additional references
            # because, if it did, then the files of all references and
            # of only the self-named reference would have the same names
            # and thus be indistinguishable by their paths.
            if name == refset:
                has_ref_named_refset = True
                logger.debug(f"Reference '{name}' had same name as {fasta}")
            if has_ref_named_refset and names:
                raise ValueError(f"Because {fasta} got a reference with the "
                                 f"same name as the file ('{name}'), it was "
                                 f"not allowed to get any other references, "
                                 f"but it also got {', '.join(names)}")
            try:
                f.write(b"".join([FASTA_NAMESYM, name.encode(), b"\n",
                                  seq, b"\n"]))
            except Exception as error:
                logger.error(
                    f"Error writing reference '{name}' to {fasta}: {error}")
            else:
                logger.debug(f"Wrote reference '{name}' of length {len(seq)}"
                             f"to {fasta}")
                names.add(name)
    logger.info(f"Wrote {len(names)} references to {fasta}")
