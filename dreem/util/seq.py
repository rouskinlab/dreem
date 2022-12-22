# Constants

BASES = b"ACGT"
COMPS = b"TGCA"
RBASE = b"ACGU"
RCOMP = b"UGCA"
BASEN = b"N"
BASES_SET = set(BASES)
BASEN_SET = set(BASES + BASEN)



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
        """
        Transcribe DNA into RNA.
        """
        return RNA(self.replace(b"T", b"U"))


class RNA(Seq):
    alph = RBASE
    comp = RCOMP
    alphaset = set(alph)
    trans = alph.maketrans(alph, comp)

    def rt(self):
        """
        Reverse transcribe RNA into DNA.
        """
        return DNA(self.replace(b"U", b"T"))
