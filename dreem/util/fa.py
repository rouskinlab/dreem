from typing import Dict, Set

from dreem.util.seq import DNA


class FastaIO(object):
    __slots__ = ["_path", "_refs"]

    defsymbol = b">"
    deftrunc = len(defsymbol)

    def __init__(self, path: str):
        self._path = path


class FastaParser(FastaIO):
    def __init__(self, path: str):
        super().__init__(path)
        self._refs: Set[bytes] = set()

    @classmethod
    def _parse_fasta_record(cls, fasta, line: bytes):
        if not line.startswith(cls.defsymbol):
            raise ValueError("FASTA definition line does not start with "
                             f"'{cls.defsymbol.decode()}'")
        name = line.split()[0][cls.deftrunc:]
        seq = bytearray()
        while (line := fasta.readline()) and not line.startswith(cls.defsymbol):
            seq.extend(line.rstrip())
        seq = DNA(bytes(seq))
        return line, name, seq

    def parse(self):
        with open(self._path, "rb") as f:
            line = f.readline()
            while line:
                line, name, seq = self._parse_fasta_record(f, line)
                if name in self._refs:
                    raise ValueError(
                        f"Duplicate entry in {self._path}: '{name.decode()}'")
                self._refs.add(name)
                yield name, seq


class FastaWriter(FastaIO):
    def __init__(self, path: str, refs: Dict[bytes, DNA]):
        super().__init__(path)
        self._refs = refs
    
    def write(self):
        with open(self._path, "wb") as f:
            for ref, seq in self._refs.items():
                f.write(b"".join((self.defsymbol, ref, b"\n", seq, b"\n")))
                

def parse_fasta(path: str):
    return FastaParser(path).parse()
