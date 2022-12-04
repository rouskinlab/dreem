import yaml, sys, os, random

def make_folder(folder):
    if not os.path.exists(folder):
        os.makedirs(folder)
    return folder

def clear_folder(folder):
    os.system('rm -fr ' + folder)
    os.makedirs(folder)

def run_cmd(cmd):
    return os.system(cmd), cmd

def make_cmd(args, module):
    cmd = 'dreem-' + module + ' '
    for key, value in args.items():
        if type(value) in (list, tuple):
            for v in value:
                cmd += '--' + key + ' ' + str(v) + ' '
        elif value is not None:
            cmd += '--' + key + ' ' + str(value) + ' '
    return cmd

class Seq(bytes):
    __slots__ = []

    alph = b""
    comp = b""
    low_qual = b"N"

    def __init__(self, seq: bytes):
        if not seq:
            raise ValueError("seq is empty")
        if any(base not in self.alph for base in seq):
            raise ValueError(f"Invalid characters in seq: '{seq.decode()}'")
        super().__init__()

    @property
    def rc(self):
        return self.__class__(
            self[::-1].translate(self.maketrans(self.alph, self.comp)))

    def __getitem__(self, item):
        return self.__class__(super().__getitem__(item))

    def __str__(self):
        return self.decode()


class DNA(Seq):
    alph = b"ACGT"
    comp = b"TGCA"

    @property
    def tr(self):
        return RNA(self.replace(b"T", b"U"))


class RNA(Seq):
    alph = b"ACGU"
    comp = b"UGCA"

    @property
    def rt(self):
        return DNA(self.replace(b"U", b"T"))


class Primer(str):
    @property
    def as_dna(self):
        return DNA(self.encode())


def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def generate_barcodes(barcode_length, n, min_barcode_distance):
    barcodes = []
    while(len(barcodes) < n):
        barcode = ''.join(random.choice('ATCG') for _ in range(barcode_length))
        if all(hamming_distance(barcode, b) > min_barcode_distance for b in barcodes):
            barcodes.append(barcode)
    return barcodes

def invert_sequence(seq):
    return ''.join([{'A':'T','T':'A','C':'G','G':'C'}[s] for s in seq])[::-1]

def next_base(base):
    return {'A':'T','T':'C','C':'G','G':'A',0:1}[base]

def create_sequence(length, bases=['A','T','C','G']):
    return ''.join([random.choice(bases) for _ in range(length)])



class FastQ:
    def __init__(self, filename:str, number_of_reads:list, mutations, sequence, insertions, deletions, barcode_start, barcodes, constructs):
        self.filename = filename  # where to write the fastq files
        self.number_of_reads = number_of_reads # number of reads to generate for each construct [list]
        self.mutations = mutations # number of mutations to introduce in each read for each construct [list(list)]
        self.sequence = sequence # sequence to use for all construct str
        self.insertions = insertions # number of insertions to introduce in each read for each construct [list(list)]
        self.deletions = deletions # number of deletions to introduce in each read for each construct [list(list)]
        self.barcode_start = barcode_start # where to start the barcode in the read [int]
        self.barcodes = barcodes  # list of barcodes to use for each construct [list]
        self.constructs = constructs # list of construct names [list]
        assert len(self.barcodes) == len(self.constructs)
        assert len(self.mutations) == len(self.constructs)
        assert len(self.insertions) == len(self.constructs)
        assert len(self.deletions) == len(self.constructs)
    
    def write(self):
        pass # TODO
                    
