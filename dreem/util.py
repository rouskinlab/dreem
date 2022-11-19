
def get_folders(args):
  root_dir = args['root_dir'] if 'root_dir' in args else ''
  root_dir = root_dir if root_dir[-1] == '/' else root_dir + '/'
  
  return {
    'demultiplexing':{
      'temp': root_dir +'temp/demultiplexing/',
      'output': root_dir +'output/demultiplexing/'
    },
    'alignment':{
      'temp': root_dir +'temp/alignment/',
      'output': root_dir +'output/alignment/'
    },
    'vectoring':{
      'temp': root_dir +'temp/vectoring/',
      'output': root_dir +'output/vectoring/'
    },
    'clustering':{
      'temp': root_dir +'temp/clustering/',
      'output': root_dir +'output/clustering/'
    },
    'aggregate':{
      'temp': root_dir +'temp/aggregate/',
      'output': root_dir +'output/aggregate/'
    },
    'post_processing':{
      'temp': root_dir +'temp/post_processing/',
      'output': root_dir +'output/post_processing/'
    }}



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

    