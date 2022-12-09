

class Construct:
    """Container object. Contains the name of the construct, the sequence, the bitvector, the read names and the read count.
    """
    
    def __init__(self, name, sequence, bv, read_names, read_count) -> None:
        self.name = name
        self.sequence = sequence
        self.bv = bv
        self.read_names = read_names
        self.read_count = read_count
        
        