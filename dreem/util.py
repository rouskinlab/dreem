import yaml, sys, os, random
import pandas as pd
import subprocess

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

def print_fasta_line(f, id, seq):
    f.write('>{}\n{}\n'.format(id, seq))

def print_fastq_line(f, id, seq, qual):
    f.write('@{}\n{}\n+\n{}\n'.format(id, seq, qual))    

def generate_fastq_files(fastq1_name, fastq2_name, number_of_reads, mutations, sequences, insertions, deletions, barcode_start, barcodes, constructs):
    """Write a fastq file with the given parameters
    
    Arguments:
        fastq1_name {str} -- where to write the fastq file for read 1
        fastq2_name {str} -- where to write the fastq file for read 2
        number_of_reads {list} -- number of reads to generate for each construct [list]
        mutations {list} -- number of mutations to introduce in each read for each construct [list(list)]
        sequences {list} -- sequence to use for each construct [list]
        insertions {list} -- number of insertions to introduce in each read for each construct [list(list)]
        deletions {list} -- number of deletions to introduce in each read for each construct [list(list)]
        barcode_start {int} -- where to start the barcode in the read [int]
        barcodes {list} -- list of barcodes to use for each construct [list]
        constructs {list} -- list of construct names [list]
    """
    assert len(barcodes) == len(constructs)
    assert len(mutations) == len(constructs)
    assert len(insertions) == len(constructs)
    assert len(deletions) == len(constructs)
    assert len(sequences) == len(constructs)

    with open(fastq1_name, 'w') as f1, open(fastq2_name, 'w') as f2:
        # write placeholder reads
        for idx, c in enumerate(constructs):
            for i in range(number_of_reads[idx]):
                sequence = sequences[idx]
                sequence = sequence[:barcode_start] + barcodes[idx] + sequence[barcode_start+len(barcodes[idx]):]
                for j in mutations[idx][i]:
                    sequence = sequence[:j] + next_base(sequence[j]) + sequence[j+1:]
                for j in insertions[idx][i]:
                    sequence = sequence[:j] + create_sequence(1) + sequence[j:]
                for j in deletions[idx][i]:
                    sequence = sequence[:j] + sequence[j+1:]
                print_fastq_line(f1, '{}:{}'.format(c, i), sequence, 'F'*len(sequence))
                print_fastq_line(f2, '{}:{}'.format(c, i), invert_sequence(sequence), 'F'*len(sequence))
    
def generate_fasta_file(filename, sequences, constructs, barcodes, barcode_start):
    """Write a fasta file with the given parameters
    
    Arguments:
        filename {str} -- where to write the fasta file
        sequences {list} -- sequence to use for each construct [list]
        constructs {list} -- list of construct names [list]
        barcodes {list} -- list of barcodes to use for each construct [list]
        barcode_start {int} -- where to start the barcode in the read [int]
    """
    assert len(barcodes) == len(constructs)
    with open(filename, 'w') as f:
        for b, c, s in zip(barcodes, constructs, sequences):
            print_fasta_line(f, c, s[:barcode_start] + b + s[barcode_start+len(b):])
            
def generate_library_file(filename, constructs, barcodes, barcode_start):
    """Write a library file with the given parameters

    Arguments:
        filename {str} -- where to write the library file
        constructs {list} -- list of construct names [list]
        barcodes {list} -- list of barcodes to use for each construct [list]
        barcode_start {int} -- where to start the barcode in the read [int]
    """
    assert len(barcodes) == len(constructs), 'Number of barcodes and constructs must be the same'
    assert sum([len(b) == len(barcodes[0]) for b in barcodes]) == len(barcodes), 'All barcodes must be the same length'
    df = pd.DataFrame({'construct':constructs, 'barcode':barcodes, 'barcode_start':barcode_start, 'barcode_end':barcode_start+len(barcodes[0])})
    df.to_csv(filename, index=False)

def generate_demultiplexed_fastq(folder, number_of_reads, mutations, sequences, insertions, deletions, barcode_start, barcodes, constructs):
    """Demultiplex a fastq file based on a library file

    Arguments:
        folder {str} -- where to write the demultiplexed fastq files
        number_of_reads {list} -- number of reads to generate for each construct [list]
        mutations {list} -- number of mutations to introduce in each read for each construct [list(list)]
        sequences {list} -- sequence to use for each construct [list]
        insertions {list} -- number of insertions to introduce in each read for each construct [list(list)]
        deletions {list} -- number of deletions to introduce in each read for each construct [list(list)]
        barcode_start {int} -- where to start the barcode in the read [int]
        barcodes {list} -- list of barcodes to use for each construct [list]
        constructs {list} -- list of construct names [list]
    """

    for idx, c in enumerate(constructs):
        generate_fastq_files(os.path.join(folder, '{}_R1.fastq'.format(c)),  os.path.join(folder, '{}_R2.fastq'.format(c)), [number_of_reads[idx]], [mutations[idx]], [sequences[idx]], [insertions[idx]], [deletions[idx]], barcode_start, [barcodes[idx]], [c])


def print_sam_header(f, construct, len_sequence):
    f.write('@HD VN:1.0 SO:unsorted\n' + '@SQ SN:{} LN:{}\n'.format(construct, len_sequence))

def print_sam_lines(f, read_name, sequence, construct, cigar):
    f.write('{}\t3\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}\tNM:i:0\tAS:i:{}\tXS:i:0\n'.format(read_name, construct, 1, 255, cigar, len(sequence), sequence))
    f.write('{}\t19\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}\tNM:i:0\tAS:i:{}\tXS:i:0\n'.format(read_name, construct, 1, 255, cigar, len(sequence), invert_sequence(sequence)))

def make_cigar(len_sequence, mutations, insertions, deletions):
    """Create a cigar string for a read with the given mutations and indels
    The sequence is assumed to be 1-based.
    The cigar string is 1-based.
    The mutations and indels are 0-based.

    The mutations and indels are assumed to be sorted. This is checked at the beginning of the code.

    Arguments:
        len_sequence {int} -- length of the sequence [int]
        mutations {list} -- list of mutations [list]
        insertions {list} -- list of insertions [list]
        deletions {list} -- list of deletions [list]
    """

    # check that the mutations and indels are sorted
    assert sorted(mutations) == mutations
    assert sorted(insertions) == insertions
    assert sorted(deletions) == deletions

    for i in range(len(mutations)-1):
        assert mutations[i]+1 != mutations[i+1], 'Consecutive mutations are not allowed'
    
    for i in range(len(insertions)-1):
        assert insertions[i]+1 != insertions[i+1], 'Consecutive insertions are not allowed'
    
    for i in range(len(deletions)-1):
        assert deletions[i]+1 != deletions[i+1], 'Consecutive deletions are not allowed'

    # initialize the cigar string
    cigar = ''

    # keep track of the current position in the sequence
    current_position = 1

    # keep track of the current position in the cigar string
    current_cigar_position = 1
    
    # loop through the mutations, insertions, and deletions to create the cigar string. If there's no mutation, insertion, or deletion at a position, then the position is a match.
    while current_position <= len_sequence:
        # if there's a mutation at the current position, then add a mutation to the cigar string
        if current_position-1 in mutations:
            cigar += '{}M1X'.format(current_cigar_position)
            current_cigar_position = 0
            current_position += 1
        # if there's an insertion at the current position, then add an insertion to the cigar string
        elif current_position-1 in insertions:
            cigar += '{}M1I'.format(current_cigar_position)
            current_cigar_position = 0
            current_position += 1
        # if there's a deletion at the current position, then add a deletion to the cigar string
        elif current_position-1 in deletions:
            cigar += '{}M1D'.format(current_cigar_position)
            current_cigar_position = 0
            current_position += 1
        # if there's no mutation, insertion, or deletion at the current position, then add a match to the cigar string
        else:
            current_cigar_position += 1
            current_position += 1
    
    # add the last match to the cigar string
    cigar += '{}M'.format(current_cigar_position)
    
        
            


def generate_sam_files(folder, temp_folder, number_of_reads, mutations, sequences, insertions, deletions, barcode_start, barcodes, constructs):
    """Generate sam files for the given parameters. 
    The sam files are written from the input parameters, not from an external software. They do not use fastq files.

    Arguments:
        folder {str} -- where to write the sam files
        number_of_reads {list} -- number of reads to generate for each construct [list]
        mutations {list} -- number of mutations to introduce in each read for each construct [list(list)]
        sequences {list} -- sequence to use for each construct [list]
        insertions {list} -- number of insertions to introduce in each read for each construct [list(list)]
        deletions {list} -- number of deletions to introduce in each read for each construct [list(list)]
        barcode_start {int} -- where to start the barcode in the read [int]
        barcodes {list} -- list of barcodes to use for each construct [list]
        constructs {list} -- list of construct names [list]
    """
    assert len(barcodes) == len(constructs)
    for idx, c in enumerate(constructs):
        with open(os.path.join(folder, '{}.sam'.format(c)), 'w') as f:
            print_sam_header(f, c, len(sequences[idx]))
            for i in range(number_of_reads[idx]):
                sequence = sequences[idx]
                sequence = sequence[:barcode_start] + barcodes[idx] + sequence[barcode_start+len(barcodes[idx]):]
                for j in mutations[idx][i]:
                    sequence = sequence[:j] + next_base(sequence[j]) + sequence[j+1:]
                for j in insertions[idx][i]:
                    sequence = sequence[:j] + create_sequence(1) + sequence[j:]
                for j in deletions[idx][i]:
                    sequence = sequence[:j] + sequence[j+1:]
                cigar = make_cigar(len(sequence), mutations[idx][i], insertions[idx][i], deletions[idx][i])
                print_sam_lines(f, '{}:{}'.format(c, str(i).zfill(len(str(max(number_of_reads))))), sequence, 'F'*len(sequence), c, 1, '{}M'.format(len(sequence)), 0)

