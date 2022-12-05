import yaml, sys, os, random
import pandas as pd
import numpy as np
import subprocess
import json

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

def generate_fastq_files(fastq1_name, fastq2_name, sample_profile):
    """Write a fastq file with the given parameters
    
    Arguments:
        fastq1_name {str} -- where to write the fastq file for read 1
        fastq2_name {str} -- where to write the fastq file for read 2
        sample_profile {dict} -- dictionary with the following keys
            'constructs' -- list of construct names [list]
                'sequences' -- sequence to use for each construct [list]
                'number_of_reads' -- number of reads to generate for each construct [list]
                'mutations' -- number of mutations to introduce in each read [list]
                'deletions' -- number of deletions to introduce in each read [list]
                'insertions' -- number of insertions to introduce in each read [list]
    """
    with open(fastq1_name, 'w') as f1, open(fastq2_name, 'w') as f2:
        # write placeholder reads
        for c, v in sample_profile.items():
            for i in range(v['number_of_reads']):
                sequence = v['sequences'][i]
                for j in v['mutations'][i]:
                    sequence = sequence[:j] + next_base(sequence[j]) + sequence[j+1:]
                for j in v['insertions'][i]:
                    sequence = sequence[:j] + create_sequence(1) + sequence[j:]
                for j in v['deletions'][i]:
                    sequence = sequence[:j] + sequence[j+1:]
                print_fastq_line(f1, '{}:{}'.format(c, i), sequence, 'F'*len(sequence))
                print_fastq_line(f2, '{}:{}'.format(c, i), invert_sequence(sequence), 'F'*len(sequence))
    
def generate_fasta_file(filename, sample_profile):
    """Write a fasta file with the given parameters
    
    Arguments:
        filename {str} -- where to write the fasta file
        sample_profile {dict} -- dictionary with the following keys
            'constructs' -- list of construct names [list]
                'sequences' -- sequence to use for each construct [list]
    """
    with open(filename, 'w') as f:
        for c, v in sample_profile.items():
            print_fasta_line(f, c, v['sequences'])
            
def generate_library_file(filename, sample_profile):
    """Write a library file with the given parameters

    Arguments:
        filename {str} -- where to write the library file
        sample_profile {dict} -- dictionary with the following keys
            constructs {list} -- list of construct names [list]
                barcodes {list} -- list of barcodes to use for each construct [list]
                barcode_start {int} -- where to start the barcode in the read [int]
    """

    df = pd.DataFrame(sample_profile).T.reset_index().rename(columns={'index':'construct'})
    content_cols = ['construct']
    if 'barcodes' in df.columns:
        df.rename(columns={'barcodes':'barcode'}, inplace=True)
        content_cols += ['barcode_start','barcode']
        df = df.explode(['barcode'])
    if 'sections' in df.columns:
        df.rename(columns={'sections':'section'}, inplace=True)
        content_cols += ['section','section_start','section_end']
        df = df.explode(['section','section_start','section_end'])[content_cols]
    df.to_csv(filename, index=False)

def generate_demultiplexed_fastq(folder, sample_profile):
    """Demultiplex a fastq file based on a library file

    Arguments:
        folder {str} -- where to write the demultiplexed fastq files
        sample_profile {dict} -- dictionary with the following keys
            'constructs' -- list of construct names [list]
                'sequences' -- sequence to use for each construct [list]
                'number_of_reads' -- number of reads to generate for each construct [list]
                'mutations' -- number of mutations to introduce in each read [list]
                'deletions' -- number of deletions to introduce in each read [list]
                'insertions' -- number of insertions to introduce in each read [list]
    """

    for c, v in sample_profile.items():
        generate_fastq_files(os.path.join(folder, '{}_R1.fastq'.format(c)),  os.path.join(folder, '{}_R2.fastq'.format(c)), {c:v})


def print_sam_header(f, construct, len_sequence):
    f.write('@HD VN:1.0 SO:unsorted\n' + '@SQ SN:{} LN:{}\n'.format(construct, len_sequence))

def print_sam_lines(f, read_name, sequence, construct, cigar):
    f.write('{}\t3\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}\n'.format(read_name, construct, 1, 255, cigar, len(sequence), sequence))
    f.write('{}\t19\t{}\t{}\t{}\t{}\t*\t0\t0\t{}\t{}\n'.format(read_name, construct, 1, 255, cigar, len(sequence), invert_sequence(sequence)))

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

    return cigar


def generate_sam_files(folder, sample_profile):
    """Generate sam files for the given parameters. 
    The sam files are written from the input parameters, not from an external software. They do not use fastq files.

    Arguments:
        folder {str} -- where to write the sam files
        sample_profile {dict} -- dictionary of parameters for each construct
            - constructs {list} -- list of construct names [list]
                - number_of_reads {list} -- number of reads to generate for each construct [list]
                - mutations {list} -- number of mutations to introduce in each read for each construct [list(list)]
                - sequences {list} -- sequence to use for each construct [list]
                - insertions {list} -- number of insertions to introduce in each read for each construct [list(list)]
                - deletions {list} -- number of deletions to introduce in each read for each construct [list(list)]
                - barcode_start {int} -- where to start the barcode in the read [int]
                - barcodes {list} -- list of barcodes to use for each construct [list]            
    """

    for c, v in sample_profile.items():
        with open(os.path.join(folder, '{}.sam'.format(c)), 'w') as f:
            print_sam_header(f, c, len(v['sequences']))
            for i in range(v['number_of_reads']):
                cigar = make_cigar(len( v['sequences'][i]), v['mutations'][i], v['insertions'][i], v['deletions'][i])
                print_sam_lines(f, '{}:{}'.format(c, str(i).zfill(len(str(v['number_of_reads'])))), v['sequences'][i], c, cigar)


def check_sample_profile(sequences, mutations, insertions, deletions):
    """ Check that the sample profile is valid.
    Make sure that:
    - the mutations, insertions, and deletions are sorted
    - the mutations, insertions, and deletions are not consecutive
    - the mutations, insertions, and deletions are not at the same position
    - the mutations, insertions, and deletions are within the sequence length    
    """
    for idx, s in enumerate(sequences):
        [mutations[idx][i].sort() for i in range(len(mutations[idx]))]
        [insertions[idx][i].sort() for i in range(len(insertions[idx]))]
        [deletions[idx][i].sort() for i in range(len(deletions[idx]))]
        for j in range(len(mutations[idx])):
            for i in range(len(mutations[idx][j])-1):
                assert mutations[idx][j][i]+1 != mutations[idx][j][i+1], 'Consecutive mutations are not allowed'
            for i in range(len(insertions[idx][j])-1):
                assert insertions[idx][j][i]+1 != insertions[idx][j][i+1], 'Consecutive insertions are not allowed'
            for i in range(len(deletions[idx][j])-1):
                assert deletions[idx][j][i]+1 != deletions[idx][j][i+1], 'Consecutive deletions are not allowed'
            for i in range(len(mutations[idx][j])):
                assert mutations[idx][j][i] < len(s[j]), 'Mutation is out of sequence length'
            for i in range(len(insertions[idx][j])):
                assert insertions[idx][j][i] < len(s[j]), 'Insertion is out of sequence length'
            for i in range(len(deletions[idx][j])):
                assert deletions[idx][j][i] < len(s[j]), 'Deletion is out of sequence length'
        for j in range(len(mutations[idx])):
            for i in range(len(insertions[idx][j])):
                assert mutations[idx][j][i] != insertions[idx][j][i], 'Mutation and insertion are at the same position'
            for i in range(len(deletions[idx][j])):
                assert mutations[idx][j][i] != deletions[idx][j][i], 'Mutation and deletion are at the same position'
            for i in range(len(insertions[idx][j])):
                assert insertions[idx][j][i] != deletions[idx][j][i], 'Insertion and deletion are at the same position'
        

def make_sample_profile(constructs, sequences, number_of_reads, mutations, insertions, deletions, barcodes=None, barcode_start=None, sections=None, section_start=None, section_end=None):
    """Make a dictionary of the sample profile.
    The dictionary is used to make sure that the sample profile is valid. 
    The lists are grouped by construct, and the number of reads is grouped by construct.
    """
    check_sample_profile(sequences, mutations, insertions, deletions)
    if barcodes is not None:
        assert len(constructs) == len(barcodes)
        for b in barcodes:
            assert len(b) == len(barcodes[0])
    
    sample_profile = {}
               
    for idx, c in enumerate(constructs):
        sample_profile[c] = {}
        sample_profile[c]['reference'] = sequences[idx][0]            
        sample_profile[c]['sequences'] = sequences[idx]            
        sample_profile[c]['number_of_reads'] = number_of_reads[idx]
        sample_profile[c]['mutations'] = mutations[idx]
        sample_profile[c]['insertions'] = insertions[idx]
        sample_profile[c]['deletions'] = deletions[idx]
        if sections is not None:
            sample_profile[c]['sections'] = sections[idx]
            sample_profile[c]['section_start'] = section_start[idx]
            sample_profile[c]['section_end'] = section_end[idx]
        if barcodes is not None:
            sample_profile[c]['barcodes'] = barcodes[idx]
            sample_profile[c]['barcode_start'] = barcode_start
            for j in range(len(sample_profile[c]['sequences'])): 
                sample_profile[c]['sequences'][j] = sample_profile[c]['sequences'][j][:barcode_start] + barcodes[idx][j] + sample_profile[c]['sequences'][j][barcode_start:]
        for i in range(sample_profile[c]['number_of_reads']):
            sequence = sample_profile[c]['sequences'][i]
            for j in sample_profile[c]['mutations'][i]:
                sequence = sequence[:j] + next_base(sequence[j]) + sequence[j+1:]
            for j in sample_profile[c]['insertions'][i]:
                sequence = sequence[:j] + create_sequence(1) + sequence[j:]
            for j in sample_profile[c]['deletions'][i]:
                sequence = sequence[:j] + sequence[j+1:]
            sample_profile[c]['sequences'][i] = sequence
    return sample_profile



def run_notebook(notebook_path):
    import dreem.util as util
    with open(notebook_path, 'r') as f:
        cells = json.load(f)['cells']
        for cell in cells:
            if cell['cell_type'] != 'markdown':
                cell['source'] = [line for line in cell['source'] if not line.startswith('%')]
                lines = ''.join(cell['source'])
                exec(lines)


def generate_bitvector_files(folder, sample_profile, library):
    """Generate a bitvector file for the sample profile. 
    
    The bitvector is an orc file generated using pandas. A bitvector is a binary file that contains the reads for a certain section of a certain construct. The sections of the construct are defined by the library. 
    For each section of each construct, we are interested in the residues inside the [section_start, section_end] 0-indexed interval. 
    The bitvector columns are the residues, and the rows are the reads. The columns are under the format '0-indexed position'+'residue'. Ex: '0A' means that the first residue is an A.

    Parameters
    ----------
    folder : str
        The folder to write the bitvector file to.
    sample_profile : dict
        The sample dictionary.
        Attributes are:
        - constructs
            - sequences
            - number_of_reads
            - mutations
            - insertions
            - deletions
    library : str
        The library to use. Columns are 'construct', 'section', 'section_start','section_end'.

    """

    library = pd.read_csv(library)

    # get through the constructs
    for construct in sample_profile:

    # create a folder under the construct's name
        construct_folder = os.path.join(folder, construct)
        if not os.path.exists(construct_folder):
            os.makedirs(construct_folder)
    # get through the sections for that construct
        for section in library[library['construct'] == construct]['section']:
            section_start = library[(library['construct'] == construct) & (library['section'] == section)]['section_start'].values[0]
            section_end = library[(library['construct'] == construct) & (library['section'] == section)]['section_end'].values[0]
            sequence = sample_profile[construct]['reference'][section_start:section_end]
            assert len(sequence) == section_end - section_start, 'Section length is not equal to the length of the sequence'
            columns = [str(i)+sequence[i] for i in range(section_end - section_start)]

            # create a dataframe with the columns
            df = pd.DataFrame(columns=columns)
            for i in range(sample_profile[construct]['number_of_reads']):
                bv = [1]*len(columns)
                for j in range(len(columns)):
                    # if the residue has a deletion, use update_bv_byte
                    if j+section_start in sample_profile[construct]['deletions'][i]:
                        bv[j] = update_bv_byte(bv[j], 'deletion')
                    if j+section_start in sample_profile[construct]['insertions'][i]:
                        bv[j] = update_bv_byte(bv[j], 'insertion_3')
                        bv[j+1] = update_bv_byte(bv[j+1], 'insertion_5')
                    base = sample_profile[construct]['sequences'][i][j+section_start]
                    if base != sequence[j]:
                        bv[j] = update_bv_byte(bv[j], 'substitution_'+base)
                df.loc[i] = ''.join([str(b) for b in bv])
            df.to_csv(os.path.join(construct_folder, section+'.orc'), index=False)

def update_bv_byte(byte, position):
    """Add a bit to a byte and return the new byte.
    00000001: match
    00000010: deletion
    00000100: ≥1 base is inserted immediately 3' of this position
    00001000: ≥1 base is inserted immediately 5' of this position
    00010000: substitution to A
    00100000: substitution to C
    01000000: substitution to G
    10000000: substitution to T
    """
    if position == 'match':
        return byte + 2**0
    elif position == 'deletion':
        return byte + 2**1
    elif position == 'insertion_3':
        return byte + 2**2
    elif position == 'insertion_5':
        return byte + 2**3
    elif position == 'substitution_A':
        return byte + 2**4
    elif position == 'substitution_C':
        return byte + 2**5
    elif position == 'substitution_G':
        return byte + 2**6
    elif position == 'substitution_T':
        return byte + 2**7
    else:
        raise ValueError('Position not recognized')            


