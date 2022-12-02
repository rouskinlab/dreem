
import pandas as pd
import random
import os
import numpy as np
import dreem.util as util
## Just a library for generating testing files

SAMPLE_NAME= 'one_mut_per_read'
DIRECTORY = 'test/test_files/'+SAMPLE_NAME
util.clear_folder(DIRECTORY)
OUTPUT_FOLDER = 'test/output'
util.clear_folder(OUTPUT_FOLDER)
NUM_READS=1000
NUM_CONSTRUCTS=10
BARCODE_LENGTH=8 
MIN_BARCODE_DISTANCE=3
BARCODE_START=10
LENGTH_SEQUENCE=200
USER = 'Yves'
DATE = '08/15/1769'
EXP_ENV = 'in_vivo'
TEMPERATURE_K = 310
INC_TIME_TOT_SECS = 10
CELL_LINE = 'my_cell_line_#'
SEQUENCE = 'CCGAAAGTCAGCTTTCATACCTATCAGACGATTCGTATTGTGTAGGTGTCTGTACTCTTTGTAGAGGTGTGGCGAAGCACGACAACCCGAACCTAATCTCCCTTCGATTCAGAAGTGTAGATTTCCTGGCAGGCCGCAAGGTCCGGACACACTCGCCGTACAACGGCCGACTGTCCGTCTAAAGTATTTTCTCCGTGTCT'#''.join(random.choice('ATCG') for _ in range(LENGTH_SEQUENCE))
##########################################
BARCODES = ['AGGATTAG', 'TTGCCCCT', 'ACCTAAGT', 'AGATGGTC', 'GAAAGCCA', 'GTTCTTGC', 'CGTGCAAT', 'CCGGGCGG', 'TCTACGTA', 'TACGAGAC']#generate_barcodes(BARCODE_LENGTH, NUM_CONSTRUCTS, MIN_BARCODE_DISTANCE)
BARCODE_START, BARCODE_END = BARCODE_START, BARCODE_START + BARCODE_LENGTH
BARCODE_START_REV, BARCODE_END_REV = LENGTH_SEQUENCE - BARCODE_END, LENGTH_SEQUENCE - BARCODE_START
BARCODES_REV = [util.invert_sequence(b) for b in BARCODES]
SEQUENCE_REV = util.invert_sequence(SEQUENCE)
CONSTRUCTS = BARCODES
NUM_SECTIONS = 2
SECTION_1_NAME = 'section1'
SECTION_2_NAME = 'section2'
SECTION_1_START, SECTION_1_END = 0, BARCODE_START
SECTION_2_START, SECTION_2_END = BARCODE_END, LENGTH_SEQUENCE-1

MUTATION_THRESHOLD = {
    10: 40,     # 10% of the reads mutate at residue 40
    2: 90,       # 50% - 10% of the reads mutate at residue 90
    None: 140   # 100% - 50% of the reads mutate at residue 140
}

def read_name(construct, read, r, mut_0_idx):
    return 'READ:{}:{}:{}:{}'.format(construct, read, r, mut_0_idx)

def add_mutation(r, this_sequence, give_idx=False):
    for k,v in MUTATION_THRESHOLD.items():
        if k is None or r < NUM_READS//k:
            break
    if give_idx:
        return v, this_sequence[:v] + util.next_base(this_sequence[v]) + this_sequence[v+1:]
    return this_sequence[:v] + util.next_base(this_sequence[v]) + this_sequence[v+1:]


def build_sequence_from_barcode(barcode, sequence, barcode_start, barcode_end):
    return sequence[:barcode_start] + barcode + sequence[barcode_end:]

def get_sequence_from_read(read):
    if read==1:
        return BARCODES, SEQUENCE, BARCODE_START, BARCODE_END
    elif read==2:
        return BARCODES_REV, SEQUENCE_REV, BARCODE_START_REV, BARCODE_END_REV
    else:
        raise ValueError('read must be 1 or 2')

def print_sam_header(f, construct):
    f.write('@HD VN:1.0 SO:unsorted\n' + '@SQ SN:{} LN:{}\n'.format(construct, LENGTH_SEQUENCE))

def print_fastq_line(f, id, seq, qual):
    f.write('@{}\n{}\n+\n{}\n'.format(id, seq, qual))    

def print_fasta_line(f, id, seq):
    f.write('>{}\n{}\n'.format(id, seq))

def print_sam_line(f, construct, read, r, mut_0_idx, seq):
    flag = 1+2+16 if read==2 else 1+2
    pre_mut_bases_count, post_mut_bases_count = mut_0_idx, LENGTH_SEQUENCE - mut_0_idx - 1 
    f.write('r{}\t{}\t{}\t1\t255\t{}M1X{}M\t=\t*\t{}\t{}\t*\n'.format(str(r).zfill(len(str(NUM_READS))), flag, construct, pre_mut_bases_count, post_mut_bases_count, LENGTH_SEQUENCE, seq))

def print_bv_txt_header(f, construct, sequence):
    f.write('@ref\t{}\t{}\tDMS\n'.format(construct,sequence))
    f.write('@coordinates:\t0,{}:{}\n'.format(LENGTH_SEQUENCE, LENGTH_SEQUENCE))
    f.write('Query_name\tBit_vector\tN_Mutations\n')

def print_bv_txt_line(f, construct, read, r, mut_0_idx, seq):
    f.write(read_name(construct, read, r, mut_0_idx) + '\t' + ''.join(['0']*mut_0_idx) + seq[mut_0_idx] + ''.join(['0']*(LENGTH_SEQUENCE-mut_0_idx-1))+ '\t1\n')

def make_samples():
    try:
        pd.DataFrame({
            'sample': SAMPLE_NAME,
            'user': USER,
            'date': DATE,
            'exp_env': EXP_ENV,
            'temperature_k': TEMPERATURE_K,
            'inc_time_tot_secs': INC_TIME_TOT_SECS,
            'cell_line': CELL_LINE
        }, index=[0]).to_csv(os.path.join(DIRECTORY, 'samples.csv'), index=False)
        return 1
    except Exception as e:
        return e

def make_library():
    try: 
        df = pd.DataFrame({'construct': 2*CONSTRUCTS, 
                        'barcode_start': 2*[BARCODE_START for i in range(NUM_CONSTRUCTS)],
                        'barcode_end': 2*[BARCODE_END for i in range(NUM_CONSTRUCTS)],  
                        'barcode': 2*BARCODES,
                        'section': [SECTION_1_NAME for i in range(NUM_CONSTRUCTS)] + [SECTION_2_NAME for i in range(NUM_CONSTRUCTS)],
                        'section_start': [SECTION_1_START for i in range(NUM_CONSTRUCTS)] + [SECTION_2_START for i in range(NUM_CONSTRUCTS)],
                        'section_end': [SECTION_1_END for i in range(NUM_CONSTRUCTS)] + [SECTION_2_END for i in range(NUM_CONSTRUCTS)],
                        }).to_csv(os.path.join(DIRECTORY, 'library.csv'), index=False)
        return df
    except:
        return 0

def make_fasta():
    with open(os.path.join(DIRECTORY, SAMPLE_NAME + '.fasta'), 'w') as f:
        for construct in CONSTRUCTS:
            print_fasta_line(f, construct, build_sequence_from_barcode(construct, SEQUENCE, BARCODE_START, BARCODE_END))
    return 1

def make_fastq(read):
    barcodes, sequence, barcode_start, barcode_end = get_sequence_from_read(read)
    with open(os.path.join(DIRECTORY, '{}_R{}.fastq'.format(SAMPLE_NAME, read)), 'w') as f:
        for barcode, construct in zip(barcodes, CONSTRUCTS):
            this_sequence = build_sequence_from_barcode(barcode, sequence, barcode_start, barcode_end)
            for r in range(NUM_READS):
                mut_idx, this_mutated_sequence = add_mutation(r, this_sequence, give_idx=True)
                print_fastq_line(f, read_name(construct, read, r, mut_idx), this_mutated_sequence, 'F'*LENGTH_SEQUENCE)
    return 1


def make_demultiplexed_fastq(read):
    directory = os.path.join(DIRECTORY,'demultiplexed_fastq')
    if not os.path.exists(directory):
        os.makedirs(directory)
    barcodes, sequence, barcode_start, barcode_end = get_sequence_from_read(read)
    for barcode, construct in zip(barcodes, CONSTRUCTS):
        with open(os.path.join(directory, '{}_R{}.fastq'.format(construct, read)), 'w') as f:
            this_sequence = build_sequence_from_barcode(barcode, sequence, barcode_start, barcode_end)
            for r in range(NUM_READS):
                mut_idx, this_mutated_sequence = add_mutation(r, this_sequence, give_idx=True)
                print_fastq_line(f, read_name(construct, read, r, mut_idx), this_mutated_sequence, 'F'*LENGTH_SEQUENCE)
    return 1

def make_sam(construct):
    directory = os.path.join(DIRECTORY, 'sam')
    if not os.path.exists(directory):
        os.makedirs(directory)
    with open(os.path.join(directory, '{}.sam'.format(construct)), 'w') as f:
        print_sam_header(f, construct)
        for r in range(NUM_READS):
            for read in [1,2]:
                barcodes, sequence, barcode_start, barcode_end = get_sequence_from_read(read)
                this_sequence = build_sequence_from_barcode(barcodes[CONSTRUCTS.index(construct)], sequence, barcode_start, barcode_end)
                mut_idx, this_mutated_sequence = add_mutation(r, this_sequence, give_idx=True)
                print_sam_line(f, construct, read, r, mut_idx, this_mutated_sequence)
    return 1
    

def make_bv_orc_from_txt(txt_path, output_path):
    with open(txt_path, 'r') as f:
        line = f.readline().split('\t')
        ref_seq = line[2]
    df = pd.read_csv(txt_path, sep='\t', skiprows=2).set_index('Query_name')
    df.drop('N_Mutations', axis=1, inplace=True)
    df = df.join(df['Bit_vector'].str.split('', expand=True)).drop(['Bit_vector',0,201], axis=1)
    df.columns = [str(c)+str(i) for i,c in enumerate(ref_seq)]
    df.to_orc(output_path, index=False)     

def make_bitvector(construct):
    output_folder = os.path.join(DIRECTORY,'bitvector')
    temp_folder = os.path.join(DIRECTORY, 'bitvector', 'temp')
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(temp_folder):
        os.makedirs(temp_folder)
    
    with open(os.path.join(temp_folder, construct)+'.txt', 'w') as f:
        sequence = build_sequence_from_barcode(BARCODES[CONSTRUCTS.index(construct)], SEQUENCE, BARCODE_START, BARCODE_END)
        print_bv_txt_header(f, construct, sequence)
        for r in range(NUM_READS):
            for read in [1,2]:
                mut_idx, this_mutated_sequence = add_mutation(r, sequence, give_idx=True)
                print_bv_txt_line(f, construct, read, r, mut_idx, this_mutated_sequence)
    make_bv_orc_from_txt(os.path.join(temp_folder, construct)+'.txt', os.path.join(output_folder, construct)+'.orc')
    return 1


def make_files():
    os.makedirs(DIRECTORY)
    make_samples()
    make_library()
    make_fasta()
    for read in [1,2]:
        make_fastq(read)
        make_demultiplexed_fastq(read)
    for construct in CONSTRUCTS:
        make_sam(construct)
        make_bitvector(construct)
