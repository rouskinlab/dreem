#!/usr/bin/env python
# coding: utf-8

# ### About this notebook
# 
# Create a test set for demultiplexing.
# 
# 1 sample, 1 construct, 2 groups of mutations, one has 2 mutations and the other has 1, 10 reads, 100nt. 

# ### Imports


import sys, os
try:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../..')))
except:
    __file__ = os.path.join(os.getcwd(),'make_alignment_files.ipynb')
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../..')))

import numpy as np
import dreem 
import dreem.util as util
import pandas as pd
import os

def get_change(n_change, bc_pos, seq_length):

    n_change_before_bc = np.round( n_change*bc_pos[0] / (seq_length - (bc_pos[1]-bc_pos[0])) ).astype(np.int64)

    change_before = np.round(np.linspace(2, bc_pos[0]-2, n_change_before_bc)).astype(np.int64)
    change_after = np.round(np.linspace(bc_pos[1]+1, seq_length-2, n_change-n_change_before_bc)).astype(np.int64)

    assert (np.diff(change_before) > 1).all()
    assert (np.diff(change_after) > 1).all()

    return [list(change_before)+list(change_after)]


# ### Create test files for `test set 1`
# - fasta file
# - pair of fastq files

sample_name = 'test_set_1'
number_of_constructs = 1
number_of_reads = [23]
mutations = [ [[]]*20+[[2, 26, 42]]+[[5, 8, 25, 35, 47]]+[[2, 8, 22, 25, 28, 31, 35, 41, 45, 48]] ]
length = 50
reads = [[util.create_sequence(length)]*number_of_reads[k] for k in range(number_of_constructs)]
insertions = [[[]]*n for n in number_of_reads]
deletions = [[[]]*n for n in number_of_reads]
constructs = ['construct_{}'.format(i) for i in range(number_of_constructs)]
barcode_start = 10
barcodes = util.generate_barcodes(10, number_of_constructs, 3)
sections_start = [[0, 25, 50, 75]]*number_of_constructs
sections_end = [[25, 50, 75, 99]]*number_of_constructs
sections = [['{}_{}'.format(ss, se) for ss,se in zip(sections_start[n], sections_end[n])] for n in range(number_of_constructs)]

sample_profile = util.make_sample_profile(constructs, reads, number_of_reads, mutations, insertions, deletions, sections=sections, section_start=sections_start, section_end=sections_end, barcodes=barcodes, barcode_start=barcode_start)
test_files_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),  '../..', 'test', 'test_files'))

inputs = ['fastq','fasta']
outputs = ['sam']
util.generate_files(sample_profile, 'alignment', inputs, outputs, test_files_dir, sample_name)


# ### Create test files for `test set 2`
# - fasta file
# # - pair of fastq files
sample_name = 'test_set_2'
seq_ls = [50, 150, 600]*3+[50]
number_of_constructs = 10
n_reads = [23, 10, 10, 10, 23, 23, 23, 100, 1, 12]
barcode_start = 10
barcode_len = 8
bc_pos = (barcode_start, barcode_start+barcode_len)

default_mut = lambda i : [[]]*20 + get_change(3, bc_pos, seq_ls[i]) + get_change(5, bc_pos, seq_ls[i]) + get_change(10, bc_pos, seq_ls[i])
default_del_insert_small = lambda i : [[]]*4 + get_change(1, bc_pos, seq_ls[i])*3 + get_change(2, bc_pos, seq_ls[i])*2 + get_change(3, bc_pos, seq_ls[i])
default_del_insert_large = lambda i : [[]]*20 + get_change(1, bc_pos, seq_ls[i]) + get_change(2, bc_pos, seq_ls[i]) + get_change(3, bc_pos, seq_ls[i])

mutations = [
    default_mut(0),
    [[]]*n_reads[1] ,
    [[]]*n_reads[2] ,
    [[]]*n_reads[3] ,
    default_mut(4),
    default_mut(5),
    default_mut(6),
    [[]]*95 + get_change(5, bc_pos, seq_ls[7])*5,
    [[]]*n_reads[8] ,
    [[]]*n_reads[9] ,
]

deletions = [
    [[]]*n_reads[0],
    default_del_insert_small(1),
    [[]]*n_reads[2] ,
    default_del_insert_small(3),
    default_del_insert_large(4),
    [[]]*n_reads[5],
    default_del_insert_large(6),
    [[]]*n_reads[7],
    [[]]*n_reads[8] ,
    [[]]*n_reads[9] ,
]

insertions = [
    [[]]*n_reads[0],
    [[]]*n_reads[1],
    default_del_insert_small(2),
    default_del_insert_small(3),
    [[]]*n_reads[4] ,
    default_del_insert_large(5),
    default_del_insert_large(6),
    [[]]*n_reads[7],
    [[]]*n_reads[8] ,
    [[]]*n_reads[9] ,
]

reads = [[util.create_sequence(seq_ls[k])]*n_reads[k] for k in range(number_of_constructs)]
constructs = ['construct_{}'.format(i) for i in range(number_of_constructs)]
barcodes = util.generate_barcodes(barcode_len, number_of_constructs, 3)
sections_start = [[0]]*number_of_constructs
sections_end = [[5]]*number_of_constructs
sections = [['{}_{}'.format(ss, se) for ss,se in zip(sections_start[n], sections_end[n])] for n in range(number_of_constructs)]

sample_profile = util.make_sample_profile(constructs, reads, n_reads, mutations, insertions, deletions, sections=sections, section_start=sections_start, section_end=sections_end, barcodes=barcodes, barcode_start=barcode_start)
test_files_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),  '../..', 'test', 'test_files'))

inputs = ['fastq','fasta']
outputs = ['sam']
util.generate_files(sample_profile, 'alignment', inputs, outputs, test_files_dir, sample_name)

