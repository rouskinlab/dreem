#!/usr/bin/env python
# coding: utf-8

# ### Imports
 
import sys, os

try:
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../..')))
except:
    __file__ = os.path.join(os.getcwd(),'make_aggregate_files.ipynb')
    sys.path.append(os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../..')))
    
import numpy as np
import dreem 
import dreem.util as util
import pandas as pd
import os


# ### Create test files for ```test_set_1.py```
# 

sample_name = 'test_set_1'
number_of_constructs = 2
number_of_reads = [10]*number_of_constructs
mutations = [[[25]]*4+[[50,75]]*(n-4) for n in number_of_reads]
length = 100
reads = [[util.create_sequence(length)]*number_of_reads[k] for k in range(number_of_constructs)]
insertions = [[[3]]*n for n in number_of_reads]
deletions = [[[]]*n for n in number_of_reads]
constructs = ['construct_{}'.format(i) for i in range(number_of_constructs)]
barcode_start = 10
barcodes = util.generate_barcodes(8, number_of_constructs, 3)
sections_start = [[0,0, 25, 50, 75]]*number_of_constructs
sections_end = [[100,25, 50, 75, 99]]*number_of_constructs
sections = [['{}_{}'.format(ss, se) for ss,se in zip(sections_start[n], sections_end[n])] for n in range(number_of_constructs)]

sample_profile = util.make_sample_profile(constructs, reads, number_of_reads, mutations, insertions, deletions, sections=sections, section_start=sections_start, section_end=sections_end, barcodes=barcodes, barcode_start=barcode_start)
test_files_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),  '../..', 'test', 'test_files'))

inputs = ['bitvector','samples_csv','library', 'clustering']
outputs = ['output']
util.generate_files(sample_profile, 'aggregate', inputs, outputs, test_files_dir, sample_name)


