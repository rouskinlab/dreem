
import os
import pandas as pd
import dreem.util as util
from dreem import vectoring
from dreem.test import files_generator
from dreem.test.files_generator import test_files_dir, input_dir, prediction_dir, output_dir
import pytest

module = 'vectoring'
sample_name = 'test_set_1'
number_of_constructs = 2
number_of_reads = [10]*number_of_constructs
mutations = [[[25]]*4+[[50,75]]*(n-4) for n in number_of_reads]
length = 100
reads = [[files_generator.create_sequence(length)]*number_of_reads[k] for k in range(number_of_constructs)]
insertions = [[[]]*n for n in number_of_reads]
deletions = [[[]]*n for n in number_of_reads]
constructs = ['construct_{}'.format(i) for i in range(number_of_constructs)]
barcode_start = 10
barcodes = files_generator.generate_barcodes(8, number_of_constructs, 3)
sections_start = [[0, 25, 50, 75]]*number_of_constructs
sections_end = [[25, 50, 75, 99]]*number_of_constructs
sections = [['{}_{}'.format(ss, se) for ss,se in zip(sections_start[n], sections_end[n])] for n in range(number_of_constructs)]

sample_profile = files_generator.make_sample_profile(constructs, reads, number_of_reads, mutations, insertions, deletions, sections=sections, section_start=sections_start, section_end=sections_end, barcodes=barcodes, barcode_start=barcode_start)

module_input = os.path.join(input_dir, module)
module_predicted = os.path.join(prediction_dir, module)
module_output =  os.path.join(output_dir, module)

inputs = ['sam','fasta','library']
outputs = ['bitvector']

def test_make_files():
    if not os.path.exists(os.path.join(test_files_dir, 'input', module)):
        os.makedirs(os.path.join(test_files_dir, 'input', module))
    files_generator.generate_files(sample_profile, module, inputs, outputs, test_files_dir, sample_name)
    files_generator.assert_files_exist(sample_profile, module, inputs, input_dir, sample_name)
    files_generator.assert_files_exist(sample_profile, module, outputs, prediction_dir, sample_name)

def test_run():
    for sample in os.listdir(module_input):
        
        vectoring.run(
            input_dir =os.path.join(module_input, sample),
            out_dir = module_output,
            fasta = os.path.join(module_input, sample, 'reference.fasta'),
            library = os.path.join(module_input, sample, 'library.csv'),
            )

def test_output_exists():        
    files_generator.assert_files_exist(sample_profile, module, outputs, output_dir, sample_name)

@pytest.mark.skip(reason="Dependencies not implemented yet")
def test_files_are_equal():
    files_generator.assert_files_exist(sample_profile, module, outputs, output_dir, sample_name)
    for sample in os.listdir(module_input):
        for construct in os.listdir(os.path.join(module_predicted,sample)):
            for section in os.listdir(os.path.join(module_predicted,sample,construct)):
                p = pd.read_orc(os.path.abspath(os.path.join(module_predicted,sample,construct,section)))
                o = pd.read_orc(os.path.abspath(os.path.join(module_output,'output','vectoring',sample,construct,section)))
                for pp, oo in zip(p.columns, o.columns):
                    assert pp == oo, 'Columns are not the same in predicted col {} and result col {} for sample {} construct {} section {}'.format(pp, oo, sample, construct, section)
                    assert (p[pp] == o[pp]).all(), 'Values are not the same in predicted {} and result {} for sample {} construct {} section {}'.format(pp, oo, sample, construct, section)
