
import pandas as pd
import dreem.util as util
import os
from dreem import aggregation
from dreem.test import files_generator
from dreem.test.files_generator import test_files_dir, input_dir, prediction_dir, output_dir
import pytest

sample_name = 'test_set_1'
module = 'aggregate'
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

inputs = ['bitvector','samples', 'library', 'fasta'] #clustering #TODO
outputs = ['output']


def test_make_files():
    if not os.path.exists(os.path.join(test_files_dir, 'input', module)):
        os.makedirs(os.path.join(test_files_dir, 'input', module), exist_ok=True)
    files_generator.generate_files(sample_profile, module, inputs, outputs, test_files_dir, sample_name)
    #files_generator.assert_files_exist(sample_profile, module, inputs, input_dir, sample_name)
    #files_generator.assert_files_exist(sample_profile, module, outputs, prediction_dir, sample_name)
    
def test_run():
    for sample in os.listdir(module_input):
        aggregation.run(
            input_dir =os.path.join(module_input, sample),
            out_dir = module_output,
            samples = os.path.join(module_input, sample_name, 'samples.csv'),
            library= os.path.join(module_input, sample_name, 'library.csv'),
            sample=sample,
            fasta = os.path.join(module_input, sample_name, 'reference.fasta')
            #clusters = os.path.join(module_input, sample, 'clustering.csv') #TODO
        )

@pytest.mark.skip(reason="Dependencies not implemented yet")
def test_output_exists():        
    files_generator.assert_files_exist(sample_profile, module, outputs, output_dir, sample_name)

@pytest.mark.skip(reason="Dependencies not implemented yet")
def test_files_are_equal():
    assert 1==0, 'not implemented'


if __name__ == '__main__':
    test_make_files()
    test_run()

