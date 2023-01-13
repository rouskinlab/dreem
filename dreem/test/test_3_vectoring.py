
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
import pandas as pd
import dreem.util as util
from dreem import vectoring
from dreem.test import files_generator
from dreem.test.files_generator import test_files_dir, input_dir, prediction_dir, output_dir
import pytest
from dreem.test.sample_profile import sample_profile, constructs, sections

module = 'vectoring'
sample_name = 'test_set_1'

module_input = os.path.join(input_dir, module)
module_expected = os.path.join(prediction_dir, module)
module_output =  os.path.join(output_dir, module)

inputs = ['bam','fasta','library']
outputs = ['bitvector']

def test_make_files():
    os.system('rm -rf {}'.format(os.path.join(test_files_dir, 'output', module)))
    os.system('rm -rf {}'.format(os.path.join(test_files_dir, 'input', module)))
    os.system('rm -rf {}'.format(os.path.join(test_files_dir, 'expected_output', module)))

    if not os.path.exists(os.path.join(test_files_dir, 'input', module)):
        os.makedirs(os.path.join(test_files_dir, 'input', module))
    files_generator.generate_files(sample_profile, module, inputs, outputs, test_files_dir, sample_name)
    files_generator.assert_files_exist(sample_profile, module, inputs, input_dir, sample_name)
    files_generator.assert_files_exist(sample_profile, module, outputs, prediction_dir, sample_name)

def test_run():
    for sample in os.listdir(module_input):
        vectoring.run(
            bam_dirs = [os.path.join(module_input, sample, f) for f in os.listdir(os.path.join(module_input, sample)) if f.endswith('.bam')],
            out_dir = os.path.join(module_output),
            fasta = os.path.join(module_input, sample, 'reference.fasta'),
            library = os.path.join(module_input, sample, 'library.csv'),
            )

def test_output_exists():        
    files_generator.assert_files_exist(sample_profile, module, outputs, output_dir, 'vectoring/'+sample_name)

def test_files_are_equal():
    files_generator.assert_files_exist(sample_profile, module, outputs, output_dir, 'vectoring/'+sample_name)
    for sample in os.listdir(module_input):
        for construct in os.listdir(os.path.join(module_expected,sample)):
            for section in os.listdir(os.path.join(module_expected,sample,construct)):
                p = pd.read_orc(os.path.abspath(os.path.join(module_expected,sample,construct,section,'0.orc')))
                o = pd.read_orc(os.path.abspath(os.path.join(module_output,'vectoring',sample,construct,section,'0.orc')))
                for pp, oo in zip(p.columns, o.columns):
                    assert pp == oo, 'Columns are not the same in expected col {} and result col {} for sample {} construct {} section {}'.format(pp, oo, sample, construct, section)
                    assert (p[pp] == o[pp]).all(), 'Values are not the same in expected {} and result {} for sample {} construct {} section {}'.format(pp, oo, sample, construct, section)

if __name__ == '__main__':
    # remove all files
    os.system('rm -rf {}'.format(os.path.join(test_files_dir, 'output', module)))
    os.system('rm -rf {}'.format(os.path.join(test_files_dir, 'input', module)))
    os.system('rm -rf {}'.format(os.path.join(test_files_dir, 'expected_output', module)))
    test_make_files()
    test_run()