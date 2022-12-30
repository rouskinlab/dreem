

import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))

import dreem.util.util as util
import pandas as pd

from dreem import demultiplexing

from dreem.test import files_generator
from dreem.test.files_generator import test_files_dir, input_dir, prediction_dir, output_dir
from dreem.test.sample_profile import sample_profile, constructs, sections

module = 'demultiplexing'
sample_name = 'test_set_1'

module_input = os.path.join(input_dir, module)
module_expected = os.path.join(prediction_dir, module)
module_output = os.path.join(output_dir, module)

inputs = ['fastq','library', 'fasta']
outputs = ['demultiplexed_fastq']

# ### Create test files for `test set 1`
def test_make_files():
    if not os.path.exists(os.path.join(test_files_dir, 'input', module)):
        os.makedirs(os.path.join(test_files_dir, 'input', module))
    files_generator.generate_files(sample_profile, module, inputs, outputs, test_files_dir, sample_name)
    files_generator.assert_files_exist(sample_profile, module, inputs, input_dir, sample_name)
    files_generator.assert_files_exist(sample_profile, module, outputs, prediction_dir, sample_name)

# ## Test code for convolution algorithm for demultiplexing

def test_run():
    for sample in os.listdir(module_input):
        demultiplexing.run(
            fastq1 = '{}/{}_R1.fastq'.format(os.path.join(module_input,sample),sample),
            fastq2 = '{}/{}_R2.fastq'.format(os.path.join(module_input,sample),sample),
            fasta= '{}/reference.fasta'.format(os.path.join(module_input,sample),sample),
            library = '{}/library.csv'.format(os.path.join(module_input,sample)),
            out_dir = os.path.join(module_output,sample)
            )

def test_output_exists():        
    files_generator.assert_files_exist(sample_profile, module, outputs, output_dir, sample_name)

def test_all_files_are_equal():
    for sample in os.listdir(module_input):
        for pred, out in zip(os.listdir(os.path.join(module_expected,sample)), [f for f in os.listdir(os.path.join(module_output,sample)) if not f.startswith('lost_reads') and f.endswith('.fastq')]):
            assert pred == out, 'The expected output and the output files are not the same'
            expected = util.fastq_to_df(os.path.join(module_expected,sample,pred))
            expected['from'] = 'expected'
            output = util.fastq_to_df(os.path.join(module_output,sample,out))
            output['from'] = 'output'
            both = pd.concat([expected,output],axis=0, ignore_index=True)
            for idx, g in both.groupby('construct'):
                if len(g) < 2:
                    assert g['from'].iloc[0] == 'expected', 'The output file is missing the construct {} for file {}/{}'.format(idx,sample,out)
                    assert g['from'].iloc[0] == 'output', 'The output file didn\'t filter out the construct {} for file {}/{}'.format(idx,sample,out)
            for idx, g in both.groupby('construct'):
                for c in both.columns:
                    if c != 'construct' and c != 'from':
                        assert g[c].unique().shape[0] == 1, 'The output file is not the same as the expected output for sample {} and construct {}'.format(sample,idx)
                        