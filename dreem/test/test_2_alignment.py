import dreem, os
import dreem.util.util as util
import pandas as pd
from dreem.test import files_generator
from dreem.test.files_generator import test_files_dir, input_dir, prediction_dir, output_dir
import pytest

sample_name = 'test_set_1'
module = 'alignment'
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

inputs = ['fastq','fasta']
outputs = ['sam']

# ### Create test files for `test set 1`
def test_make_files():
    if not os.path.exists(os.path.join(test_files_dir, 'input', module)):
        os.makedirs(os.path.join(test_files_dir, 'input', module))
    files_generator.generate_files(sample_profile, module, inputs, outputs, test_files_dir, sample_name)
    files_generator.assert_files_exist(sample_profile, module, inputs, input_dir, sample_name)
    files_generator.assert_files_exist(sample_profile, module, outputs, prediction_dir, sample_name)

#@pytest.mark.skip(reason="Dependencies not implemented yet")
def test_run():
    for sample in os.listdir(module_input):
        dreem.alignment.run(            
            fastq1 = '{}/{}_R1.fastq'.format(os.path.join(module_input,sample),sample),\
            fastq2 = '{}/{}_R2.fastq'.format(os.path.join(module_input,sample),sample),\
            fasta = '{}/reference.fasta'.format(os.path.join(module_input,sample)),\
            out_dir = module_output,\
            demultiplexed = False
            )
        
#def test_copy_prediction_as_results():
#    files_generator.copy_prediction_as_results(module_predicted, os.path.join(module_output,'output'))

#@pytest.mark.skip(reason="Dependencies not implemented yet")
def test_files_exists():        
    files_generator.assert_files_exist(sample_profile, module, outputs, output_dir, sample_name)

#@pytest.mark.skip(reason="Dependencies not implemented yet")
def test_all_files_are_equal():
    files_generator.assert_files_exist(sample_profile, module, outputs, output_dir, sample_name)
    for sample in os.listdir(module_input):
        for pred, out in zip(os.listdir(os.path.join(module_predicted,sample)), os.listdir(os.path.join(module_output,'output','alignment',sample))):
            if not pred.endswith('.sam'):
                continue
            p, o = util.sam_to_df(os.path.join(module_predicted,sample,pred)), util.sam_to_df(os.path.join(module_output,'output','alignment',sample,out))
            both = pd.concat([p,o], ignore_index=True).reset_index(drop=True)
            for (r, f), g in both.groupby(['QNAME','FLAG']):
                assert len(g) == 2, 'Read {} with flag {} is missing'.format(r,f)
                for col in g.columns:
                    if col not in ['RNEXT', 'PNEXT']:
                        assert g[col].iloc[0] == g[col].iloc[1], 'Read {} with flag {} has different {} values'.format(r,f,col)
 