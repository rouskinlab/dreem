
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
import dreem, os
import dreem.util as util
import pandas as pd
from dreem.test import files_generator
from dreem.test.files_generator import test_files_dir, input_dir, prediction_dir, output_dir
import pytest
import dreem.pipeline
import shutil
from dreem.util.files_sanity import compare_fields
import json
from dreem.test.sample_profile import sample_profile, constructs, sections

sample_name = 'test_set_1'
module = 'pipeline'


inputs = ['fastq','fasta','samples','library']
outputs = ['output']

module_input = os.path.join(input_dir, module)
module_expected = os.path.join(prediction_dir, module)
module_output =  os.path.join(output_dir, module)

# ### Create test files for `test set 1`
def test_make_files():
    for folder in ['input', 'output', 'expected_output']:
        os.system('rm -rf {}'.format(os.path.join(test_files_dir, folder, module)))
        os.makedirs(os.path.join(test_files_dir, folder, module), exist_ok=True)
    files_generator.generate_files(sample_profile, module, inputs, outputs, test_files_dir, sample_name, rnastructure_config={'path':'/Users/ymdt/src/RNAstructure/exe', 'temperature':False, 'fold_args':'', 'dms':False, 'dms_min_unpaired_value':0.01, 'dms_max_paired_value':0.06, 'partition':False, 'probability':False })
    files_generator.assert_files_exist(sample_profile, module, inputs, input_dir, sample_name)
    files_generator.assert_files_exist(sample_profile, module, outputs, prediction_dir, sample_name)
    global expected
    expected = json.load(open(os.path.join(module_expected, sample_name+'.json'), 'r'))

def test_run():
    for sample in os.listdir(module_input):
        
        dreem.pipeline.run(            
            fastq1 = '{}/{}_R1.fastq'.format(os.path.join(module_input,sample),sample),\
            fastq2 = '{}/{}_R2.fastq'.format(os.path.join(module_input,sample),sample),\
            fasta = '{}/reference.fasta'.format(os.path.join(module_input,sample)),\
            library = '{}/library.csv'.format(os.path.join(module_input,sample)),\
            samples = '{}/samples.csv'.format(os.path.join(module_input,sample)),\
            out_dir = module_output,
            rnastructure_path='/Users/ymdt/src/RNAstructure/exe'
            )
    global output
    output = json.load(open(os.path.join(module_output, 'output','aggregation',sample_name+'.json'), 'r'))
        
def test_output_exists():        
    files_generator.assert_files_exist(sample_profile, module, outputs, output_dir, sample_name)

def test_output_exists():        
    files_generator.assert_files_exist(sample_profile, module, outputs, output_dir, sample_name)

@pytest.mark.parametrize('attr',['date', 'sample','user','temperature_k','inc_time_tot_secs','exp_env','DMS_conc_mM'])
def test_sample_attributes(attr):
    compare_fields(expected, output, [attr])

@pytest.mark.parametrize('construct', constructs)
@pytest.mark.parametrize('attr',['sequence', 'barcode', 'barcode_start', 'some_random_attribute'])
def test_library_attributes(construct,attr):
    compare_fields(expected, output, [construct, attr])

@pytest.mark.parametrize('construct', constructs)
@pytest.mark.parametrize('attr',['num_reads', 'num_aligned','skips_short_reads','skips_too_many_muts','skips_low_mapq'])
def test_alignment_attributes(construct,attr):
    compare_fields(expected, output, [construct, attr])

@pytest.mark.parametrize('construct,section', [(c,s) for c in constructs for s in sections[constructs.index(c)]])
@pytest.mark.parametrize('attr',['section_start','section_end','sequence','structure','deltaG'])
def test_sections_attributes(construct,section,attr):
    compare_fields(expected, output, [construct, section, attr])

@pytest.mark.parametrize('construct', constructs)
@pytest.mark.parametrize('attr', ['worst_cov_bases','mut_bases', 'del_bases', 'ins_bases', 'cov_bases', 'info_bases', 'mut_rates', 'mod_bases_A', 'mod_bases_C', 'mod_bases_G', 'mod_bases_T','poisson_high','poisson_low'] ) #'num_of_mutations'
def test_mp_pop_avg(construct,attr):
    for section in sections[constructs.index(construct)]:
        compare_fields(expected, output, [construct, section, 'pop_avg', attr])

@pytest.mark.parametrize('construct', constructs)
def test_section_idx(construct):
    for section in sections[constructs.index(construct)]:
        ss, se = output[construct][section]['section_start'], output[construct][section]['section_end']
        assert len(output[construct][section]['sequence']) == se-ss+1, 'section length is not correct: {} != {}-{}+1'.format(len(output[construct][section]['sequence']), se, ss)
        assert output[construct]['sequence'].index(output[construct][section]['sequence']) == ss-1, 'section start is not correct'
        assert output[construct]['sequence'][ss-1:se] == output[construct][section]['sequence'], 'section sequence is not correct'
        
if __name__ == '__main__':
    # remove test files
    os.system('rm -rf {}'.format(os.path.join(test_files_dir, 'output', module)))
    test_make_files()
    test_run()
