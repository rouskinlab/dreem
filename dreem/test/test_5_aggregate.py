
import pandas as pd
import dreem.util as util
import os
from dreem import aggregation
from dreem.test import files_generator
from dreem.test.files_generator import test_files_dir, input_dir, prediction_dir, output_dir
import pytest
from dreem.util.files_sanity import compare_fields
import json

sample_name = 'test_set_1'
module = 'aggregate'
number_of_constructs = 2
number_of_reads = [10]*2
mutations = [ [[]]+[[25]]+[[35]]+[[]]*4+[[37]]+[[32]]+[[97,99]] for n in number_of_reads ]
insertions = [ [[]]*3+[[11]]+[[10, 21]]+[[]]*2+[[15]]+[[]]*2 for n in number_of_reads ]
deletions = [ [[]]*5+[[2]]+[[4, 6]]+[[]]+[[3]]+[[]] for n in number_of_reads ]

length = [150, 150]
sequences = [[files_generator.create_sequence(length[k])]*number_of_reads[k] for k in range(number_of_constructs)]
constructs = ['construct_{}'.format(i) for i in range(number_of_constructs)]
barcode_start = 30
len_barcode = 10
barcodes = files_generator.generate_barcodes(len_barcode, number_of_constructs, min_barcode_distance=3)
sections_start = [[0, 25],[0, 25, 50, 75]]
sections_start = [[0, 25, 50, 75], [0, 25, 50, 75]]
sections_end = [[25, 49],[25, 50, 75, 99]]
sections_end = [[25, 50, 75, 100],   [25, 50, 75, 100]]
sections = [['{}_{}'.format(ss, se) for ss,se in zip(sections_start[n], sections_end[n])] for n in range(number_of_constructs)]
sample_profile = files_generator.make_sample_profile(constructs, sequences, number_of_reads, mutations, insertions, deletions, sections=sections, section_start=sections_start, section_end=sections_end, barcodes=barcodes, barcode_start=barcode_start)

module_input = os.path.join(input_dir, module)
module_predicted = os.path.join(prediction_dir, module)
module_output =  os.path.join(output_dir, module)

inputs = ['bitvector','samples', 'library', 'fasta', 'sam'] #clustering #TODO
outputs = ['output']


def test_make_files():
    os.system('rm -rf {}'.format(module_output))
    os.system('rm -rf {}'.format(module_input))
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
            fasta = os.path.join(module_input, sample_name, 'reference.fasta'),
            rnastructure_path='/Users/ymdt/src/RNAstructure/exe',
            #clusters = os.path.join(module_input, sample, 'clustering.csv') #TODO
        )

def test_output_exists():        
    files_generator.assert_files_exist(sample_profile, module, outputs, output_dir, sample_name)

if os.path.isfile(os.path.join(module_predicted, sample_name+'.json')):
    predicted = json.load(open(os.path.join(module_predicted, sample_name+'.json'), 'r'))
if os.path.isfile(os.path.join(module_output, sample_name+'.json')):
    output = json.load(open(os.path.join(module_output, sample_name+'.json'), 'r'))

@pytest.mark.parametrize('attr',['date', 'sample','user','temperature_k','inc_time_tot_secs','exp_env','DMS_conc_mM'])
def test_sample_attributes(attr):
    compare_fields(predicted, output, [attr])

@pytest.mark.parametrize('construct', constructs)
@pytest.mark.parametrize('attr',['sequence', 'barcode', 'barcode_start', 'some_random_attribute'])
def test_library_attributes(construct,attr):
    compare_fields(predicted, output, [construct, attr])

@pytest.mark.parametrize('construct', constructs)
@pytest.mark.parametrize('attr',['num_reads', 'num_aligned','worst_cov_bases','skips_short_reads','skips_too_many_muts','skips_low_mapq'])
def test_alignment_attributes(construct,attr):
    compare_fields(predicted, output, [construct, attr])

@pytest.mark.parametrize('construct,section', [(c,s) for c in constructs for s in sections[constructs.index(c)]])
@pytest.mark.parametrize('attr',['section_start','section_end','sequence'])
def test_sections_attributes(construct,section,attr):
    compare_fields(predicted, output, [construct, section, attr])

@pytest.mark.parametrize('construct', constructs)
@pytest.mark.parametrize('attr', ['num_of_mutations', 'mut_bases', 'del_bases', 'ins_bases', 'cov_bases', 'info_bases', 'mut_rates', 'mod_bases_A', 'mod_bases_C', 'mod_bases_G', 'mod_bases_T'])
def test_mp_pop_avg(construct,attr):
    for section in sections[constructs.index(construct)]:
        compare_fields(predicted, output, [construct, section, 'pop_avg', attr])

if __name__ == '__main__':
    test_make_files()
    test_run()
    for attr in ['num_of_mutations', 'mut_bases', 'del_bases', 'ins_bases', 'cov_bases', 'info_bases', 'mut_rates', 'mod_bases_A', 'mod_bases_C', 'mod_bases_G', 'mod_bases_T']:
        test_mp_pop_avg(attr)