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

sample_name = 'test_set_1'
module = 'pipeline'

number_of_constructs = 2
number_of_reads = [10]*2
mutations = [ [[]]+[[24]]+[[35]]+[[]]*4+[[37]]+[[32]]+[[33,36]] for n in number_of_reads ]
insertions = [ [[]]*3+[[11]]+[[10, 21]]+[[]]*2+[[15]]+[[]]*2 for n in number_of_reads ]
deletions = [ [[]]*5+[[2]]+[[4, 6]]+[[]]+[[3]]+[[]] for n in number_of_reads ]
no_info = [ [[]]*2+[[2]]+[[4, 6]]+[[]]+[[3]]+[[]]*5 for n in range(number_of_constructs) ] # 0-based

length = [50, 150]
sequences = [[files_generator.create_sequence(length[k])]*number_of_reads[k] for k in range(number_of_constructs)]
constructs = ['construct_{}'.format(i) for i in range(number_of_constructs)]
barcode_start = 30
barcodes = files_generator.generate_barcodes(10, number_of_constructs, 3)
sections_start = [[0, 25],[0, 25, 50, 75]]
sections_end = [[25, 49],[25, 50, 75, 99]]
sections = [['{}-{}'.format(ss+1, se) for ss,se in zip(sections_start[n], sections_end[n])] for n in range(number_of_constructs)] # 0-based

mode = 'light' # only uses the first sample

reads_partition = [
    [5000], [10000], [20000],
    [25000, 25000], [20000,80000], [10000,90000],
    [10000,20000,70000], [33333,33333,33334], [20000,40000,40000] 
]

half_sequence_length = [50, 100, 150]

unpaired_bases = {k:[int(k*l) for l in  [0.2, 0.4, 0.6]] for k in half_sequence_length}

shared_bases = {}
for k, v in unpaired_bases.items():
    shared_bases[k] = []
    for l in v:
        shared_bases[k].append([int(u*l) for u in  [0.2, 0.4, 0.6]])
        
samples_clustering = {}
for r in reads_partition:
    for ac in half_sequence_length:
        for uc in unpaired_bases[ac]:
            for sc in shared_bases[ac][unpaired_bases[ac].index(uc)]:
                samples_clustering['r{}_sl{}_ub{}_sb{}'.format(r, ac, uc, sc)] = {
                    'n_reads': r,
                    'n_AC': ac,
                    'n_unpaired': uc,
                    'n_shared': sc,
                    'path_bv': os.path.join(test_files_dir, 'input', module, sample_name, 'r{}_sl{}_ub{}_sb{}.orc'.format(r, ac, uc, sc)),
                    'path_json': os.path.join(test_files_dir, 'output', module, 'r{}_sl{}_ub{}_sb{}'.format(r, ac, uc, sc))
                }

if mode == 'light':
    samples_clustering = {list(samples_clustering.keys())[15]:samples_clustering[list(samples_clustering.keys())[15]]}

inputs = ['fastq','fasta','samples','library']
outputs = ['output']

sample_profile = files_generator.make_sample_profile(constructs, sequences, number_of_reads, mutations, insertions, deletions, no_info, sections=sections, section_start=sections_start, section_end=sections_end, barcodes=barcodes, barcode_start=barcode_start)

module_input = os.path.join(input_dir, module)
module_expected = os.path.join(prediction_dir, module)
module_output =  os.path.join(output_dir, module)

# ### Create test files for `test set 1`
def test_make_files():
    os.system('rm -rf {}'.format(os.path.join(test_files_dir, 'output', module)))

    if not os.path.exists(os.path.join(test_files_dir, 'input', module)):
        os.makedirs(os.path.join(test_files_dir, 'input', module))
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
@pytest.mark.parametrize('attr', ['num_of_mutations', 'worst_cov_bases','mut_bases', 'del_bases', 'ins_bases', 'cov_bases', 'info_bases', 'mut_rates', 'mod_bases_A', 'mod_bases_C', 'mod_bases_G', 'mod_bases_T','poisson_high','poisson_low'])
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
