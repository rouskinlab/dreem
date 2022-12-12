import dreem, os
import dreem.util as util
import pandas as pd
from dreem.test import files_generator
from dreem.test.files_generator import test_files_dir, input_dir, prediction_dir, output_dir
import pytest
import dreem.pipeline


sample_name = 'test_set_1'
module = 'pipeline'

number_of_constructs = 2
number_of_reads = [10]*2
mutations = [ [[]]+[[25]]+[[35]]+[[]]*4+[[37]]+[[32]]+[[33,36]] for n in number_of_reads ]
insertions = [ [[]]*3+[[11]]+[[10, 21]]+[[]]*2+[[15]]+[[]]*2 for n in number_of_reads ]
deletions = [ [[]]*5+[[2]]+[[4, 6]]+[[]]+[[3]]+[[]] for n in number_of_reads ]

length = [50, 150]
sequences = [[files_generator.create_sequence(length[k])]*number_of_reads[k] for k in range(number_of_constructs)]
constructs = ['construct_{}'.format(i) for i in range(number_of_constructs)]
barcode_start = 30
barcodes = files_generator.generate_barcodes(10, number_of_constructs, 3)
sections_start = [[0, 25],[0, 25, 50, 75]]
sections_end = [[25, 49],[25, 50, 75, 99]]
sections = [['{}_{}'.format(ss, se) for ss,se in zip(sections_start[n], sections_end[n])] for n in range(number_of_constructs)]


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
                    'path_bv': os.path.join(test_files_dir, 'input', module, 'r{}_sl{}_ub{}_sb{}.orc'.format(r, ac, uc, sc)),
                    'path_json': os.path.join(test_files_dir, 'output', module, 'r{}_sl{}_ub{}_sb{}'.format(r, ac, uc, sc))
                }

if mode == 'light':
    samples_clustering = {list(samples_clustering.keys())[15]:samples_clustering[list(samples_clustering.keys())[15]]}

inputs = ['fastq','fasta','samples','library', 'clustering']
outputs = ['output']

sample_profile = files_generator.make_sample_profile(constructs, sequences, number_of_reads, mutations, insertions, deletions, sections=sections, section_start=sections_start, section_end=sections_end, barcodes=barcodes, barcode_start=barcode_start)

module_input = os.path.join(input_dir, module)
module_predicted = os.path.join(prediction_dir, module)
module_output =  test_files_dir

# ### Create test files for `test set 1`
#@pytest.mark.skip(reason="Dependencies not implemented yet")
def test_make_files():
    if not os.path.exists(os.path.join(test_files_dir, 'input', module)):
        os.makedirs(os.path.join(test_files_dir, 'input', module))
    files_generator.generate_files(sample_profile, module, [i for i in inputs if i != 'clustering'], outputs, test_files_dir, sample_name)
    files_generator.generate_files(samples_clustering, module, ['clustering'], [], test_files_dir, sample_name)
    files_generator.assert_files_exist(sample_profile, module, [i for i in inputs if i != 'clustering'], input_dir, sample_name)
    files_generator.assert_files_exist(samples_clustering, module, ['clustering'], input_dir, sample_name)
    files_generator.assert_files_exist(sample_profile, module, outputs, prediction_dir, sample_name)


@pytest.mark.skip(reason="Dependencies not implemented yet")
def test_run():
    for sample in os.listdir(module_input):
        
        dreem.pipeline.run(            
            fastq = '{}/{}_R1.fastq'.format(os.path.join(module_input,sample),sample),\
            fastq2 = '{}/{}_R2.fastq'.format(os.path.join(module_input,sample),sample),\
            fasta = '{}/reference.fasta'.format(os.path.join(module_input,sample)),\
            library = '{}/library.csv'.format(os.path.join(module_input,sample)),\
            samples = '{}/samples.csv'.format(os.path.join(module_input,sample)),\
            
            )
        

def test_copy_prediction_as_results():
    files_generator.copy_prediction_as_results(module_predicted, os.path.join(module_output,'output'))
        
def test_output_exists():        
    files_generator.assert_files_exist(sample_profile, module, outputs, output_dir, sample_name)

@pytest.mark.skip(reason="Dependencies not implemented yet")
def test_all_files_are_equal():
    assert 1 == 0, 'Not implemented yet'

