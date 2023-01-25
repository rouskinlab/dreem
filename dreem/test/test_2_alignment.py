
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
import dreem
import dreem.util.util as util
import pandas as pd
from dreem.test import files_generator
from dreem.test.files_generator import test_files_dir, input_dir, prediction_dir, output_dir
import pytest
import numpy as np

# Function creating lists of mutations, insertions or deletions
def get_change(n_change, bc_pos, seq_length):

    n_change_before_bc = np.round( n_change*bc_pos[0] / (seq_length - (bc_pos[1]-bc_pos[0])) ).astype(np.int64)

    change_before = np.round(np.linspace(2, bc_pos[0]-2, n_change_before_bc)).astype(np.int64)
    change_after = np.round(np.linspace(bc_pos[1]+1, seq_length-2, n_change-n_change_before_bc)).astype(np.int64)

    assert (np.diff(change_before) > 1).all()
    assert (np.diff(change_after) > 1).all()

    return [list(change_before)+list(change_after)]


# Code for test set 1
sample_name_1 = 'test_set_1'
module = 'alignment'
number_of_constructs = 1
number_of_reads = [23]
mutations = [ [[]]*20+[[2, 26, 42]]+[[5, 8, 25, 35, 47]]+[[2, 8, 22, 25, 28, 31, 35, 41, 45, 48]] ]
length = 50
reads = [[files_generator.create_sequence(length)]*number_of_reads[k] for k in range(number_of_constructs)]
insertions = [[[]]*number_of_reads[k] for k in range(number_of_constructs)]
deletions = [[[]]*number_of_reads[k] for k in range(number_of_constructs)]
no_info = [[[]]*number_of_reads[k] for k in range(number_of_constructs)]

constructs = ['construct_{}'.format(i) for i in range(number_of_constructs)]
barcode_start = 10
barcodes = files_generator.generate_barcodes(10, number_of_constructs, 3)
sections_start = [[0]]*number_of_constructs
sections_end = [[5]]*number_of_constructs
sections = [['{}-{}'.format(ss+1, se) for ss,se in zip(sections_start[n], sections_end[n])] for n in range(number_of_constructs)]

sample_profile_1 = files_generator.make_sample_profile(constructs, reads, number_of_reads, mutations, insertions, deletions, no_info, sections=sections, section_start=sections_start, section_end=sections_end, barcodes=barcodes, barcode_start=barcode_start)



# Code for test set 2
sample_name_2 = 'test_set_2'
seq_ls = [50, 150, 600]*3+[50]
number_of_constructs = 10
n_reads = [23, 10, 10, 10, 23, 23, 23, 100, 1, 12]
barcode_start = 10
barcode_len = 8
bc_pos = (barcode_start, barcode_start+barcode_len)

# Generic mutation lists
default_mut = lambda i : [[]]*20 + get_change(3, bc_pos, seq_ls[i]) + get_change(5, bc_pos, seq_ls[i]) + get_change(10, bc_pos, seq_ls[i])
default_del_insert_small = lambda i : [[]]*4 + get_change(1, bc_pos, seq_ls[i])*3 + get_change(2, bc_pos, seq_ls[i])*2 + get_change(3, bc_pos, seq_ls[i])
default_del_insert_large = lambda i : [[]]*20 + get_change(1, bc_pos, seq_ls[i]) + get_change(2, bc_pos, seq_ls[i]) + get_change(3, bc_pos, seq_ls[i])

mutations = [
    default_mut(0),
    [[]]*n_reads[1] ,
    [[]]*n_reads[2] ,
    [[]]*n_reads[3] ,
    default_mut(4),
    default_mut(5),
    default_mut(6),
    [[]]*95 + get_change(5, bc_pos, seq_ls[7])*5,
    [[]]*n_reads[8] ,
    [[]]*n_reads[9] ,
]

deletions = [
    [[]]*n_reads[0],
    default_del_insert_small(1),
    [[]]*n_reads[2] ,
    default_del_insert_small(3),
    default_del_insert_large(4),
    [[]]*n_reads[5],
    default_del_insert_large(6),
    [[]]*n_reads[7],
    [[]]*n_reads[8] ,
    [[]]*n_reads[9] ,
]

insertions = [
    [[]]*n_reads[0],
    [[]]*n_reads[1],
    default_del_insert_small(2),
    default_del_insert_small(3),
    [[]]*n_reads[4] ,
    default_del_insert_large(5),
    default_del_insert_large(6),
    [[]]*n_reads[7],
    [[]]*n_reads[8] ,
    [[]]*n_reads[9] ,
]

no_info = [
    [[]]*n_reads[0],
    [[]]*n_reads[1],
    default_del_insert_small(2),
    default_del_insert_small(3),
    [[]]*n_reads[4] ,
    default_del_insert_large(5),
    [[]]*n_reads[6] ,
    [[]]*n_reads[7],
    [[]]*n_reads[8] ,
    [[]]*n_reads[9] ,
]

for i in range(number_of_constructs):
    for j in range(n_reads[i]):
        for m in mutations[i][j]:
            if m in deletions[i][j]:
                deletions[i][j].remove(m)
            if m in insertions[i][j]:
                insertions[i][j].remove(m)
            if m in no_info[i][j]:
                no_info[i][j].remove(m)
        for d in deletions[i][j]:
            if d in insertions[i][j]:
                insertions[i][j].remove(d)
            if d in no_info[i][j]:
                no_info[i][j].remove(d)

reads = [[files_generator.create_sequence(seq_ls[k])]*n_reads[k] for k in range(number_of_constructs)]
constructs = ['construct_{}'.format(i) for i in range(number_of_constructs)]
barcodes = files_generator.generate_barcodes(barcode_len, number_of_constructs, 3)
sections_start = [[0]]*number_of_constructs
sections_end = [[5]]*number_of_constructs
sections = [['{}-{}'.format(ss+1, se) for ss,se in zip(sections_start[n], sections_end[n])] for n in range(number_of_constructs)]

sample_profile_2 = files_generator.make_sample_profile(constructs, reads, n_reads, mutations, insertions, deletions, no_info, sections=sections, section_start=sections_start, section_end=sections_end, barcodes=barcodes, barcode_start=barcode_start)


sample_profiles = [sample_profile_1, sample_profile_2]
sample_names = [sample_name_1, sample_name_2]

module_input = os.path.join(input_dir, module)
module_expected = os.path.join(prediction_dir, module)
module_output =  os.path.join(output_dir, module)

inputs = ['fastq','fasta']
outputs = ['bam']

# ### Create test files for `test set 1`
def test_make_files():
    os.system('rm -rf {}'.format(os.path.join(test_files_dir, 'output', module)))
    os.system('rm -rf {}'.format(os.path.join(test_files_dir, 'input', module)))
    os.system('rm -rf {}'.format(os.path.join(test_files_dir, 'expected_output', module)))
    if not os.path.exists(os.path.join(test_files_dir, 'input', module)):
        os.makedirs(os.path.join(test_files_dir, 'input', module))
    
    for sample_profile, sample_name in zip(sample_profiles, sample_names):
        files_generator.generate_files(sample_profile, module, inputs, outputs, test_files_dir, sample_name)
        files_generator.assert_files_exist(sample_profile, module, inputs, input_dir, sample_name)
        files_generator.assert_files_exist(sample_profile, module, outputs, prediction_dir, sample_name)

#@pytest.mark.skip(reason="Dependencies not implemented yet")
def test_run():
    for sample in os.listdir(module_input):
        dreem.alignment.run(
            out_dir = os.path.dirname(output_dir),
            fastq = '{}/{}_R1.fastq'.format(os.path.join(module_input,sample),sample),
            fastq2 = '{}/{}_R2.fastq'.format(os.path.join(module_input,sample),sample),
            fasta = '{}/reference.fasta'.format(os.path.join(module_input,sample)),
            demultiplexing = False,
            nextseq_trim=False
        )
        
#def test_copy_prediction_as_results():
#    files_generator.copy_prediction_as_results(module_expected, os.path.join(module_output,'output'))

#@pytest.mark.skip(reason="Dependencies not implemented yet")
def test_output_exists(): 
    for sample_profile, sample_name in zip(sample_profiles, sample_names):       
        files_generator.assert_files_exist(sample_profile, module, outputs, output_dir, sample_name)

#@pytest.mark.skip(reason="Dependencies not implemented yet")
def test_all_files_are_equal():
    for sample_profile, sample_name in zip(sample_profiles, sample_names):
        files_generator.assert_files_exist(sample_profile, module, outputs, output_dir, sample_name)
    for sample in os.listdir(module_input):
        for sam in os.listdir(os.path.join(module_expected,sample)):
            if not sam.endswith('.sam'):
                continue
            files_generator.bam_to_sam(os.path.join(module_output,sample,sam.replace('.sam','.bam')), os.path.join(module_output,sample,sam))
            p, o = util.sam_to_df(os.path.join(module_expected,sample,sam)), util.sam_to_df(os.path.join(module_output,sample,sam))
            both = pd.concat([p,o], ignore_index=True).reset_index(drop=True)
            for (r, f), g in both.groupby(['QNAME','FLAG']):
                if len(g['SEQ'].iloc[0]) == 50 and g['CIGAR'].iloc[0].count('X') > 3 or \
                   len(g['SEQ'].iloc[0]) == 50 and g['CIGAR'].iloc[0].count('D') > 2 or \
                    len(g['SEQ'].iloc[0]) == 150 and g['CIGAR'].iloc[0].count('X') > 5:
                    assert len(g) == 1, 'SAM: {}: Read {} with flag {} and CIGAR {} is here and shouldn\'t be'.format(sam, r,f,g['CIGAR'].iloc[0])
                    continue
                
                assert len(g) == 2, 'SAM: {}: Read {} with flag {} is missing'.format(sam, r,f)
                for col in g.columns:
                    if col not in ['RNEXT', 'PNEXT', 'MAPQ','TLEN','CIGAR']:
                        val_pred = g[col].iloc[0]
                        val_out = g[col].iloc[1]
                        assert val_out == val_pred, 'SAM: {}: Read {} with flag {} has different {} values: \n{} \n{}'.format(sam,r,f,col, val_pred, val_out)

if __name__ == '__main__':
    test_make_files()
    test_run()
    test_output_exists()
    test_all_files_are_equal()
