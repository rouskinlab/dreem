
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
import os
import pytest
from dreem.main import run
from dreem.util.files_sanity import compare_fields
import json

sample = 'my_test_sample'
test_files = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), 'test_files')
expected_file = json.load(open(os.path.join(test_files, '{}.json'.format(sample)), 'r'))
top_dir = test_files.replace('test_files', 'test_output')
temp_dir = test_files.replace('test_files', 'test_temp')
output_file_path = os.path.join(top_dir, '{}.json'.format(sample))

FIELDS_FROM_ALIGNMENT_REPORT = ['num_reads','skips_low_mapq','skips_short_reads','skips_too_many_muts']

def output_file():
    return json.load(open(output_file_path, 'r'))

def references_from_json(json_file):
    return [k for k in json_file.keys() if type(json_file[k]) == dict]

def attribute_from_reference(json_file):
    reference = references_from_json(json_file)[0]
    return json_file[reference].keys()

def section_from_reference(json_file):
    reference = references_from_json(json_file)[0]
    return [k for k in json_file[reference].keys() if type(json_file[reference][k]) == dict]

def attribute_from_section(json_file):
    reference = references_from_json(json_file)[0]
    section = section_from_reference(json_file)[0]
    return json_file[reference][section].keys()

def attribute_from_pop_avg(json_file):
    reference = references_from_json(json_file)[0]
    section = section_from_reference(json_file)[0]
    return json_file[reference][section]['pop_avg'].keys()

def test_run():        
    run(
        out_dir= top_dir,      
        temp_dir = temp_dir,
        fastq1 = ('{}/{}_R1.fastq'.format(test_files,sample),),
        fastq2 = ('{}/{}_R2.fastq'.format(test_files,sample),),
        fasta = '{}/{}.fasta'.format(test_files,sample),
        library = '{}/library.csv'.format(test_files,sample),
        samples = '{}/samples.csv'.format(test_files,sample),
        rerun=True,
        rnastructure_path = '/Users/ymdt/src/RNAstructure/exe',
        )
        
def test_output_exists():        
    assert os.path.exists(os.path.join(os.getcwd(),'test_output','my_test_sample.json'))

@pytest.mark.parametrize('attr', expected_file.keys())
def test_sample_attributes(attr):
    if not type(expected_file[attr]) == dict:
        compare_fields(expected_file, output_file(), [attr])

@pytest.mark.parametrize('reference', references_from_json(expected_file))
@pytest.mark.parametrize('attr', attribute_from_reference(expected_file))
def test_library_attributes(reference,attr):
    if not type(expected_file[reference][attr]) == dict and not attr in FIELDS_FROM_ALIGNMENT_REPORT:
        compare_fields(expected_file, output_file(), [reference, attr])

@pytest.mark.parametrize('reference,section', [(c,s) for c in list(references_from_json(expected_file)) for s in section_from_reference(expected_file)])
@pytest.mark.parametrize('attr', attribute_from_section(expected_file))
def test_sections_attributes(reference,section,attr):
    if not type(expected_file[reference][section][attr]) == dict and not attr in FIELDS_FROM_ALIGNMENT_REPORT:
        compare_fields(expected_file, output_file(), [reference, section, attr])

@pytest.mark.parametrize('reference', references_from_json(expected_file))
@pytest.mark.parametrize('attr', attribute_from_pop_avg(expected_file) ) 
def test_mp_pop_avg(reference,attr):
    for section in section_from_reference(expected_file):
        compare_fields(expected_file, output_file(), [reference, section, 'pop_avg', attr])

@pytest.mark.parametrize('reference', references_from_json(expected_file))
def test_section_idx(reference):
    output = output_file()
    for section in section_from_reference(expected_file):
        ss, se = output[reference][section]['section_start'], output[reference][section]['section_end']
        assert len(output[reference][section]['sequence']) == se-ss+1, 'section length is not correct: {} != {}-{}+1'.format(len(output[reference][section]['sequence']), se, ss)
        assert output[reference]['sequence'].index(output[reference][section]['sequence']) == ss-1, 'section start is not correct'
        assert output[reference]['sequence'][ss-1:se] == output[reference][section]['sequence'], 'section sequence is not correct'
        
if __name__ == '__main__':
    # remove test files
    os.system('rm {}'.format(output_file_path))
    test_run()
