import os, sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import os
import pytest
from dreem.main import run
from dreem.util.files_sanity import compare_fields
import json
import shutil

dreem_root = os.path.dirname(os.path.dirname(__file__))
test_files = os.path.join(dreem_root, 'test_files')
top_dir = test_files.replace('test_files', 'test_output')
temp_dir = test_files.replace('test_files', 'test_temp')
samples = ['my_python_sample', 'my_cli_sample']

FIELDS_FROM_ALIGNMENT_REPORT = ['num_reads', 'skips_low_mapq', 'skips_short_reads', 'skips_too_many_muts']

def get_expected_file(sample):
    return json.load(open(os.path.join(test_files, sample, '{}.json'.format(sample)), 'r'))

def get_output_file_path(sample):
    return os.path.join(top_dir, '{}.json'.format(sample))

def output_file(sample):
    assert os.path.exists(get_output_file_path(sample))
    return json.load(open(get_output_file_path(sample), 'r'))


def references_from_json(json_file):
    return [k for k in json_file.keys() if type(json_file[k]) == dict]


def attribute_from_reference(json_file, ref):
    return [k for k in json_file[ref].keys() if type(json_file[ref][k]) != dict]


def section_from_reference(json_file, ref):
    return [k for k in json_file[ref].keys() if type(json_file[ref][k]) == dict]


def attribute_from_section(json_file):
    reference = references_from_json(json_file)[0]
    section = section_from_reference(json_file, reference)[0]
    return [k for k in json_file[reference][section].keys() if type(json_file[reference][section][k]) != dict]


def attribute_from_pop_avg(json_file):
    reference = references_from_json(json_file)[0]
    section = section_from_reference(json_file, reference)[0]
    return json_file[reference][section]['pop_avg'].keys()


def test_clear_output():
    if os.path.exists(top_dir):
        shutil.rmtree(top_dir)
    os.makedirs(top_dir)
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    os.makedirs(temp_dir)

@pytest.mark.parametrize('sample', ['my_python_sample'])
def test_run_python(sample):
    run(
        out_dir=top_dir,
        temp_dir=temp_dir,
        fastq1=(os.path.join(test_files, sample, '{}_R1.fastq'.format(sample)),),
        fastq2=(os.path.join(test_files, sample, '{}_R2.fastq'.format(sample)),),
        fasta=os.path.join(test_files, sample, '{}.fasta'.format(sample)),
        library=os.path.join(test_files, sample, 'library.csv'),
        samples=os.path.join(test_files, sample, 'samples.csv'),
        rerun=True,
        flat=True,
        #rnastructure_path='/Users/ymdt/src/RNAstructure/exe/',
    )

@pytest.mark.parametrize('sample', ['my_cli_sample'])
def test_run_cli(sample):
    os.system('dreem --out-dir {} --temp-dir {} --fastq1 {} --fastq2 {} --fasta {} --library {} --samples {} --section roi --section ms2 --section 1-50 --rerun '.format(
        top_dir,
        temp_dir,
        os.path.join(test_files, sample, '{}_R1.fastq'.format(sample)),
        os.path.join(test_files, sample, '{}_R2.fastq'.format(sample)),
        os.path.join(test_files, sample, '{}.fasta'.format(sample)),
        os.path.join(test_files, sample, 'library.csv'),
        os.path.join(test_files, sample, 'samples.csv'),
    ))

@pytest.mark.parametrize('sample', samples)
def test_output_exists(sample):
    assert os.path.exists(os.path.join(os.getcwd(), 'test_output', sample+'.json'))


@pytest.mark.parametrize('sample,attr', [(s,k) for s in samples for k in get_expected_file(s).keys() if type(get_expected_file(s)[k]) != dict])
def test_sample_attributes(sample, attr):
    if not type(get_expected_file(sample)[attr]) == dict:
        compare_fields(get_expected_file(sample), output_file(sample), [attr])


@pytest.mark.parametrize('sample,reference,attr', [(s, c, a) for s in samples \
                                                             for c in list(references_from_json(get_expected_file(s))) 
                                                             for a in attribute_from_reference(get_expected_file(s), c)])
def test_library_attributes(sample, reference, attr):
    if attr in FIELDS_FROM_ALIGNMENT_REPORT:
        pytest.skip('skipped')
    else:
        compare_fields(get_expected_file(sample), output_file(sample), [reference, attr])


@pytest.mark.parametrize('sample,reference,section,attr', [(ss, c, s, a) for ss in samples \
                                                                 for c in list(references_from_json(get_expected_file(ss))) 
                                                                 for s in section_from_reference(get_expected_file(ss), c)
                                                                 for a in attribute_from_section(get_expected_file(ss))])
def test_sections_attributes(sample, reference, section, attr):
    if attr in FIELDS_FROM_ALIGNMENT_REPORT:
        pytest.skip('skipped')
    else:
        compare_fields(get_expected_file(sample), output_file(sample), [reference, section, attr])


@pytest.mark.parametrize('sample,reference,attr', [(s, c, a) for s in samples \
                                                             for c in list(references_from_json(get_expected_file(s)))
                                                             for a in attribute_from_pop_avg(get_expected_file(s))])
def test_mp_pop_avg(sample, reference, attr):
    for section in section_from_reference(get_expected_file(sample), reference):
        compare_fields(get_expected_file(sample), output_file(sample), [reference, section, 'pop_avg', attr])


@pytest.mark.parametrize('sample,reference', [(s,c) for s in samples for c in list(references_from_json(get_expected_file(s)))])
def test_section_idx(sample, reference):
    output = output_file(sample)
    for section in section_from_reference(get_expected_file(sample), reference):
        ss, se = output[reference][section]['section_start'], output[reference][section]['section_end']
        assert 'sequence' in output[reference][section], 'ref/section/sequence is not found'
        assert len(output[reference][section][
                       'sequence']) == se - ss + 1, 'section length is not correct: {} != {}-{}+1'.format(
            len(output[reference][section]['sequence']), se, ss)


@pytest.mark.parametrize('sample,reference,section', [(ss, c, s) for ss in samples \
                                                                 for c in list(references_from_json(get_expected_file(ss))) 
                                                                 for s in section_from_reference(get_expected_file(ss), c)])
@pytest.mark.parametrize('plot', ['mutation_fraction', 'mutation_fraction_identity', 'base_coverage'])
def test_draw(sample, reference, section, plot):
    assert os.path.isfile(figure := os.path.join(top_dir, 'draw', sample, '{}__{}__{}.html'.format(reference, section,
                                                                                                   plot))), 'file {} does not exist'.format(
        figure)


if __name__ == '__main__':
    # remove test files
    test_clear_output()
    test_run_python('my_python_sample')
    test_run_cli('my_cli_sample')