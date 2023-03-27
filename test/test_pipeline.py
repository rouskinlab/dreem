import os, sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import os
import pytest
from dreem.main import run
from dreem.util.files_sanity import compare_fields
import json
import shutil

sample = 'my_test_sample'
test_files = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'test_files')
expected_file = json.load(open(os.path.join(test_files, '{}.json'.format(sample)), 'r'))
top_dir = test_files.replace('test_files', 'test_output')
temp_dir = test_files.replace('test_files', 'test_temp')
output_file_path = os.path.join(top_dir, '{}.json'.format(sample))

FIELDS_FROM_ALIGNMENT_REPORT = ['num_reads', 'skips_low_mapq', 'skips_short_reads', 'skips_too_many_muts']

def output_file():
    assert os.path.exists(output_file_path)
    return json.load(open(output_file_path, 'r'))


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


def test_run():
    if os.path.exists(output_file_path):
        os.remove(output_file_path)
    run(
        out_dir=top_dir,
        temp_dir=temp_dir,
        fastq1=('{}/{}_R1.fastq'.format(test_files, sample),),
        fastq2=('{}/{}_R2.fastq'.format(test_files, sample),),
        fasta='{}/{}.fasta'.format(test_files, sample),
        library='{}/library.csv'.format(test_files, sample),
        samples='{}/samples.csv'.format(test_files, sample),
        rerun=True,
        flat=True,
        rnastructure_path = '/Users/ymdt/src/rnastructure/exe'
    )


def test_output_exists():
    assert os.path.exists(os.path.join(os.getcwd(), 'test_output', 'my_test_sample.json'))


@pytest.mark.parametrize('attr', [k for k in expected_file.keys() if type(expected_file[k]) != dict])
def test_sample_attributes(attr):
    if not type(expected_file[attr]) == dict:
        compare_fields(expected_file, output_file(), [attr])


@pytest.mark.parametrize('reference,attr', [(c, a) for c in list(references_from_json(expected_file)) for a in
                                            attribute_from_reference(expected_file, c)])
def test_library_attributes(reference, attr):
    if attr in FIELDS_FROM_ALIGNMENT_REPORT:
        pytest.skip('skipped')
    else:
        compare_fields(expected_file, output_file(), [reference, attr])


@pytest.mark.parametrize('reference,section', [(c, s) for c in list(references_from_json(expected_file)) for s in
                                               section_from_reference(expected_file, c)])
@pytest.mark.parametrize('attr', attribute_from_section(expected_file))
def test_sections_attributes(reference, section, attr):
    if attr in FIELDS_FROM_ALIGNMENT_REPORT:
        pytest.skip('skipped')
    else:
        compare_fields(expected_file, output_file(), [reference, section, attr])


@pytest.mark.parametrize('reference', references_from_json(expected_file))
@pytest.mark.parametrize('attr', attribute_from_pop_avg(expected_file))
def test_mp_pop_avg(reference, attr):
    for section in section_from_reference(expected_file, reference):
        compare_fields(expected_file, output_file(), [reference, section, 'pop_avg', attr])


@pytest.mark.parametrize('reference', references_from_json(expected_file))
def test_section_idx(reference):
    output = output_file()
    for section in section_from_reference(expected_file, reference):
        ss, se = output[reference][section]['section_start'], output[reference][section]['section_end']
        assert 'sequence' in output[reference][section], 'ref/section/sequence is not found'
        assert len(output[reference][section][
                       'sequence']) == se - ss + 1, 'section length is not correct: {} != {}-{}+1'.format(
            len(output[reference][section]['sequence']), se, ss)


@pytest.mark.parametrize('reference,section', [(c, s) for c in list(references_from_json(expected_file)) for s in
                                               section_from_reference(expected_file, c)])
@pytest.mark.parametrize('plot', ['mutation_fraction', 'mutation_fraction_identity', 'base_coverage'])
def test_draw(reference, section, plot):
    assert os.path.isfile(figure := os.path.join(top_dir, 'draw', sample, '{}__{}__{}.html'.format(reference, section,
                                                                                                   plot))), 'file {} does not exist'.format(
        figure)


if __name__ == '__main__':
    # remove test files
    if os.path.exists(output_file_path):
        os.remove(output_file_path)
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    test_run()
