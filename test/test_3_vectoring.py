
import os
import pandas as pd
import dreem.util as util

path_input = os.path.join(os.getcwd(),'test','test_files','input','vectoring')
path_predicted = os.path.join(os.getcwd(),'test','test_files','predicted_output','vectoring')
path_output = os.path.join(os.getcwd(),'test','test_files')


def test_files_exists():
    for sample in os.listdir(path_input):
        cmd = 'dreem-vectoring '\
            +'—-id {} '.format(os.path.join(path_input, sample))\
            +'--o ' + path_output \
            +'-fa ' + os.path.join(path_input, sample, 'reference.fasta')
        util.run_cmd(cmd)

        assert os.path.exists(os.path.join(path_output,'output','vectoring',sample)), 'Output folder for sample {} doesn\'t exist'.format(os.path.join(path_output,'output','vectoring',sample))
        for bv_file in os.listdir(os.path.join(path_predicted,sample)):
            assert os.path.exists(os.path.join(path_output,'output','vectoring',sample,bv_file)), 'File {} is missing'.format(bv_file)

def test_files_are_equal():
    for sample in os.listdir(path_input):
        for construct in os.listdir(os.path.join(path_predicted,sample)):
            for section in os.listdir(os.path.join(path_predicted,sample,construct)):
                p = pd.read_orc(os.path.abspath(os.path.join(path_predicted,sample,construct,section)))
                o = pd.read_orc(os.path.abspath(os.path.join(path_output,'output','vectoring',sample,construct,section)))
                for pp, oo in zip(p.columns, o.columns):
                    assert pp == oo, 'Columns are not the same in predicted col {} and result col {} for sample {} construct {} section {}'.format(pp, oo, sample, construct, section)
                    assert (p[pp] == o[pp]).all(), 'Values are not the same in predicted {} and result {} for sample {} construct {} section {}'.format(pp, oo, sample, construct, section)
