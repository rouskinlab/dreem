
import pandas as pd
import dreem.util as util
import os

path_input = os.path.join(os.getcwd(),'test','test_files','input','aggregate')
path_predicted = os.path.join(os.getcwd(),'test','test_files','predicted_output','aggregate')
path_output = os.path.join(os.getcwd(),'test','test_files')

def test_files_exists():
    for sample in os.listdir(path_input):
        cmd = 'dreem-aggregate '\
            +' —id {} '.format(os.path.join(path_input, sample))\
            +' -o ' + path_output \
            +' -fa ' + os.path.join(path_input, sample, 'reference.fasta')\
            +' --clusters ' + os.path.join(path_input, sample, 'clustering.csv')
        util.run_cmd(cmd)

        assert os.path.exists(os.path.join(path_output,'output','aggregate',sample)), 'Output folder for sample {} doesn\'t exist'.format(os.path.join(path_output,'output','aggregate',sample))
        for bv_file in os.listdir(os.path.join(path_predicted,sample)):
            assert os.path.isfile(os.path.join(path_output,'output','aggregate',sample,bv_file)), 'File {} is missing'.format(bv_file)


