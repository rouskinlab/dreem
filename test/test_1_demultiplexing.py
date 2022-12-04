
import os
import dreem.util as util
import pandas as pd

path_input = os.path.join(os.getcwd(),'test','test_files','input','demultiplexing')
path_predicted = os.path.join(os.getcwd(),'test','test_files','predicted_output','demultiplexing')
path_output = os.path.join(os.getcwd(),'test','test_files')

def test_written():
    assert 0, 'Test isn\'t written yet'

def test_make_output_folder():
    util.clear_folder(os.path.join(os.getcwd(),'test/output'))

def test_files_exists():
    for sample in os.listdir(path_input):
        cmd = 'dreem-demultiplexing '\
            +'—-fastq1 \'{}/{}_R1.fastq\' '.format(os.path.join(path_input,sample),sample)\
            +'--fastq2 \'{}/{}_R2.fastq\' '.format(os.path.join(path_input,sample),sample)\
            +'—-library \'{}/library.csv\' '.format(os.path.join(path_input,sample))\
            +'-o ' + path_output
        util.run_cmd(cmd)
        
        assert os.path.exists(os.path.join(path_output,'demultiplexing',sample)), 'Output folder for sample {} doesn\'t exist'.format(sample) 
        for demultiplexed_file in os.listdir(os.path.join(path_output,'output','demultiplexing',sample)):
            assert os.path.isfile(os.path.join(path_predicted,sample,demultiplexed_file)), 'File {} is missing'.format(demultiplexed_file)
