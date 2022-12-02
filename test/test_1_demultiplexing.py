
import os
import dreem.util as util
import pandas as pd

path_input = os.path.join(os.getcwd(),'test','test_files')
path_output = os.path.join(os.getcwd(),'test','output')

def test_written():
    assert 1, 'Test isn\t done yet'

def test_make_output_folder():
    util.clear_folder(os.path.join(os.getcwd(),'test/output'))

def test_files_names():
    for sample in os.listdir(os.path.join(os.getcwd(),'test','test_files')):
        cmd = 'dreem-demultiplexing '\
            +'—-fastq1 \'{}/{}_R1.fastq\' '.format(os.path.join(path_input,sample),sample)\
            +'--fastq2 \'{}/{}_R2.fastq\' '.format(os.path.join(path_input,sample),sample)\
            +'—-library \'{}/library.csv\' '.format(os.path.join(path_input,sample))\
            +'-o ' + path_output
        util.run_cmd(cmd)
        predicted_files_list = os.listdir(os.path.join(os.path.join(path_input,sample),'demultiplexed_fastq'))
        output_files_list  = os.listdir(os.path.join(os.path.join(path_output,sample),'demultiplexed_fastq'))
        assert set(output_files_list) == set(predicted_files_list)