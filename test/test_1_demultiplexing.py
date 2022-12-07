
import os
import dreem.util as util
import pandas as pd

from dreem.demultiplexing.run import run as demultiplexing_run

path_input = os.path.join(os.getcwd(),'test','test_files','input','demultiplexing')
path_predicted = os.path.join(os.getcwd(),'test','test_files','predicted_output','demultiplexing')
path_output = os.path.join(os.getcwd(),'test','test_files')

def test_files_exists():
    for sample in os.listdir(path_input):
        
        args = {
            'fastq1': '{}/{}_R1.fastq'.format(os.path.join(path_input,sample),sample),
            'fastq2': '{}/{}_R2.fastq'.format(os.path.join(path_input,sample),sample),
            'library': '{}/library.csv'.format(os.path.join(path_input,sample)),
            'out_dir': path_output
        }
        
        #demultiplexing_run(args)
        
        cmd = 'dreem-demultiplexing '\
            +'—-fastq1 \'{}/{}_R1.fastq\' '.format(os.path.join(path_input,sample),sample)\
            +'--fastq2 \'{}/{}_R2.fastq\' '.format(os.path.join(path_input,sample),sample)\
            +'—-library \'{}/library.csv\' '.format(os.path.join(path_input,sample))\
            +'-o ' + path_output
        util.run_cmd(cmd)
        print(cmd)
        
        assert os.path.exists(os.path.join(path_output,'output','demultiplexing',sample)), 'Output folder for sample {} doesn\'t exist'.format(os.path.join(path_output,'output','demultiplexing',sample))
        for demultiplexed_file in os.listdir(os.path.join(path_predicted,sample)):
            assert os.path.isfile(os.path.join(path_output,'output','demultiplexing',sample,demultiplexed_file)), 'File {} is missing'.format(demultiplexed_file)

def test_all_files_are_equal():
    for sample in os.listdir(path_input):
        for pred, out in zip(os.listdir(os.path.join(path_predicted,sample)), os.listdir(os.path.join(path_output,'output','demultiplexing',sample))):
            assert pred == out, 'The predicted output and the output files are not the same'
            predicted = util.fastq_to_df(os.path.join(path_predicted,sample,pred))
            predicted['from'] = 'predicted'
            output = util.fastq_to_df(os.path.join(path_output,'output','demultiplexing',sample,out))
            output['from'] = 'output'
            both = pd.concat([predicted,output],axis=0, ignore_index=True)
            for idx, g in both.groupby('construct'):
                if len(g) < 2:
                    assert g['from'] == 'predicted', 'The output file is missing the construct {} for sample {}'.format(idx,sample)
                    assert g['from'] == 'output', 'The output file didn\'t filter out the construct {} for sample {}'.format(idx,sample)
            for idx, g in both.groupby('construct'):
                for c in both.columns:
                    if c != 'construct' and c != 'from':
                        assert g[c].unique().shape[0] == 1, 'The output file is not the same as the predicted output for sample {} and construct {}'.format(sample,idx)
                        
