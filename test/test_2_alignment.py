import dreem, os
import dreem.util as util
import pandas as pd

path_input = os.path.join(os.getcwd(),'test','test_files','input','alignment')
path_predicted = os.path.join(os.getcwd(),'test','test_files','predicted_output','alignment')
path_output = os.path.join(os.getcwd(),'test','test_files')

def test_files_exists():
    for sample in os.listdir(path_input):
        cmd = 'dreem-alignment '\
            +'—-fastq1 \'{}/{}_R1.fastq\' '.format(os.path.join(path_input,sample),sample)\
            +'--fastq2 \'{}/{}_R2.fastq\' '.format(os.path.join(path_input,sample),sample)\
            +'—-fasta \'{}/reference.fasta\' '.format(os.path.join(path_input,sample))\
            +'-o ' + path_output\
            +'--sample '+sample
        util.run_cmd(cmd)
        
        assert os.path.exists(os.path.join(path_output,'output','alignment',sample)), 'Output folder for sample {} doesn\'t exist'.format(os.path.join(path_output,'output','alignment',sample))
        for bam_file in os.listdir(os.path.join(path_predicted,sample)):
            assert os.path.isfile(os.path.join(path_output,'output','alignment',sample,bam_file)), 'File {} is missing'.format(bam_file)

