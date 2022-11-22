
from dreem.demultiplexing.demultiplexing import demultiplex
import os
import dreem.util as util
import pandas as pd
from create_test_demultiplexing import generate_demultiplexing_test_samples

## Parameters
length_sequence = 200
barcode_length = 8
min_barcode_distance = 3
directory = 'test/demultiplexing'
num_reads = 1000
num_constructs = 10
barcode_start = 10
sample_name = 'test_reference'
##

def demultiplexed_files_test_fun(df, output_folder, sample_name):
   assert os.path.exists(os.path.join(output_folder, sample_name))
   for i, row in df.iterrows():
      assert os.path.exists(os.path.join(directory, sample_name, row['construct']+'_R1.fastq'))
      assert len(open(os.path.join(directory, sample_name, row['construct']+'_R1.fastq')).readlines()) == num_reads*4
      with open(os.path.join(directory, sample_name, row['construct']+'_R1.fastq')) as f:
         assert f.readline()[:6+barcode_length] == '@READ:'+row['construct']
         line = f.readline()
         assert len(line) == length_sequence+1
         assert line[row['barcode_start']:row['barcode_start']+barcode_length] == row['barcode_sequence']
         assert f.readline() == '+\n'
         assert len(f.readline()) == length_sequence+1
         f.close()
      assert os.path.exists(os.path.join(directory, sample_name, row['construct']+'_R2.fastq'))
      assert len(open(os.path.join(directory, sample_name, row['construct']+'_R2.fastq')).readlines()) == num_reads*4

def test_generate_test_sample():
   # generate a test sample
   os.system('rm -fr {}'.format(directory))
   os.mkdir(directory)
   os.mkdir(os.path.join(directory, sample_name))
   assert generate_demultiplexing_test_samples(directory=directory, num_reads=num_reads, num_constructs=num_constructs, barcode_length=barcode_length, min_barcode_distance=min_barcode_distance, barcode_start=barcode_start, length_sequence=length_sequence, sample_name=sample_name),\
      'Failed to generate test sample'
   assert os.path.exists(os.path.join(directory, f'{sample_name}_R1.fastq'))
   assert os.path.exists(os.path.join(directory, 'library.csv'))
   df = pd.read_csv(directory+'/library.csv')
   assert len(df) == num_constructs
   demultiplexed_files_test_fun(df, directory, sample_name)
   assert os.path.exists(os.path.join(directory, '{}_R2.fastq'.format(sample_name)))
   with open(os.path.join(directory, '{}_R1.fastq'.format(sample_name))) as f:
      assert len(f.readlines()) == num_reads*num_constructs*4
   with open(os.path.join(directory, '{}_R2.fastq'.format(sample_name))) as f:
      assert len(f.readlines()) == num_reads*num_constructs*4


def test_demultiplexing():
   fastq1, fastq2 = os.path.join(directory, '{}_R1.fastq'.format(sample_name)), os.path.join(directory, '{}_R2.fastq'.format(sample_name))
   library = pd.read_csv(os.path.join(directory, 'library.csv'))
   output_folder = os.path.join(directory, 'test_output')
   temp_folder = os.path.join(directory, 'test_temp')
   os.system('rm -fr {}'.format(output_folder))
   os.mkdir(output_folder)
   os.system('rm -fr {}'.format(temp_folder))
   os.mkdir(temp_folder)
   assert demultiplex(fastq1, fastq2, library, output_folder, temp_folder), 'Demultiplexing failed'
   demultiplexed_files_test_fun(library, directory, 'test_output')

   
