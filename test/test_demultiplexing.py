import os, sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import os
import pytest
import dreem

test_files = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'test_files')
test_files_sample = os.path.join(test_files, 'my_sample_to_demultiplex')
demultiplexed_fastq_dir = os.path.join(test_files_sample,'my_sample_to_demultiplex')
library = os.path.join(test_files_sample, 'library.csv')
fasta =  os.path.join(test_files_sample, 'my_sample_to_demultiplex.fasta')
fastq1 = os.path.join(test_files_sample, 'my_sample_to_demultiplex_R1.fastq')
fastq2 = os.path.join(test_files_sample, 'my_sample_to_demultiplex_R2.fastq')
out_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'test_output','demultiplexing')
temp_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'test_temp','demultiplexing')

def test_demultiplexing_python():
    dreem.demultiplex.run(
        library=library,
        fastq1=(fastq1,),
        fastq2=(fastq2,),
        clipped=2,
        index_tolerance=10,
        mismatch_tolerence=2,
        out_dir=os.path.join(out_dir, 'python'),
        temp_dir=os.path.join(out_dir, 'python'),
        fasta=fasta,
    )

def test_demultiplexing_cli():
    os.system(f"dreem demultiplex --library {library} --fasta {fasta} --fastq1 {fastq1} --fastq2 {fastq1} --out-dir {os.path.join(out_dir, 'cli')} --temp-dir {os.path.join(temp_dir, 'cli')} --index-tolerance {10} --mismatch-tolerence {2} ")
    
@pytest.mark.parametrize('sample', ['cli', 'python'])
def test_results(sample):
    directory = os.path.join(out_dir, sample, 'my_sample_to_demultiplex')
    for file in os.listdir(directory):
        lines = open(os.path.join(directory, file)).readlines()
        assert len(lines) == 4*9 or len(lines) == 4*9+1, 'File {} has {} lines, expected {}'.format(file, len(lines), 4*10) # 4 lines per read, 10 reads per sample, the 10th read is garbage
    
if __name__ == '__main__':
    test_demultiplexing_python()
    test_demultiplexing_cli()