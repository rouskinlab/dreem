import os, sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import os
import dreem

test_files = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'test_files')
test_files_sample = os.path.join(test_files, 'my_sample_to_demultiplex')
demultiplexed_fastq_dir = os.path.join(test_files_sample,'my_sample_to_demultiplex')
library = os.path.join(test_files_sample, 'library.csv')
fasta =  os.path.join(test_files_sample, 'my_sample_to_demultiplex.fasta')
fastq1 = os.path.join(test_files_sample, 'my_sample_to_demultiplex_R1.fastq')
fastq2 = os.path.join(test_files_sample, 'my_sample_to_demultiplex_R2.fastq')
out_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'test_output','demultiplexing','my_sample_to_demultiplex')

def test_demultiplexing_python():
    dreem.demultiplex.run(
        library=library,
        fastq1=fastq1,
        fastq2=fastq2,
        out_dir=out_dir,
        fasta=fasta,
    )

def test_demultiplexing_cli():
    os.system('dreem demultiplex --library {} --fasta {} --fastq1 {} --fastq2 {} --out-dir {} --rerun'.format(
        library,
        fasta,
        fastq1,
        fastq2,
        out_dir,
    ))
    
if __name__ == '__main__':
    test_demultiplexing_python()
    test_demultiplexing_cli()