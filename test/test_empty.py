import os, sys
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import os
import pytest
from dreem.main import run
from dreem.util.files_sanity import compare_fields

test_files = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'test_files')
test_files_sample = os.path.join(test_files, 'my_empty_sample')
library = os.path.join(test_files_sample, 'library.csv')
fasta =  os.path.join(test_files_sample, 'my_empty_sample.fasta')
fastq1 = os.path.join(test_files_sample, 'my_empty_sample_R1.fastq')
fastq2 = os.path.join(test_files_sample, 'my_empty_sample_R2.fastq')
out_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'test_output','empty')



def test_run_cli():
    os.system("dreem --library {} --fasta {} --fastq1 {} --fastq2 {} --out-dir {} --rerun".format(
        library,
        fasta,
        fastq1,
        fastq2,
        os.path.join(out_dir, 'cli'),
    ))


def test_run_python():
    run(
        library=library,
        fasta=fasta,
        fastq1=(fastq1,), 
        fastq2=(fastq2,),
        out_dir=os.path.join(out_dir, 'python'),
    )


# what should be the output of an empty sample?
@pytest.mark.skip(reason="I don't know what the output of an empty sample is.")
def test_no_output():  
    pass


if __name__ == '__main__':
    test_run_cli()
    test_run_python()
    test_no_output()