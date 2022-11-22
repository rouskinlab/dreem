
import os

def run():
    os.system("""source venv/bin/activate 
                make
               dreem -cd 1 -fa 'test/test_sample/ref.fasta' -dx True -cl True -fq1 'test/test_sample/my_sample_1_R1.fastq' -fq2 'test/test_sample/my_sample_1_R2.fastq' -l 'test/test_sample/library.csv' -s 'test/test_sample/samples.csv' -v True""")


if __name__ == '__main__':
    run()
