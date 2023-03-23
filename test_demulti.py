
import os
import sys
 
#sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
import os
from dreem.main import run


if __name__ == '__main__':

    run(
        fasta="/Users/scottgrote/Documents/merging_matty/test_files/ref.fasta",
        fastq1="/Users/scottgrote/Documents/absolutely_final_repo/lauren473_S4_R1.fastq",
        fastq2="/Users/scottgrote/Documents/absolutely_final_repo/lauren473_S4_R2.fastq",
        library="/Users/scottgrote/Documents/absolutely_final_repo/library_test.csv",
        clipped=0,
        mismatch_tolerence=0,
        index_tolerance=0,
        parallel_demultiplexing=True,
        demult_on=True,
        demulti_overwrite=True
        )