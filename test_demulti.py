
import os
import sys
 
#sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
import os
import pytest
from dreem.main import run


if __name__ == '__main__':
    run(
        fastq1="/Users/scottgrote/Documents/hopefully_final_repo/test_files/lauren473_S4_R1.fastq",
        fastq2="/Users/scottgrote/Documents/hopefully_final_repo/test_files/lauren473_S4_R2.fastq",
        library="/Users/scottgrote/Documents/hopefully_final_repo/test_files/library_test.csv",
        clipped=0,
        mismatch_tolerence=0,
        index_tolerance=0,
        parallel_demultiplexing=False,
        demult_on=True,
        )

"""
run(
        fastq1="/Users/scottgrote/Documents/hopefully_final_repo/test_files/lauren473_S4_R1_001.fastq",
        fastq2="/Users/scottgrote/Documents/hopefully_final_repo/test_files/lauren473_S4_R2_001.fastq",
        library="/Users/scottgrote/Documents/hopefully_final_repo/test_files/library.csv",
        clipped=0,
        mismatch_tolerence=0,
        index_tolerance=0,
        parallel_demultiplexing=True,
        )

given this input, when the fastq1 given doesnt exist it gives, Segment '' failed to match pattern, 

but to make it more clear we should also check if the files exists


"""
