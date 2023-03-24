
import os
import sys

#sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
import os
from dreem.main import run


if __name__ == '__main__':

    run(
        fasta="/n/data1/hms/microbiology/rouskin/lab/projects/TheJuicyOrange/final_dreem_repo/test_files/all_ref.fasta",
        fastq1=sys.argv[1],
        fastq2=sys.argv[2],
        library="/n/data1/hms/microbiology/rouskin/lab/projects/TheJuicyOrange/final_dreem_repo/test_files/library_full.csv",
        clipped=2,
        mismatch_tolerence=2,
        index_tolerance=20,
        parallel_demultiplexing=True,
        demult_on=True,
        demulti_overwrite=False,
	out_dir="/n/data1/hms/microbiology/rouskin/lab/projects/TheJuicyOrange/final_dreem_repo/test_files/"
        )
