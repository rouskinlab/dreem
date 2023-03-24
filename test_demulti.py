
import os
import sys
 
#sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
import os
from dreem.main import run


if __name__ == '__main__':
    print("/Users/scottgrote/Documents/absolutely_final_repo/dreem_data/short_R1.fastq".split("/")[-1].split("_R")[0])
    """run(
        fasta="/Users/scottgrote/Documents/absolutely_final_repo/dreem_data/new_ref.fasta",
        fastq1="/Users/scottgrote/Documents/absolutely_final_repo/dreem_data/short_R1.fastq",
        fastq2="/Users/scottgrote/Documents/absolutely_final_repo/dreem_data/short_R2.fastq",
        library="/Users/scottgrote/Documents/absolutely_final_repo/dreem_data/updated_library.csv",
        clipped=2,
        mismatch_tolerence=2,
        index_tolerance=20,
        parallel_demultiplexing=True,
        demult_on=True,
        demulti_overwrite=False,
        )"""