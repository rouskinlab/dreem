import dreem, os
import pandas as pd
import sys

if __name__ == '__main__':
    # Parameters
    assert len(sys.argv) == 2, 'Please provide a sample name (e.g. 01_1_S22_reads)'
    sample = sys.argv[1]
    
    platform = 'O2' # 'Yves' or 'O2'
    
    if platform == 'Yves':
        # Yves' computer
        root_dir = '/Users/ymdt/src/data_Lauren/'
        samples = root_dir + 'samples.csv'
        rnastructure_path = '/Users/ymdt/src/RNAstructure/exe'

    if platform == 'O2':
        # O2
        root_dir = '/n/data1/hms/microbiology/rouskin/lab/projects/TheJuicyOrange/Sarah_demulti/to_run/dreem_runs/'
        samples = root_dir + 'samples.csv'
        rnastructure_path = '/n/data1/hms/microbiology/rouskin/lab/lib/RNAstructure/exe'    
    
    out_dir = root_dir
    fasta = root_dir +'reference.fasta'
    library = root_dir + 'library.csv'

    verbose = True
    fastq1 = [os.path.join(root_dir + sample, f) for f in os.listdir(root_dir + sample) if f.endswith('R1.fastq')]
    fastq2 = [os.path.join(root_dir + sample, f) for f in os.listdir(root_dir + sample) if f.endswith('R2.fastq')]
    fastq1.sort(), fastq2.sort()    
    
    def verbose_print(*args):
        print(*args) if verbose else None

    # make output and temp folders
    for folder in ['output', 'temp']:
        os.makedirs(os.path.join(out_dir, folder), exist_ok=True)

    # Run DREEM
    verbose_print("""

    ========================================

                RUNNING   DREEM

    ========================================

    """)
    

    ## Alignment: 
    # -----------------------------------------------------------------------------------------------------------------------
    verbose_print('\nalignment \n----------------')
    for f1, f2 in zip(fastq1, fastq2):
        dreem.alignment.run(
                        out_dir=os.path.join(out_dir),#, 'output','alignment'),
                        fasta=fasta,
                        fastq=f1,
                        fastq2=f2,
                        demultiplexed=True
                        )
        # -----------------------------------------------------------------------------------------------------------------------

    ## Vectoring
    # -----------------------------------------------------------------------------------------------------------------------
    verbose_print('\nvectoring \n------------------')
    path_to_bam = os.path.join(out_dir, 'output', 'alignment', sample)
    bam_files = [os.path.join(path_to_bam, f) for f in os.listdir(path_to_bam) if f.endswith('.bam')]
    dreem.vectoring.run(
                    out_dir= os.path.join(out_dir, 'output'), #TODO
                    bam_files= bam_files,
                    fasta=fasta,
                    library=library,
                    )
    # -----------------------------------------------------------------------------------------------------------------------

    ## Aggregate
    # -----------------------------------------------------------------------------------------------------------------------
    verbose_print('\naggregation \n------------------')
    vect_out_dir = os.path.join(out_dir, 'output', 'vectoring', sample)
    constructs = [c for c in os.listdir(vect_out_dir) if os.path.isdir(os.path.join(vect_out_dir, c))]
    sections = [[s for s in os.listdir(os.path.join(vect_out_dir, c)) if os.path.isdir(os.path.join(vect_out_dir, c, s))] for c in constructs]
    dreem.aggregation.run(
                    bv_files= [os.path.join(out_dir,'output','vectoring', sample, c, s) for c in constructs for s in sections[constructs.index(c)]],
                    fasta=fasta,
                    library= library,
                    sample = sample, 
                    samples= samples,
                    out_dir= os.path.join(out_dir),
                    rnastructure_path=rnastructure_path)
                #  poisson=True)

    verbose_print('Done!')
    # -----------------------------------------------------------------------------------------------------------------------
