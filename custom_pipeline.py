import dreem, os
import pandas as pd
import sys
import glob

if __name__ == '__main__':
    # Parameters
    #assert len(sys.argv) == 2, 'Please provide a sample name (e.g. 01_1_S22_reads)'
    #sample = sys.argv[1]
    
    platform = 'Matty' # 'Yves' or 'O2'
    
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
    if platform == 'Scott':
        # Yves' computer
        root_dir = '/Users/scottgrote/Documents/ultimate_dreem_repo/dreem/dreem_data/'
        samples = root_dir + 'samples.csv'
        rnastructure_path = '/Applications/RNAstructure/exe'
    if platform == "Matty":
        # Matty's computer
        root_dir = "/Users/mfa/git/dreem/test_data/"
        samples = root_dir + "samples.csv"
    
    out_dir = root_dir
    lib = root_dir + "library.csv"
    library = root_dir + 'library.csv'
    fasta = root_dir + 'reference.fasta'
    sample="lauren474_S5"
    
    verbose = True

    big_fastq1="/Users/scottgrote/Documents/ultimate_dreem_repo/dreem/dreem_data/test_inputs/lauren474_S5_R1.fastq"#sys.argv[1]
    big_fastq2="/Users/scottgrote/Documents/ultimate_dreem_repo/dreem/dreem_data/test_inputs/lauren474_S5_R2.fastq"#sys.argv[2]
    

    def verbose_print(*args):
        print(*args) if verbose else None

    # make output and temp folders
    for folder in ['output', 'temp']:

        os.makedirs(os.path.join(out_dir, folder), exist_ok=True)

    run_demulti=False
    run_alignment=True
    run_vectoring=True
    run_aggregrate=True

    # Run DREEM
    verbose_print("""

    ========================================

                RUNNING   DREEM

    ========================================

    """)
    ## demulti: 
    # -----------------------------------------------------------------------------------------------------------------------
    if run_demulti: 
        verbose_print('\ndemulti \n----------------')
        #print(vars(dreem.alignment._))
        #print(dir(dreem.demultiplex.main))
        demultiplexed_reads_dir=dreem.demultiplexing.main.run(
                            top_dir=os.path.join(out_dir),#
                            fastq1=big_fastq1,
                            fastq2=big_fastq2,
                            library=lib
                            )

        fastq1 = glob.glob(demultiplexed_reads_dir+"*_R1.fastq")
        fastq2 = glob.glob(demultiplexed_reads_dir+"*_R2.fastq")
        fastq1.sort(), fastq2.sort()
        # -----------------------------------------------------------------------------------------------------------------------
    else:
        pass
        #fastq1 = glob.glob(demultiplexed_reads_dir+"*_R1.fastq")
        #fastq2 = glob.glob(demultiplexed_reads_dir+"*_R2.fastq")
        #fastq1.sort(), fastq2.sort()
    demultiplexed_reads_dir= root_dir + "output/demultiplexing/" + sample
    
    ## Alignment: 
    # -----------------------------------------------------------------------------------------------------------------------


    if run_alignment:
        verbose_print('\nalignment \n----------------')
        
        #print(f1)
        #print(f2)
        dreem.alignment.run(
                            top_dir=os.path.join(out_dir),#, 'output','alignment'),
                            fasta=fasta,
                            fastq12_dir=demultiplexed_reads_dir,
                            )
            # -----------------------------------------------------------------------------------------------------------------------                       

    ## Vectoring
    # -----------------------------------------------------------------------------------------------------------------------
    if run_vectoring:
        verbose_print('\nvectoring \n------------------')
        path_to_bam = os.path.join(out_dir, 'output', 'alignment', sample)
        bam_dirs = [os.path.join(path_to_bam, f) for f in os.listdir(path_to_bam) if f.endswith('.bam')]
        dreem.vectoring.run(
                        out_dir= os.path.join(out_dir, 'output'), #TODO
                        bam_dirs= bam_dirs,
                        fasta=fasta,
                        library=library,
                        parallel="profiles"
                        )
    # -----------------------------------------------------------------------------------------------------------------------



    ## Aggregate
    # -----------------------------------------------------------------------------------------------------------------------
    if run_aggregrate:
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
                    #  poisson=True)"""

    verbose_print('Done!')
    # -----------------------------------------------------------------------------------------------------------------------
