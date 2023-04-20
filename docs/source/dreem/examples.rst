
Examples
========

    Here are some examples to get started with the dreem pipeline:
    .. note::
        It is assumed that dreem has been installed in your conda enviroment.


CLI
---------
    The following command will run dreem on some testfiles provided in the repo:

        ``dreem --fasta docker-data/3509.fasta --fastq1 docker-data/3509-O-flank_1=bi1-ms2-DB_R1.fastq --fastq2 docker-data/3509-O-flank_1=bi1-ms2-DB_R2.fastq --verbose --autosect --out-dir docker-data/out --temp-dir docker-data/temp``
    
    Command explained:

        * **dreem** - conda's entry point executable for the dreem program
        * **--fasta** - path to the reference fasta (reference sequences)
        * **--fastq1** - path to 1 of 2 of two files containing sequnece data in the form if "reads" 
        * **--fastq2** - path to 1 of 2 of two files containing sequnece data in the form if "reads" 
        * **--verbose** - forces dreem to output info on real-time processes
        * **--out-dir** - path for which dreem will write outputs of each step of analysis
        * **--temp-dir** - path for which dreem will write outputs of each step of analysis ( see documentation of modules I/O for more information on what files each module stores in either the temp or out dir)
    
    Since these test files are directly in the repo this command should be able to be run as is, if run from the repo directory



Python
------------

