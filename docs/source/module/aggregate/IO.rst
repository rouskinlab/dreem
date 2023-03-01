
I/O files
++++++++++++++++++++++++

*This module needs the files to be sorted into a specific folder structure. The folder structure is described below.*

**Input files structure**

The input file type is a batch of :ref:`bitvector` files, one per section of a reference. The folder structure is as follows::

    sample_1/              # <=> a fastq file
        |- report.txt
        |- reference_1/    # <=> a single reference in the fasta file
            |- section_1/  # <=> a section (a sub-sequence of a reference) 
                |- 0.orc   # <=> a batch of bitvectors
                |- 1.orc
            |- section_2  # <=> another section for this reference
                |- 0.orc
        |- reference_2    # <=> another reference
            |- section_3
                |- 0.orc
                |- 1.orc
                |- 2.orc



**Output file structure**

:: 

    sample_1.json

