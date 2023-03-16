
I/O files
++++++++++++++++++++++++

*This module needs the files to be sorted into a specific folder structure. The folder structure is described below.*

**Input file structure**

The input is one or more files of aligned reads in binary alignment map (BAM) format. Each file must be named after the reference sequence to which it was aligned, end with the file extension ``.bam``, and reside in a directory named after the sample from which the reads came::

    |- sample_1/
        |- reference_1.bam
        |- reference_2.bam
    |- sample_2/
        |- reference_1.bam
        |- reference_3.bam


**Output file structure**

The output is structured as follows::

    |- sample_1/
        |- reference_1/
            |- section_1/
                |- 0.orc
                |- 1.orc
                |- report.txt
            |- section_2
                |- 0.orc
                |- report.txt
        |- reference_2
            |- section_3
                |- 0.orc
                |- 1.orc
                |- 2.orc
                |- report.txt
