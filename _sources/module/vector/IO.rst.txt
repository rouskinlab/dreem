
I/O files
++++++++++++++++++++++++

*This module needs the files to be sorted into a specific folder structure. The folder structure is described below.*

**Input file structure**

The module takes one or more files of aligned reads in binary alignment map (BAM) format. Each file must be named after the reference sequence to which it was aligned, end with the file extension ``.bam``, and reside in a directory named after the sample from which the reads came::

    |- sample_1/
        |- reference_1.bam
        |- reference_2.bam
    |- sample_2/
        |- reference_1.bam
        |- reference_3.bam


**Output file structure**

The module outputs one report file (`JSON format <https://www.json.org/json-en.html>`_) and one or more mutation vector batch files (`Apache ORC format <https://orc.apache.org/docs>`_) for each section of each reference for each sample. Each sample corresponds to one directory that contains one sub-directory for each reference to which the sample was aligned. Each reference directory contains one sub-directory for each section of the reference for which mutation vectors were computed. Each section directory contains one or more mutation vector batch files (numbered ``0.orc``, ``1.orc``, ...) and a report file (``report.json``) that contains metadata about the vectors::

    |- sample_1/
        |- reference_1/
            |- section_1.1/
                |- 0.orc
                |- 1.orc
                |- report.json
            |- section_1.2
                |- 0.orc
                |- report.json
        |- reference_2
            |- section_2.1
                |- 0.orc
                |- 1.orc
                |- 2.orc
                |- report.json
    |- sample_2/
        |- reference_1/
            |- section_1.1/
                |- 0.orc
                |- report.json
            |- section_1.2
                |- 0.orc
                |- 1.orc
                |- report.json
        |- reference_3/
            |- section_3.1/
                |- 0.orc
                |- report.json
