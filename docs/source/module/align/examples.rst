Examples
++++++++

CLI
---

Basic usage of the align module, showing how to align a pair of FASTQ files containing paired-end reads
(``DMS_100mM_R1.fq.gz`` and ``DMS_100mM_R2.fq.gz``) to a set of reference sequences in ``human_transcriptome.fasta``::

    dreem align -1 DMS_100mM_R1.fq.gz -2 DMS_100mM_R2.fq.gz -f human_transcriptome.fasta


Python
------

Basic usage of the Python API, showing the same example as with the CLI:

>>> from dreem import align
>>> bam_files = align.main.run(fastq1=("DMS_100mM_R1.fq.gz",),
...                            fastq2=("DMS_100mM_R2.fq.gz",),
...                            fasta="human_transcriptome.fasta")
