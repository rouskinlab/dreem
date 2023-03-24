
Examples
++++++++

CLI
---

Basic usage showing how to create mutation vectors from the BAM file ``DMS_100mM/MAPT-3R.bam``
over section 231 - 417 of the transcript named "MAPT-3R",
with the reference sequence defined in ``human_transcriptome.fasta``::

    dreem vector -b DMS_100mM/MAPT-3R.bam -f human_transcriptome.fasta -c MAPT-3R 231 417


Python
------

Basic usage of the Python API to perform the same operation as shown in the CLI example.

>>> from dreem import vector
>>> report_files = vector.main.run(bamf=("DMS_100mM/MAPT-3R.bam",),
...                                fasta="human_transcriptome.fasta",
...                                coords=(("MAPT-3R", 231, 417),))
