
Key algorithm(s)
++++++++++++++++

Create SAM files to vectorize
-----------------------------

Purpose: In order to turn BAM files into mutation vectors, the BAM files must be converted into SAM files containing only reads that align to the references and sections of interest.

1. For each pair of 5'/3' coordinates (``--coords``) or forward/reverse primers (``--primers``) given for every reference, a *section* (i.e. subsequence between two coordinates) of its reference is defined.
2. If the option ``--autosect`` (CLI) or ``autosect=True`` (API) is given, then for every reference for which no coordinates/primers were given, a section equal to the full reference sequence is defined automatically.
3. For each given BAM file (i.e. alignment of a sample to a reference sequence), and for each section defined for its reference sequence, a temporary SAM file of only the reads aligning to the section is generated using ``samtools``.

In the following steps, each temporary SAM file is then processed independently, either in parallel (with Python's ``multiprocessing`` module) or in series.

Define batches of reads
-----------------------

Purpose: memory and parallelization

1.

Parse CIGAR strings into vectors
--------------------------------

Find ambiguous mutations
------------------------


Write output files
------------------
