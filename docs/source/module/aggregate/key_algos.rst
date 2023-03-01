
Key algorithms
++++++++++++++

Counting mutations for the entire bitvector
-------------------------------------------

The following table describes the aggregation of the bitvector into a set of per-residue counts (except for ``num_of_mutations`` which is per-read). 

Please refer to :ref:`bitvector` for the definition of the bits and the bitvector format.


======================= =========================================================================================================================
**Aggregation**          **Description**
----------------------- -------------------------------------------------------------------------------------------------------------------------
``cov_bases``            Per-residue count of covered bases
``del_bases``            Per-residue count of deleted bases
``ins_bases``            Per-residue count of inserted bases
``info_bases``           Per-residue count of bases that are substituted for another base or non-mutated. Doesn't include deleted bases. 
``mod_bases_A``          Per-residue count of bases that substituted to A
``mod_bases_C``          Per-residue count of bases that substituted to C
``mod_bases_G``          Per-residue count of bases that substituted to G
``mod_bases_T``          Per-residue count of bases that substituted to T
``mod_bases_N``          Per-residue count of bases that substituted to any base
``mut_bases``            Per-residue count of bases that substituted to any base
``mut_rates``            mut_bases / info_bases for each residue
``num_of_mutations``     Per-read count of mutations histogram (i.e. number of reads with 0 mutations, 1 mutation, 2 mutations, etc.)
======================= =========================================================================================================================



**Example**:


.. note::
    
    @Matty please read this and validate it!

.. code-block:: bash

    # example.orc
                    # Bitvector 
                    C    A    C    A    A    T    G    T    G   # reference sequence 
                    9    5    1    128  1    0    1    1    2   # read 1
                    1    1    128  1    0    1    16   2    1   # read 2 
                    64   1    2    1    1    1    16   1    1   # read 3
    # Aggregation ---------------------------------------------
    cov_bases       3    3    3    3    2    2    3    3    3
    del_bases       0    0    0    0    0    0    0    1    1
    ins_bases       0    1    0    0    0    0    0    0    0  
    info_bases      3    3    2    3    2    2    3    2    2
    mod_bases_A     0    0    0    0    0    0    2    0    0
    mod_bases_C     0    0    0    0    0    0    0    0    0
    mod_bases_G     1    0    0    0    0    0    0    0    0
    mod_bases_T     0    0    1    1    0    0    0    0    0
    mod_bases_N     1    0    1    1    0    0    2    0    0
    mut_bases       1    0    1    1    0    0    2    0    0
    mut_rates       0.33 0.00 0.5  0.33 0.00 0.00 0.66 0.00 0.00

    # This one is per read, not per base
    num_of_mutations   0  2  1  0  0  0  0  0  0  0  



Counting mutations for a certain cluster
----------------------------------------

[To be implemented]


Predicting the structure and the free energy
--------------------------------------------

We use the following script:

.. code-block:: text

    echo ">ref\nGGCGACACAGTCGACGGTTTTCACA">GGCGACACAGTCGACGGTTTTCACA.fasta
    Fold GGCGACACAGTCGACGGTTTTCACA.fasta GGCGACACAGTCGACGGTTTTCACA.ct
    ct2dot GGCGACACAGTCGACGGTTTTCACA.ct 1 GGCGACACAGTCGACGGTTTTCACA_dot.txt
    cat GGCGACACAGTCGACGGTTTTCACA_dot.txt

.. note::

    The files are named after the sequence, so you can use the same result files for the same sequence amongst different runs.    