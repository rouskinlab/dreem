

Bitvector
++++++++++++++++

The bitvector is storage format for the aligned reads. The first line of the bitvector is the read sequence. The following lines are the aligned reads, each line being a read.

Each base if the read is represented by a byte. 
The meaning of each byte is as follows:


 ========== ========================================================= 
  Bit        Meaning                                                  
 ========== ========================================================= 
  00000001   match                                                    
  00000010   deletion                                                 
  00000100   insertion of ≥1 base(s) immediately 3' of this position  
  00001000   insertion of ≥1 base(s) immediately 5' of this position  
  00010000   substitution to A                                        
  00100000   substitution to C                                        
  01000000   substitution to G                                        
  10000000   substitution to T                                        
 ========== ========================================================= 

The bitvector is stored under Apache's .orc format. 
See more about the vectorization algorithm :ref:`here <vector_algos>`.

**Example**

Using the following reference and a single aligned read:

.. code-block:: text

    Reference:  CACCATTCG
    Read:      CTACTA?TC

Bitvector:

.. code-block:: text

        C   A   C   C   A   T   C   T   G  # reference sequence 
        9   5   1   128 1   0   1   1   2  # bitvector


More detailed computation of the bitvector:

.. code-block:: text

    Reference:   C   A   C   C   A   T   C   T   G  
    Read:        C T A   C   T   A   ?   C   T   -  (? = couldn't align, - = deletion)   
    Bitvector:   9   5   1   128 1   0   1   1   2

