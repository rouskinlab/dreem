
Get started
++++++++++++++++


Let's turn the bitvectors of ``sample_1`` into the :ref:`dreem_output` format. 
You will need the reference genome under the `fasta format <https://en.wikipedia.org/wiki/FASTA_format>`_ ``reference.fasta``.

Store the bitvectors as follow:

::

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




CLI
---------

Run the command:

::
    
    dreem-aggregate --input_dir path/to/sample_1 —fa path/to/reference.fasta 


You get the following output:
 
::

    sample_1.json


See the :ref:`aggregate_reference` for more details.


Python
------------

Run the following code:

:: 

    from dreem import aggregate
    aggregate.run(
        input_dir='path/to/sample_1', 
        fasta='path/to/reference.fasta'
        )
        

You get the following output:
 
::

    sample_1.json


See the :ref:`aggregate_reference` for more details.
