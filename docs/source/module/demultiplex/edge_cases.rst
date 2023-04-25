
Edge case handling
+++++++++++++++++++

**Index Tolerence:**
--------------------

Giving a index tolerence on index too close to the start or end of the sequence will result in a error if not handled. 
This currently is handled by taking the value that would go below zero and instead adding it to the oppostie side of the index tolerence

**Fastq cases:**
----------------

Giving two Fastqs that do not have the same number of reads
    * this can be because they are corrupted
    * erroronoeus input 
        - fastqs from different samples
        - incorrect file path

    * these issues are dealth with with in the object super_fastq

**Mutliple Barcode:**
----------------
If barcodes are not desgined well enough and they appear mutliple times in the sequence, it will cause issues



