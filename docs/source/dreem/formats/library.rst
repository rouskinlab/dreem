
Library
+++++++++++++++

Add information related to a specific reference from the fasta file by writing a library file, which is a csv file with the following columns:

 ========================== =============================================================================== 
  column                     description                                                                    
 ========================== =============================================================================== 
  construct                  reference name                                                             
  barcode_start (optional)   start position of the barcode (1-based)                                    
  barcode_end (optional)     end position of the barcode (1-based inclusive)                            
  section (optional)         section of the construct                                                   
  section_start (optional)   start position of the section (1-based)                                    
  section_end (optional)      end position of the section  (1-based inclusive)                           
  [any other columns]        any other columns will be added to corresponding construct in the output file  
 ========================== =============================================================================== 


.. note::

    If one construct has multiple sections, the library file should contain multiple rows for the same construct.
    Make sure that a construct doesn't have multiple values for the other columns.

**Example: library.csv**

+-----------------------+---------------+-------------+----------+---------------+-------------+--------+
| construct             | barcode_start | barcode_end | section  | section_start | section_end | family |
+=======================+===============+=============+==========+===============+=============+========+
| 3042-O-flank_1=hp1-DB | 140           | 151         | section1 | 69            | 86          | hp1    |
+-----------------------+---------------+-------------+----------+---------------+-------------+--------+
| 3042-O-flank_1=hp1-DB |               |             | section2 | 20            | 42          | hp1    |
+-----------------------+---------------+-------------+----------+---------------+-------------+--------+