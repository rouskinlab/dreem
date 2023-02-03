
# Library

The library is a csv file containing per-construct information. The columns are:
- construct: the construct name
- barcode (optional): the barcode sequence
- barcode_start (optional): the start position of the barcode (0-based)
- section (optional): the section of the construct
- section_start (optional): the start position of the section (0-based)
- section_end (optional): the end position of the section  (0-based, exclusive)
- [any other columns]: any other columns will be added to corresponding construct in the output file

If one construct has multiple sections, the library file should contain multiple rows for the same construct.
In that case, make sure that a construct doesn't have multiple inputs for the other columns (i.e multiple barcodes).

