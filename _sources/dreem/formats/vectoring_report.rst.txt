

Vectoring Report
++++++++++++++++

The vectoring report file summarizes the results of vectoring and enables DREEM to locate all mutation vector batch files and verify their integrity.
The report contains the following fields:

- ``Sample Name``: name of the sample from which the reads came
- ``Reference Name``: name of the reference sequence to which the reads were aligned
- ``Section 5' End``: coordinate of the first (5'-most) base in the section, with respect to the full reference sequence
- ``Section 3' End``: coordinate of the last (3'-most) base in the section, with respect to the full reference sequence
- ``Section Length``: number of bases in the section: ``"Section 3' End" - "Section 5' End" + 1``
- ``Section Sequence``: sequence of the section between coordinates ``Section 5' End`` and ``Section 3' End``
- ``Section is Full Reference``: whether the section is the same as the full reference sequence to which the reads were aligned
- ``Reads Vectorized``: number of reads (or pairs of paired-end reads) that were successfully converted into mutation vectors and written to an ORC file
- ``Reads with Errors``: number of reads (or pairs of paired-end reads) that failed to be converted into mutation vectors
- ``Fraction Vectorized``: fraction of total reads that were vectorized successfully: ``"Reads Vectorized" / ("Reads Vectorized" + "Reads with Errors")``
- ``Batches``: number of ORC files written, each containing a batch of mutation vectors 
- ``MD5 Checksums``: comma-separated list of the MD5 checksum of each ORC file (for verifying file integrity)
- ``Began``: local time at which vectoring began, formatted as ``YYYY-MM-DD hh:mm:ss.µs``
- ``Ended``: local time at which vectoring ended, formatted as ``YYYY-MM-DD hh:mm:ss.µs``
- ``Time Taken (s)``: time taken to complete vectorization, in seconds: ``"Ended" - "Began"``
- ``Speed (1/s)``: average number of vectors written per second: ``"Reads Vectorized" / "Time Taken (s)"``

**Example: report.json**

Example contents of a vectoring report file::

    {
     "Sample Name": "DMS_100mM",
     "Reference Name": "SARS2",
     "Section 5' End": "13369",
     "Section 3' End": "13398",
     "Section Length": "30",
     "Section Sequence": "AGTCTGTACCGTCTGCGGTATGTGGAAAGG",
     "Section is Full Reference": "False",
     "Reads Vectorized": "272638",
     "Reads with Errors": "0",
     "Fraction Vectorized": "1.0",
     "Batches": "2",
     "MD5 Checksums": "9f4ffa03d172df618fa664358cc9dc8e, a8f23d7cf38e8a15bcb298037bd23f5c",
     "Began": "2023-03-22 20:00:05.712658",
     "Ended": "2023-03-22 20:00:28.427814",
     "Time taken (s)": "22.715",
     "Speed (1/s)": "12002.6"
    }
