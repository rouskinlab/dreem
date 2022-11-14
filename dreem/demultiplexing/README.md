input:
   - config['demultiplexing']['use']
   - path fasta
   - path fastq1
   - path fastq2
   - library.csv with :
      - barcode_start
      - barcode_stop 

Notes:
 - if use i False, just move the files into the output folder

output:
  - per-sequence fasta/fastq in a folder
