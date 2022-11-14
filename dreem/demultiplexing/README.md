input:
   - path fasta
   - path fastq1
   - path fastq2
   - library.csv with :
      - barcode_start
      - barcode_stop 

output:
  - per-sequence fasta/fastq in a folder
