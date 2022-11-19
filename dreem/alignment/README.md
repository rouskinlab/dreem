# DREEM Alignment Module
Contributor: Matty Allan

## Purpose
- Convert sequencing files into BAM files

## Interface

### Input Files
- [=1] ```my_fasta.fasta```. Sequence record file that contains the names and sequences of all references.
- [=1] ```my_fastq_R1.fastq```. Sequence alignment file(s) containing one or several sequences. 
- [=1] ```my_fastq_R2.fastq```. Sequence alignment file(s) containing one or several sequences. 

### Output Files
- [≥1] ```{construct}.bam```. One BAM file per construct (sequence name) in `my_fastq_R1.fastq` and  `my_fastq_R2.fastq`.  

### Command-line usage

```dreem-alignment --fastq1 [path to file] --fastq2 [path to file] --fasta [path to file]```

- ```dreem-alignment```: Wrapper for ```run``` function in ```dreem/alignment/run.py```. 
- [=1] `--fasta` : ```my_fasta.fasta```
- [=1] `--fastq1`: ```my_fastq_R1.fastq```
- [=1] `--fastq2`: ```my_fastq_R2.fastq```
- [≥1] `output`: output repository.
- `[bowtie2 args]`
- `[cutadapt args]`
