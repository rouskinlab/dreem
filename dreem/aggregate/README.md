# DREEM Alignment Module
Contributor: Yves Martin

## Purpose
- Turns the bitvector into a csv file.

## Interface

### Input Files
- [≥1] `{bit_vector}.orc`. Mutation vector stored in Apache ORC format.
- [≥1] `{bit_vector}.orc`. Mutation vector stored in Apache ORC format.
- [=1] `{reference}.fasta`. Fasta file containing the reference for each sequence of the bitvector. 

### Output Files
- [≥1] ```{sample}.bam```.

### Command-line usage

```dreem-aggregate -bv [file] —fasta [file] —output [file] —per_mp_file [True/False]```

- ```dreem-aggregate```: Wrapper for ```run``` function in ```dreem/aggregate/run.py```. 
- `--bit_vector`: `{bit_vector}.orc`
- `--fasta` : ```{reference}.fasta```
- `--output`: name of the output file(s)
