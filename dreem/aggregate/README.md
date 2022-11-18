# DREEM Aggregate Module
Contributor: Yves Martin

## Purpose
- Turns the bitvector into a csv file.

## Interface

### Input Files
- [≥1] `{construct}.orc`. Mutation vector stored in Apache ORC format.
- [=1] `reference.fasta`. Fasta file containing the reference for each sequence of the bitvector. 
- [≤1] `{sample}_clustering.json`. JSON file containing the clustering likelihood for each read of each bitvector.

### Output Files
- [≥1] ```{sample}.csv```.

### Command-line usage

```dreem-aggregate -bv [file] —fasta [file] —output [file] —per_mp_file [True/False]```

- ```dreem-aggregate```: Wrapper for ```run``` function in ```dreem/aggregate/run.py```. 
- [≥1] `--bit_vector`: `{construct}.orc`
- [=1] `--fasta` : ```{reference}.fasta```
- [≤1] `--clusters`: `{sample}_clustering.json`.
- [=1] `--output`: name of the output file, in this context the sample name.

