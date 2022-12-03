# DREEM Aggregate Module
Contributor: Yves Martin

## Purpose
Counts the bitvector: 
- mutation per residue across reads
- number of mutations per read
- coverage per residue across reads
- insertions per residue across reads
- deletions per residue across reads
- ...

## Interface

### Input Files
- [≥1] `/{sample}`. A directory containing mutation vector stored in Apache ORC format.
```
/{sample}
  - {construct_1}.orc
  - {construct_2}.orc
  - ...
```
- [=1] `reference.fasta`. Fasta file containing the reference for each sequence of the bitvector. 
- [≤1] `clustering.json`. JSON file containing the clustering likelihood for each read of each bitvector.

### Output Files
- [≥1] ```{out_dir}/output/aggregate/{sample}.json```.

### Command-line usage

```dreem-aggregate -bv [file] —fasta [file]  --out_dir [sample] —per_mp_file [True/False]```

- ```dreem-aggregate```: Wrapper for ```run``` function in ```dreem/aggregate/run.py```. 
- [≥1] `--sample`: `/{sample}`
- [=1] `--fasta` : ```{reference}.fasta```
- [≤1] `--clusters`: `clustering.json`.
- [=1] `--out_dir`: name of the output directory.

