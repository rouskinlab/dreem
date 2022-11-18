
# DREEM Clustering Module
Contributors: Matty Allan, Scott Grote, Yves Martin

## Purpose
- Cluster the bitvector using EM clustering

## Interface

### Input files
- [≥1] `{bit_vector}.orc`. Mutation vector stored in Apache ORC format.
- [≥1] `{sample}.fasta`. Fasta file containing the reference for each sequence of the bitvector. 

### Output files
- [≥1] `{bit_vector}_{cluster_number}.orc`

### Command-line usage

```
dreem-clustering -bv [file] —-fasta [file] —-output [dir]
```

- `dreem-clustering`: wrapper for function run in dreem.clustering.run.
- [≥1] `-bv`: path to bitvector
- [≥1] `--fasta`: path to fasta file
- [=1] `--output`: output directory
- [≤1] `--N_clusters`: number of clusters
- [≤1] `--max_N_clusters`: use the optimal number of clusters below or equal to this value
- [≤1] `--signal_thresh`: signal threshold #TODO, float in [0,1]
- [≤1] `--include_G_U`: include G and U bases 
- [≤1] `--include_del`: include deleted bases
- [≤1] `--info_thresh`: #TODO
- [≤1] `--min_reads`: minimum amount of reads for a sequence to be clustered
- [≤1] `--convergence_cutoff`: float #TODO
- [≤1] `--num_runs`: number of iterations with random parameters initialisation
- [≥0] `-c`: 3 arguments: reference sequence name, start, stop. Multiple `-c` arguments are possible.
- [≥0] `-p`: 3 arguments: reference sequence name, forward primer, reverse primer. Multiple `-p` arguments are possible.
