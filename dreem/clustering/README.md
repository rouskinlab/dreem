
# DREEM Clustering Module
Contributors: Matty Allan, Scott Grote, Yves Martin

## Purpose
- Cluster the bitvector using EM clustering

## Interface

### Input files
- [≥1] `{construct}.orc`. Mutation vector stored in Apache ORC format.
- [=1] `reference.fasta`. Fasta file containing the reference for each sequence of the bitvectors. 
- [≤1] `library.csv`. CSV file containing the following columns:
  - `construct`: name of a sequence in the fasta file, corresponding to a bitvector name.
  - `section_name`: name of a sub-sequence of the construct's sequence, to cluster. If this cell is empty, default value is '"section_start-section_stop"`.
  - `section_start`: 0-index of the start of this sub-sequence w.r.t the global sequence.
  - `section_stop`: 0-index of the end of this sub-sequence w.r.t the global sequence, not included.

### Output files
- [=1] `{prefix}_clustering.json`. 

This json file is structured as follow:
  - `{construct 1}`: name of the bitvector file
    - `{section 1}`: name of the clustered section
      - `{read 1}`: read number
        - `cluster_1`: likelihood that this read belongs to cluster_1
        - `cluster_2`: likelihood that this read belongs to cluster_2
        - ...
      - `{read 2}`: 
        - ...
    - `{section 2}`
      - ...
   - `{construct 2}`
     - ...
        
### Command-line usage

```
dreem-clustering -bv [file] —-fasta [file] —-output [dir]
```

- `dreem-clustering`: wrapper for function run in dreem.clustering.run.
- [≥1] `-bv`: path to bitvectors `{construct}.orc`
- [=1] `--fasta`: path to `reference.fasta` fasta file
- [=1] `--output`: output directory
- [=1] `--prefix`: name of the prefix for the output file prefix (the sample in the context of the dreem module)
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
