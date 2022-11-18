
## Command line instruction


```
dreem-clustering -bv [file] —-fasta [file] —-output [dir]
```

### Mandatory arguments:

`-bv`: path to bitvector

`--fasta`: path to fasta file

`--output`: output directory

### Optional arguments

`--N_clusters`: number of clusters

`--max_N_clusters`: use the optimal number of clusters below or equal to this value

`--signal_thresh`: signal threshold #TODO, float in [0,1]

`--include_G_U`: include G and U bases 

`--include_del`: include deleted bases

`--info_thresh`: #TODO

`--min_reads`: minimum amount of reads for a sequence to be clustered

`--convergence_cutoff`: float #TODO

`--num_runs`: number of iterations with random parameters initialisation

`-c`: 3 arguments: reference sequence name, start, stop. Multiple `-c` arguments are possible.

`-p`: 3 arguments: reference sequence name, forward primer, reverse primer. Multiple `-p` arguments are possible.
