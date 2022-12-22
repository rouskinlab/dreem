

# Samples.csv

The samples file is a csv file containing per-sample information. The columns are:
- sample: the sample name, corresponding to the prefix of the fastq file (ex: sample: 'sample1' corresponds to the input file 'sample1_R1.fastq')
- user: the user name
- date: the date of the experiment
- exp_env: the experimental environment, in_vivo or in_vitro
- temperature_k: the temperature in Kelvin
- inc_time_tot_secs: the total incubation time in seconds
- DMS_conc_mM: the concentration of DMS in mM
- buffer (if exp_env is in_vitro): the exact buffer including Mg, eg 300mM Sodium Cacodylate , 3mM Mg
- cell_line (if exp_env is in_vivo): the cell line
- [any other columns]: any other columns will be added to corresponding sample in the output file
