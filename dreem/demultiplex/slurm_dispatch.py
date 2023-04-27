import os
import math
import glob
import random

"""
params:
    library csv format
    fastq_dir or fastq1 and fastq2
    outdir
    speed slow or fast or both
    discard or dont discard

"""
B_TO_GB=1073741824

def write_script_as_str(
        fasta="",out_dir="",fq1="",fq2="",
    days="",hours="",memory="",cores="",
    index_tol="",clipped_value="",partition="",
    sample_name="",p_demulti="",library="",
    barcode_start="",barcode_length=""

    ):

    library_opt=f" --library {library} " if library!="" else ""
    barcode_start_opt=f" --barcode-start {barcode_start} " if barcode_start!="" else ""
    barcode_length_opt=f" --barcode-length {barcode_length} " if barcode_length!="" else ""
    os.makedirs(out_dir+"/slurm_logs/",exist_ok=True)

    return f"""
#!/bin/bash 
#SBATCH -c {cores} 
#SBATCH -t 0{days}-{hours}:00
#SBATCH -p {partition}
#SBATCH --mem {memory}G
#SBATCH -o {out_dir}temp/%A.out
#SBATCH -e {out_dir}temp/%A.err


# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/n/app/miniconda3/4.10.3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then 
    eval "$__conda_setup"
else
    if [ -f "/n/app/miniconda3/4.10.3/etc/profile.d/conda.sh" ]; then 
        . "/n/app/miniconda3/4.10.3/etc/profile.d/conda.sh"
    else
        export PATH="/n/app/miniconda3/4.10.3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

conda activate /n/data1/hms/microbiology/rouskin/lab/conda/envs/dreem

dreem --fasta {fasta} \\
    --fastq1 {fq1} \\
    --fastq2 {fq2} \\ 
    --out-dir {out_dir}{sample_name}/out --temp-dir {out_dir}{sample_name}/temp \\ 
    --no-fastqc --clipped {clipped_value} --index-tolerence {index_tol} \\ 
    --parallel_demultiplexing {p_demulti} --demult-on \\ 
    {library_opt}{barcode_start_opt}{barcode_length_opt} 
    """

def slurm_dispatch(list_of_jobs:list,delete_job:bool=True):
    job_folder=f"job_files_random_int_{random.randrange(0,1000)}/"
    job_files=[]
    os.makedirs(job_folder)
    for i in range(len(list_of_jobs)):
        job_file=job_folder+f"job_{i}.sh"
        job_files.append(job_file)
        with open(job_file,"wt") as j_file:
            j_file.write(list_of_jobs[i])

    for j in job_files:
        cmd=f"sbatch {j}"
        print(cmd)
        #os.system("")
    
    if delete_job:
        os.system(f"rm -r {job_folder}")
    



def write_fq1_fq2_script(fasta="",out_dir="",fq1="",fq2="",
    index_tol="",clipped_value="",p_demulti="",barcode_length="",barcode_start="",library=""):

    if(out_dir==""):
            dir_split=fq1.split("/")
            the_dir=""
            for i in range(len(dir_split)-1):
                the_dir+=dir_split[i]

    fasta_len=len(open(fasta,"rt").readlines())/2

    fq1_size=round((os.stat(fq1).st_size)/B_TO_GB)
    sample_name=fq1.split("/")[-1].split("_R1")[0]

    memory_val = fq1_size*2 + (.4*fq1_size)
    
    time = 3+ ((fq1_size)*1.5) +(clipped_value * 5) + (index_tol *5)  +(3*fasta_len/600)
    print(time)

    if(time <= 12.0):
        part="short"
        hours=math.ceil(time)
        days=0
    else: 
        part="medium"
        if time < 23.0:
            hours= time
            days= 0
        elif time < 24.0 and time>=23.0:
            hours=0
            days=1
        else:
            hours=0
            days=math.ceil(time/24)
    

    if (p_demulti):
        cores=15
    else:
        cores=4

    script=write_script_as_str(fasta=fasta,out_dir=out_dir, fq1=fq1,fq2=fq2,days=days,hours=hours,
                               clipped_value=clipped_value,partition=part,
                               sample_name=sample_name,index_tol=index_tol,p_demulti=True,memory=memory_val,
                               barcode_start=barcode_start,barcode_length=barcode_length,library=library)
    
    return script

def dispatch_jobs(fasta,fq1="",fq2="",fastq_dir="",clipped=0,index_tol=0,library="",p_demulti=True,barcode_start="",barcode_length="",out_dir="",delete_jobs=True):

    jobs=[]
    
    if(barcode_start=="" and barcode_length=="" and library==""):
        raise Exception("no barcode locations given")

    if(fq1=="" and fq2=="" and fastq_dir==""):
        raise Exception("no fastqs given")
    
    if(fastq_dir==""):
        
            


        jobs.append(write_fq1_fq2_script(fasta=fasta,out_dir=out_dir,fq1=fq1,fq2=fq2,
                                         index_tol=index_tol,clipped_value=index_tol,
                                         p_demulti=p_demulti,barcode_length=barcode_length,
                                         barcode_start=barcode_start))
    else:
        fq1s=glob.glob(fastq_dir+"*_R1_*")
        fq2s=glob.glob(fastq_dir+"*_R2_*")

        fq1s.sort()
        fq2s.sort()

        for i in range(len(fq1s)):
            jobs.append(write_fq1_fq2_script(fasta=fasta,out_dir=out_dir,fq1=fq1s[i],fq2=fq2s[i],
                                         index_tol=index_tol,clipped_value=index_tol,
                                         p_demulti=p_demulti,barcode_length=barcode_length,
                                         barcode_start=barcode_start))
            
    slurm_dispatch(list_of_jobs=jobs,delete_job=delete_jobs)


dispatch_jobs(fasta="/Users/scottgrote/Documents/absolutely_final_repo/new_ref.fasta",fastq_dir="/Users/scottgrote/Documents/gabe_dreem/dreem_data/test_inputs/",library="fake_library.csv",delete_jobs=False,out_dir="/Users/scottgrote/Documents/absolutely_final_repo/josh_data/")