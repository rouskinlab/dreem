from click import command, pass_obj

from .align import run_steps_fqs
from .reads import FastqUnit
from ..util.cli import (DreemCommandName, dreem_command,
                        opt_fasta,
                        opt_fastqs, opt_fastqi, opt_fastq1, opt_fastq2,
                        opt_fastqs_dir, opt_fastqi_dir, opt_fastq12_dir,
                        opt_phred_enc,
                        opt_out_dir, opt_temp_dir,
                        opt_rerun, opt_resume, opt_save_temp,
                        opt_parallel, opt_max_procs,
                        opt_fastqc, opt_fastqc_extract,
                        opt_cutadapt,
                        opt_cut_a1, opt_cut_g1, opt_cut_a2, opt_cut_g2,
                        opt_cut_o, opt_cut_e, opt_cut_q1, opt_cut_q2,
                        opt_cut_m, opt_cut_indels,
                        opt_cut_discard_trimmed, opt_cut_discard_untrimmed,
                        opt_cut_nextseq,
                        opt_bt2_local, opt_bt2_unal,
                        opt_bt2_discordant, opt_bt2_mixed,
                        opt_bt2_dovetail, opt_bt2_contain,
                        opt_bt2_i, opt_bt2_x, opt_bt2_score_min,
                        opt_bt2_s, opt_bt2_l, opt_bt2_d, opt_bt2_r,
                        opt_bt2_gbar, opt_bt2_dpad, opt_bt2_orient,
                        opt_rem_buffer)
from ..util import docdef


# Parameters for command line interface
params = [
    # Input files
    opt_fasta,
    opt_fastqs,
    opt_fastqi,
    opt_fastq1,
    opt_fastq2,
    opt_fastqs_dir,
    opt_fastqi_dir,
    opt_fastq12_dir,
    # FASTQ options
    opt_phred_enc,
    # Output directories
    opt_out_dir,
    opt_temp_dir,
    # File generation
    opt_rerun,
    opt_resume,
    opt_save_temp,
    # Parallelization
    opt_parallel,
    opt_max_procs,
    # FASTQC options
    opt_fastqc,
    opt_fastqc_extract,
    # Cutadapt options
    opt_cutadapt,
    opt_cut_a1,
    opt_cut_g1,
    opt_cut_a2,
    opt_cut_g2,
    opt_cut_o,
    opt_cut_e,
    opt_cut_q1,
    opt_cut_q2,
    opt_cut_m,
    opt_cut_indels,
    opt_cut_discard_trimmed,
    opt_cut_discard_untrimmed,
    opt_cut_nextseq,
    # Bowtie2 options
    opt_bt2_local,
    opt_bt2_discordant,
    opt_bt2_mixed,
    opt_bt2_dovetail,
    opt_bt2_contain,
    opt_bt2_unal,
    opt_bt2_i,
    opt_bt2_x,
    opt_bt2_score_min,
    opt_bt2_s,
    opt_bt2_l,
    opt_bt2_gbar,
    opt_bt2_d,
    opt_bt2_r,
    opt_bt2_dpad,
    opt_bt2_orient,
    opt_rem_buffer
]


@command(DreemCommandName.ALIGN.value, params=params)
# Pass context object.
@pass_obj
# Turn into DREEM command.
@dreem_command(imports=("fastqs_dir", "fastqi_dir", "fastq12_dir"),
               exports=("fasta", "phred_enc"),
               result_key="bamf")
def cli(**kwargs):
    return run(**kwargs)


@docdef.auto()
def run(*,
        fastqs: tuple[str],
        fastqi: tuple[str],
        fastq1: tuple[str],
        fastq2: tuple[str],
        fastqs_dir: tuple[str],
        fastqi_dir: tuple[str],
        fastq12_dir: tuple[str],
        phred_enc: int,
        **kwargs):
    """
    Run the alignment module.

    Align the reads to the set of reference sequences and output one BAM file
    for each sample aligned to each reference in the directory 'output'.
    Temporary intermediary files are written in the directory 'temp' and then
    deleted after they are no longer needed.

    {output_dir}/
        {sample_1}/
            {ref_1}.bam
            {ref_2}.bam
            ...
        {sample_2}/
        ...

    {temp_dir}/
        {step_1}/
            {sample_1}/
                {ref_1}.file
                {ref_2}.file
                ...
            {sample_2}/
            ...
        {step_2}/
        ...
    """
    
    # FASTQ files of read sequences may come from up to seven different
    # sources (i.e. each argument beginning with "fq_unit"). This step
    # collects all of them into one list (fq_units) and also bundles
    # together pairs of FASTQ files containing mate 1 and mate 2 reads.
    fq_units = list(FastqUnit.from_strs(fastqs=fastqs,
                                        fastqi=fastqi,
                                        fastq1=fastq1,
                                        fastq2=fastq2,
                                        fastqs_dir=fastqs_dir,
                                        fastqi_dir=fastqi_dir,
                                        fastq12_dir=fastq12_dir,
                                        phred_enc=phred_enc,
                                        no_dup_samples=True))

    # Run the alignment pipeline on every FASTQ.
    return run_steps_fqs(fq_units, **kwargs)
