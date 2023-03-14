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
    # Inputs
    opt_fasta,
    opt_fastqs,
    opt_fastqi,
    opt_fastq1,
    opt_fastq2,
    opt_fastqs_dir,
    opt_fastqi_dir,
    opt_fastq12_dir,
    opt_phred_enc,
    # Outputs
    opt_out_dir,
    opt_temp_dir,
    opt_rerun,
    opt_resume,
    opt_save_temp,
    # Parallelization
    opt_parallel,
    opt_max_procs,
    # FASTQC
    opt_fastqc,
    opt_fastqc_extract,
    # Cutadapt
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
    # Bowtie2
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
    # Post-processing
    opt_rem_buffer
]


@command(DreemCommandName.ALIGN.value, params=params)
# Pass context object.
@pass_obj
# Turn into DREEM command.
@dreem_command(imports=("fasta", "fastqs_dir", "fastqi_dir", "fastq12_dir"),
               exports=("fasta", "phred_enc"),
               result_key="bamf")
def cli(**kwargs):
    return run(**kwargs)


@docdef.auto()
def run(*,
        # Inputs
        fasta: str,
        fastqs: tuple[str],
        fastqi: tuple[str],
        fastq1: tuple[str],
        fastq2: tuple[str],
        fastqs_dir: tuple[str],
        fastqi_dir: tuple[str],
        fastq12_dir: tuple[str],
        phred_enc: int,
        # Outputs
        out_dir: str,
        temp_dir: str,
        save_temp: bool,
        rerun: bool,
        resume: bool,
        # Parallelization
        max_procs: int,
        parallel: bool,
        # FASTQC
        fastqc: bool,
        fastqc_extract: bool,
        # Cutadapt
        cut: bool,
        cut_q1: int,
        cut_q2: int,
        cut_g1: str,
        cut_a1: str,
        cut_g2: str,
        cut_a2: str,
        cut_o: int,
        cut_e: float,
        cut_indels: bool,
        cut_nextseq: bool,
        cut_discard_trimmed: bool,
        cut_discard_untrimmed: bool,
        cut_m: int,
        # Bowtie2
        bt2_local: bool,
        bt2_discordant: bool,
        bt2_mixed: bool,
        bt2_dovetail: bool,
        bt2_contain: bool,
        bt2_unal: bool,
        bt2_score_min: str,
        bt2_i: int,
        bt2_x: int,
        bt2_gbar: int,
        bt2_l: int,
        bt2_s: str,
        bt2_d: int,
        bt2_r: int,
        bt2_dpad: int,
        bt2_orient: str,
        # Post-processing
        rem_buffer: int):
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
                                        phred_enc=phred_enc))

    # Run the alignment pipeline on every FASTQ.
    return run_steps_fqs(fq_units=fq_units,
                         fasta=fasta,
                         out_dir=out_dir,
                         temp_dir=temp_dir,
                         save_temp=save_temp,
                         rerun=rerun,
                         resume=resume,
                         max_procs=max_procs,
                         parallel=parallel,
                         fastqc=fastqc,
                         fastqc_extract=fastqc_extract,
                         cut=cut,
                         cut_q1=cut_q1,
                         cut_q2=cut_q2,
                         cut_g1=cut_g1,
                         cut_a1=cut_a1,
                         cut_g2=cut_g2,
                         cut_a2=cut_a2,
                         cut_o=cut_o,
                         cut_e=cut_e,
                         cut_indels=cut_indels,
                         cut_nextseq=cut_nextseq,
                         cut_discard_trimmed=cut_discard_trimmed,
                         cut_discard_untrimmed=cut_discard_untrimmed,
                         cut_m=cut_m,
                         bt2_local=bt2_local,
                         bt2_discordant=bt2_discordant,
                         bt2_mixed=bt2_mixed,
                         bt2_dovetail=bt2_dovetail,
                         bt2_contain=bt2_contain,
                         bt2_unal=bt2_unal,
                         bt2_score_min=bt2_score_min,
                         bt2_i=bt2_i,
                         bt2_x=bt2_x,
                         bt2_gbar=bt2_gbar,
                         bt2_l=bt2_l,
                         bt2_s=bt2_s,
                         bt2_d=bt2_d,
                         bt2_r=bt2_r,
                         bt2_dpad=bt2_dpad,
                         bt2_orient=bt2_orient,
                         rem_buffer=rem_buffer)
