from dreem.align.align import run_steps_fqs
from dreem.align.reads import FastqUnit
from dreem.util.logio import set_verbosity


def run(fasta: str,
        phred_enc: int,
        fastqs: tuple[str],
        fastqi: tuple[str],
        fastq1: tuple[str],
        fastq2: tuple[str],
        fastqs_dir: tuple[str] = (),
        fastqi_dir: tuple[str] = (),
        fastq12_dir: tuple[str] = (),
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

    Parameters
    ----------
    fasta: str,
        FASTA file containing all reference sequences
    fastqs: tuple[str] †
        FASTQ files of single-end reads
    fastqi: tuple[str] †
        FASTQ files of interleaved, paired-end reads
    fastq1: tuple[str] †
        FASTQ files of mate 1 paired-end reads. If given, fastq2 must
        also be given, and the sample names (in order) must match
    fastq2: tuple[str] †
        FASTQ files of mate 2 paired-end reads. If given, fastq1 must
        also be given, and the sample names (in order) must match
    fastqs_dir: tuple[str] ‡
        Directories of demultiplexed FASTQ files of single-end reads
    fastqi_dir: tuple[str] ‡
        Directories of demultiplexed FASTQ files of interleaved,
        paired-end reads
    fastq12_dir: tuple[str] ‡
        Directories of demultiplexed FASTQ files of mate 1 and mate 2
        paired-end reads
    phred_enc: int
        The ASCII encoding offset of the Phred scores. For example, if
        ```phred_enc = 33```, then a Phred score of 30 will be encoded
        as the ASCII character corresponding to 30 + 33 = 63, which is
        '?'; and in a FASTQ or SAM file, the character 'F' (ASCII value
        70) denotes a Phred score of 70 - 33 = 37. This encoding of 33
        is used by most modern Illumina sequencers and is the default.
    **kwargs
        Additional keyword arguments to pass to the alignment function

    † The file name (minus its extension) becomes the sample name.
    ‡ The file name (minus its extension) becomes the reference name;
      the name of the parent directory becomes the sample name.
    """

    # FASTQ files of read sequences may come from up to seven different
    # sources (i.e. each argument beginning with "fastq"). This step
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
    return fasta, run_steps_fqs(fasta=fasta, fq_units=fq_units, **kwargs)
