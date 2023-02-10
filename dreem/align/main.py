from dreem.align import align
from dreem.util.reads import FastqUnit


def run(parallel: str,
        max_cpus: int,
        top_dir: str,
        fasta: str,
        phred_enc: int,
        fastqs: tuple[str] = (),
        fastqi: tuple[str] = (),
        fastq1: tuple[str] = (),
        fastq2: tuple[str] = (),
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

    {top_dir}/
        output/
            {sample_1}/
                {ref_1}.bam
                {ref_2}.bam
                ...
            {sample_2}/
            ...
        temp/
            {step_1}/
                {sample_1}/
                    {ref_1}.file
                    {ref_2}.file
                    ...
                {sample_2}/
                ...
            {step_2}/
            ...
        ...

    Parameters
    ----------
    top_dir: str
        Path to the top-level folder into which all final outputs and temporary
        files are written. Must already exist.
    fasta: str
        Path to the reference FASTA file.
    fastqs: tuple[str] †
        Paths to all FASTQ files of single-end reads.
    fastqi: tuple[str] †
        Paths to all FASTQ files of interleaved, paired-end reads.
    fastq1: tuple[str] †
        Paths to all FASTQ files of mate 1 paired-end reads. If given, fastq2
        must also be given, and the sample names (in order) must match.
    fastq2: tuple[str] †
        Paths to all FASTQ files of mate 2 paired-end reads. If given, fastq1
        must also be given, and the sample names (in order) must match.
    fastqs_dir: tuple[str] ‡
        Paths to all directories containing demultiplexed FASTQ files of
        single-end reads.
    fastqi_dir: tuple[str] ‡
        Paths to all directories containing demultiplexed FASTQ files of
        interleaved paired-end reads.
    fastq12_dir: tuple[str] ‡
        Paths to all directories containing demultiplexed FASTQ files of
        mate 1 and mate 2 paired-end reads.
    phred_enc: int
        The ASCII encoding offset of the Phred scores. For example, if
        ```phred_enc = 33```, then a Phred score of 30 will be encoded as the
        ASCII character corresponding to 30 + 33 = 63, which is '?'; and in a
        FASTQ or SAM file, the character 'F' (ASCII value 70) denotes a Phred
        score of 70 - 33 = 37. This encoding of 33 is used by most modern
        Illumina sequencers and is the default.
    **kwargs
        Additional keyword arguments to pass to the alignment function.

    † The file name (minus the extension) will be used as the sample name.
    ‡ The file name (minus the extension) will be used as the reference name,
      and the name of the parent directory will be used as the sample name.
    """

    # Get all FASTQ units that need to be aligned.
    fq_units = list(FastqUnit.from_strs(phred_enc=phred_enc,
                                        fastqs=fastqs,
                                        fastqi=fastqi,
                                        fastq1=fastq1,
                                        fastq2=fastq2,
                                        fastqs_dir=fastqs_dir,
                                        fastqi_dir=fastqi_dir,
                                        fastq12_dir=fastq12_dir))

    # Run the alignment pipeline on each FASTQ.
    return align.run_steps_fqs(top_dir, fasta, fq_units, parallel, max_cpus,
                               **kwargs)
