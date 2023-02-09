from collections import Counter
import pandas as pd
import numpy as np
from scipy import signal
import datetime
from dreem.util.files_sanity import check_library
from dreem.util import path
from dreem.util.reads import FastqUnit

def demultiplex(fq_unit: FastqUnit,
                fasta: path.RefsetSeqInFilePath,
                out_dir: path.ModuleDirPath,
                library: str,
                max_barcode_mismatches: int):
    """Demultiplex a pair of FASTQ files.

    Publishes to `output_folder` a pair of FASTQ files for each construct, named {construct}_R1.fastq and {construct}_R2.fastq.

    Parameters
    ----------
    fq_unit: FastqUnit
        One single-end or interleaved paired-end FASTQ file, or paired-end
        reads with mates 1 and 2 in two separate FASTQ files.
    fasta: RefsetSeqInFilePath
        FASTA file containing the reference sequences.
    library: pd.DataFrame
        Columns are (non-exclusively): ['construct', 'barcode_start', 'barcode']
    out_dir: str
        Where to output the results.
    max_barcode_mismatches: int
        Maximum number of mutations allowed on the barcode.
        
    Returns
    -------
    dict[str, FastqUnit]
        Dictionary mapping construct names to the demultiplexed FASTQ files.
    """

    out_dir.path.mkdir(parents=True, exist_ok=True)

    # Remove the report file if it exists
    report_path = out_dir.path.joinpath('report.txt')
    if report_path.is_file():
        report_path.unlink()
    
    library = check_library(pd.read_csv(library), str(fasta.path))
    constructs = library['construct'].unique()
    barcodes = library['barcode'].unique()

    assert len(constructs) == len(barcodes)
    for construct, barcode in zip(constructs, barcodes):
        # check if the barcode and the construct are on the same row
        assert library.loc[library['construct']==construct, 'barcode'].values[0] == barcode

    lost_reads = fq_unit.trans(path.ReadsInToReadsOut, ref="lost_reads",
                               **out_dir.dict())

    construct_fastqs: dict[str, FastqUnit] = dict()

    # copy the reads from the fastq files that contain the barcode in a fastq file named after the construct in the output folder
    for second, fq in enumerate(fq_unit.inputs):
        
        # infos for the report
        perfect_matches_count = 0
        off_matches_count = {k:0 for k in range(1, max_barcode_mismatches+1)}
        lost_reads_count = 0
        count_per_construct = {construct:0 for construct in constructs}
        barcode_shifts = []
                        
        with open(fq.path, 'r') as f:
            while True:
                header, sequence, quality = read_fastq_line(f)
                if not header:
                    break

                for construct, barcode in zip(constructs, barcodes):
                    
                    minimal_corr_score = worst_matching_score(barcode, max_barcode_mismatches)
                    
                    if second:
                        barcode = reverse_complement(barcode)
                    
                    corr = compute_correlation(embed_sequence_as_binary(barcode), embed_sequence_as_binary(sequence))
                    is_match = barcode_in_read(corr, minimal_corr_score)
                    if is_match:
                        if max(corr) == 1:
                            perfect_matches_count += 1
                        else:
                            for k in range(1, max_barcode_mismatches+1):
                                if max(corr) >= worst_matching_score(barcode, k):
                                    off_matches_count[k] += 1
                                    break
                            
                        barcode_start = library.loc[library['construct']==construct, 'barcode_start'].values[0]
                        if second:
                            barcode_shifts.append(len(sequence) - np.argmax(corr) - barcode_start - len(barcode))
                        else:
                            barcode_shifts.append(np.argmax(corr) - barcode_start)
                            
                        count_per_construct[construct] += 1

                        # TODO: Open each construct's file a minimum number of times (not once per line)
                        try:
                            with open(construct_fastqs[construct].paths[second], "a") as g:
                                write_fastq_line(g, header, sequence, quality)
                        except KeyError:
                            construct_fq_unit = fq_unit.trans(path.ReadsInToReadsOut,
                                                              ref=construct,
                                                              **out_dir.dict())
                            construct_fastqs[construct] = construct_fq_unit
                            with open(construct_fq_unit.paths[second], "w") as g:
                                write_fastq_line(g, header, sequence, quality)
                        break
                else:
                    lost_reads_count += 1
                    with open(lost_reads.paths[second], 'a') as g:
                        write_fastq_line(g, header, sequence, quality)
        write_report(fq.path, report_path, perfect_matches_count, off_matches_count, lost_reads_count, barcode_shifts, count_per_construct)
    return construct_fastqs

def write_report(fastq, report_path, perfect_matches_count, off_matches_count, lost_reads_count, barcode_shifts, count_per_construct):
    """Write a report of the demultiplexing process for the given fastq file."""
    with open(report_path, 'a') as f:
        f.write("Time: " + str(datetime.datetime.now()) + "\n")
        f.write('\n'+'='*len('Demultiplexing report for ' + fastq) + '\n')
        f.write('Demultiplexing report for ' + fastq + '\n')
        f.write('='*len('Demultiplexing report for ' + fastq) + '\n')
        f.write('Count of perfect matches: ' + str(perfect_matches_count) + '\n')
        for k in off_matches_count:
            f.write('Count of ' + str(k) + '-off matches: ' + str(off_matches_count[k]) + '\n')
        f.write('Count of lost reads: ' + str(lost_reads_count) + '\n')
        f.write('Count of reads per barcode position: ' + str(bin_positions(barcode_shifts)) + '\n')
        f.write('\nCount of reads per construct: ' + '\n' + '-'*len('Count of reads per construct:') + '\n')
        for construct in count_per_construct:
            f.write(construct + ': ' + str(count_per_construct[construct]) + '\n')
        f.write('='*len('Demultiplexing report for ' + fastq) + '\n')
        
def worst_matching_score(barcode, max_muts=1):
    return 1. - float(max_muts)/len(barcode) -1E-9
                                    
def bin_positions(positions):
    """Turns a list of positions into a dictionary of bins."""
    bins = {}
    for pos in positions:
        if pos not in bins:
            bins[pos] = 0
        bins[pos] += 1
    bins = {k: bins[k] for k in sorted(bins)}
    return bins

def read_fastq_line(f):
    """Read the line of a fastq file and return a tuple (header, sequence, quality)"""
    header = f.readline().strip()
    sequence = f.readline().strip()
    f.readline()
    quality = f.readline().strip()
    return header, sequence, quality

def write_fastq_line(f, header, sequence, quality):
    """Write a fastq file line."""
    f.write(header + '\n')
    f.write(sequence + '\n')
    f.write('+\n')
    f.write(quality + '\n')

def embed_sequence_as_binary(sequence):
    """Each sequence is represented as 4 binary vector of length len(sequence), one per base A C T G."""
    return np.array([np.array([base == 'A', base == 'C', base == 'T', base == 'G']) for base in sequence], dtype=np.int8).T

def compute_correlation(barcode, read):
    """ Use the correlation between the barcode and the read to determine if the barcode is in the read."""
    return signal.correlate(read,barcode, mode="valid", method="auto").squeeze()/barcode.shape[1]

def barcode_in_read(corr, min_corr_score):
    """Return True if the correlation score is above the threshold, False otherwise."""
    return max(corr) > min_corr_score

def hamming_distance(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def reverse_complement(seq):
    return seq[::-1].translate(str.maketrans('ATCG', 'TAGC'))

def next_base(base):
    return {'A':'T','T':'C','C':'G','G':'A',0:1}[base]

def run(top_dir: str, fasta: str, phred_enc: int,
        fastqs: tuple[str], fastqi: tuple[str],
        fastq1: tuple[str], fastq2: tuple[str],
        library: str, max_barcode_mismatches: int):
    """Run the demultiplexing pipeline.

    Demultiplexes the reads and outputs one fastq file per construct in the directory `output_path`, using `temp_path` as a temp directory.

    Parameters from args:
    -----------------------
    fastq1: str
        Path to the FASTQ file or list of paths to the FASTQ files, forward primer.
    fastq2: str
        Path to the FASTQ file or list of paths to the FASTQ files, reverse primer.
    fasta: str
        Path to the FASTA file.
    interleaved: bool
        If True, the FASTQ files are interleaved.
    library: str
        Path to the library file. Columns are (non-excusively): ['construct', 'barcode_start', 'barcode']
    out_dir: str
        Name of the output directory.
    max_barcode_mismatches: int
        Maximum number of mutations allowed on the barcode.
    verbose: bool
        Print progress to stdout (default: no).
        
    Returns
    -------
    dict[str, dict[str, FastqUnit]]
        Dictionary mapping sample names to

    """

    fasta_path = path.RefsetSeqInFilePath.parse_path(fasta)
    # Make the folders
    out_dir = path.ModuleDirPath(top=top_dir,
                                 partition=path.Partition.OUTPUT,
                                 module=path.Module.DEMULT)

    # Demultiplex.
    demultiplexed: dict[str, dict[str, FastqUnit]] = dict()
    fq_units = FastqUnit.from_strs(fastqs=fastqs, fastqi=fastqi,
                                   fastq1=fastq1, fastq2=fastq2,
                                   phred_enc=phred_enc, demult=True)
    # Ensure that no sample names are duplicated.
    if dups := [sample for sample, count
                in Counter(fq.sample for fq in fq_units).items()
                if count > 1]:
        raise ValueError(f"Got duplicate sample names: {', '.join(dups)}")

    # TODO: Parallelize with multiprocessing.Pool.starmap
    for fq_unit in fq_units:
        demultiplexed[fq_unit.sample] = demultiplex(
            fq_unit=fq_unit,
            fasta=fasta_path,
            out_dir=out_dir,
            library=library,
            max_barcode_mismatches=max_barcode_mismatches)
    return demultiplexed
