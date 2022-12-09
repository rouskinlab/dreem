import dreem.util as util
import os  
import pandas as pd
import numpy as np
from scipy import signal
import datetime

def demultiplex(f1, f2, library, output_folder, max_mutations_on_barcode=1):
    """Demultiplex a pair of FASTQ files.

    Publishes to `output_folder` a pair of FASTQ files for each construct, named {construct}_R1.fastq and {construct}_R2.fastq.

    Parameters
    ----------
    f1: str
        Path to the FASTQ file, forward primer.
    f2: str
        Path to the FASTQ file, reverse primer.
    library: pd.DataFrame
        Columns are (non-exclusively): ['construct', 'barcode_start', 'barcode']
    output_folder: str
        Where to output the results.
    max_mutations_on_barcode: int
        Maximum number of mutations allowed on the barcode.

    returns
    -------
    1 if successful, 0 otherwise.
    """
    
    library = pd.read_csv(library)
    constructs = library['construct'].unique()
    barcodes = library['barcode'].unique()

    assert len(constructs) == len(barcodes)
    for construct, barcode in zip(constructs, barcodes):
        # check if the barcode and the construct are on the same row
        assert library.loc[library['construct']==construct, 'barcode'].values[0] == barcode

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)        

    # copy the reads from the fastq files that contain the barcode in a fastq file named after the construct in the output folder
    for primer in [1,2]:
        
        fq = f1 if primer == 1 else f2
        
        # infos for the report
        perfect_matches_count = 0
        off_matches_count = {k:0 for k in range(1, max_mutations_on_barcode+1)}
        lost_reads_count = 0
        count_per_construct = {construct:0 for construct in constructs}
        barcode_shifts = []  
        
        construct_fastq_was_written = {construct:False for construct in constructs}
                        
        with open(fq, 'r') as f:
            while True:
                header, sequence, quality = read_fastq_line(f)
                if not header:
                    break
                
                flag_match = False
                for construct, barcode in zip(constructs, barcodes):
                    
                    minimal_corr_score = worst_matching_score(barcode, max_mutations_on_barcode)
                    
                    if primer == 2:
                        barcode = reverse_complement(barcode)
                    
                    corr = compute_correlation(embed_sequence_as_binary(barcode), embed_sequence_as_binary(sequence))
                    is_match = barcode_in_read(corr, minimal_corr_score)
                    if is_match:
                        if max(corr) == 1:
                            perfect_matches_count += 1
                        else:
                            for k in range(1, max_mutations_on_barcode+1):
                                if max(corr) >= worst_matching_score(barcode, k):
                                    off_matches_count[k] += 1
                                    break
                            
                        barcode_start = library.loc[library['construct']==construct, 'barcode_start'].values[0]
                        if primer == 2:
                            barcode_shifts.append(len(sequence) - np.argmax(corr) - barcode_start - len(barcode))
                        else:
                            barcode_shifts.append(np.argmax(corr) - barcode_start)
                            
                        count_per_construct[construct] += 1
                        
                        if not construct_fastq_was_written[construct]:
                            with open(os.path.join(output_folder, construct + '_R' + str(primer) + '.fastq'), 'w') as g:
                                write_fastq_line(g, header, sequence, quality)
                            construct_fastq_was_written[construct] = True
                        else:
                            with open(os.path.join(output_folder, construct + '_R' + str(primer) + '.fastq'), 'a') as g:
                                write_fastq_line(g, header, sequence, quality)
                            
                        flag_match = True
                        break
                    
                if not flag_match:
                    lost_reads_count += 1
                    with open(os.path.join(output_folder, 'lost_reads_R' + str(primer) + '.fastq'), 'a') as g:
                        write_fastq_line(g, header, sequence, quality)
                        
        write_report(fq, output_folder, perfect_matches_count, off_matches_count, lost_reads_count, barcode_shifts, count_per_construct)
    return 1

def write_report(fastq, output_folder, perfect_matches_count, off_matches_count, lost_reads_count, barcode_shifts, count_per_construct):
    """Write a report of the demultiplexing process for the given fastq file."""
    with open(os.path.join(output_folder, 'report.txt'), 'a') as f:
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

def run(**args):
    """Run the demultiplexing pipeline.

    Demultiplexes the reads and outputs one fastq file per construct in the directory `output_path`, using `temp_path` as a temp directory.

    Parameters from args:
    -----------------------
    library: str
        Path to the library file. Columns are (non-excusively): ['construct', 'barcode_start', 'barcode']
    fastq1: str
        Path to the FASTQ file or list of paths to the FASTQ files, forward primer.
    fastq2: str
        Path to the FASTQ file or list of paths to the FASTQ files, reverse primer.
    out_dir: str
        Name of the output directory.
    coords: tuple
        coordinates for reference: '-c ref-name first last'
    primers: tuple
        primers for reference: '-p ref-name fwd rev'
    fill: bool
        Fill in coordinates of reference sequences for which neither coordinates nor primers were given (default: no).
        
    Returns
    -------
    1 if successful, 0 otherwise.

    """

    assert args['library'] is not None, "Library file does not exist"

    
    # Get the paths
    output_folder = args['out_dir']
    fastq1 = args['fastq1'] if type(args['fastq1']) == list else [args['fastq1']]
    fastq2 = args['fastq2'] if type(args['fastq2']) == list else [args['fastq2']]

    # Make the folders
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Remove the report file if it exists
    report_path = os.path.join(output_folder, 'report.txt')
    if os.path.exists(report_path):
        os.remove(report_path)
    
    # Demultiplex
    for f1, f2 in zip(fastq1, fastq2):
        assert demultiplex(f1, f2, args['library'], output_folder), "Demultiplexing failed"
    return 1