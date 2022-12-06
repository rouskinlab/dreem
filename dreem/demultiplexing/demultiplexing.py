import dreem.util as util
import os  
import pandas as pd
import numpy as np

def demultiplex(f1, f2, library, output_folder, temp_folder):
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
    temp_folder: str
        Use as a temporary folder. 

    returns
    -------
    1 if successful, 0 otherwise.
    """

    try:
        return __demultiplex(f1, f2, library, output_folder, temp_folder)
    except Exception as e:
        print(e)
        return 0

def hamming_distance(s1, s2):
    assert len(s1) == len(s2)
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def reverse_complement(seq):
    return seq[::-1].translate(str.maketrans('ATCG', 'TAGC'))

def barcode_in_read(barcode, read, max_hamming_distance=1):
    if len(read) < len(barcode):
        return False
    for i in range(len(read)-len(barcode)+1):
        if hamming_distance(barcode, read[i:i+len(barcode)]) <= max_hamming_distance:
            return True
    return False

def __demultiplex(f1, f2, library, output_folder, temp_folder):
    """Demultiplex a pair of FASTQ files."""

    constructs = library['construct'].unique()
    barcodes = library['barcode'].unique()

    assert len(constructs) == len(barcodes)
    for construct, barcode in zip(constructs, barcodes):
        # chec if the barcode andthe construct are on the same row
        assert library.loc[library['construct']==construct, 'barcode'].values[0] == barcode

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # copy the reads from the fastq files that contain the barcode in a fastq file named after the construct in the output folder
    for primer in [1,2]:
        for construct, barcode in zip(constructs, barcodes):
            if primer == 2:
                barcode = reverse_complement(barcode)
            
            df = f1 if primer == 1 else f2
            
            with open(output_folder+construct+'_R'+str(primer)+'.fastq', 'w') as f:
                for _, r in df.iterrows():
                    if barcode_in_read(barcode, r['rseq']):
                        f.write(r['rname']+'\n'+r['rseq']+'\n+\n'+r['qual']+'\n')
    return 1