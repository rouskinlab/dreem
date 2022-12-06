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



def __demultiplex(f1, f2, library, output_folder, temp_folder):
    """Demultiplex a pair of FASTQ files."""
    
    library = pd.read_csv(library)
    f1 = pd.read_csv(np.loadtxt('data.txt').reshape(3, -1).T, header=None, names=['construct', 'read'])

    return 1