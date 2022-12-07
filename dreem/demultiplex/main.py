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
    __demultiplex(f1, f2, library, output_folder, temp_folder)
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
    
    library = pd.read_csv(library)
    f1 = pd.DataFrame(np.loadtxt(f1, dtype=str).reshape(-1, 4), columns=['rname','rseq','+','qual'])
    f2 = pd.DataFrame(np.loadtxt(f2, dtype=str).reshape(-1, 4), columns=['rname','rseq','+','qual'])
    
    for df in [f1,f2]:
        df['construct'] = df['rname'].str.split(':').str[0].str[1:]
        
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
            
            with open(output_folder+'/'+construct+'_R'+str(primer)+'.fastq', 'w') as f:
                for _, r in df.iterrows():
                    if barcode_in_read(barcode, r['rseq']):
                        f.write(r['rname']+'\n'+r['rseq']+'\n+\n'+r['qual']+'\n')
    return 1

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

    Returns
    -------
    1 if successful, 0 otherwise.

    """
    # Get the paths
    root = args['out_dir']
    temp_folder = os.path.join(root,'temp','demultiplexing')
    output_folder = os.path.join(root,'output','demultiplexing')
    fastq1 = args['fastq1'] if type(args['fastq1']) == list else [args['fastq1']]
    fastq2 = args['fastq2'] if type(args['fastq2']) == list else [args['fastq2']]

    # Make the folders
    util.make_folder(output_folder)
    util.make_folder(temp_folder)

    # Demultiplex
    for f1 in fastq1:
        for f2 in fastq2:
            if f1[:-len('_R1.fastq')] == f2[:-len('_R2.fastq')]:
                sample = f1.split('/')[-1][:-len('_R1.fastq')]
                util.make_folder(os.path.join(output_folder, sample))
                assert demultiplex(f1, f2, args['library'], os.path.join(output_folder, sample), os.path.join(temp_folder, sample)), "Demultiplexing failed"
    return 1