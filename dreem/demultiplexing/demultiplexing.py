import dreem.util as util
import os  

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
        Columns are (non-exclusively): ['construct', 'barcode_start', 'barcode_end', 'barcode_sequence']
    output_folder: str
        Where to output the results.
    temp_folder: str
        Use as a temporary folder. 

    returns
    -------
    1 if successful, 0 otherwise.
    """


    # As a placeholder, just copy the files to the output folder
    util.run_cmd(f"cp {f1} {os.path.join(output_folder,'mttr-6-alt-h3_R1.fastq')}")
    util.run_cmd(f"cp {f2} {os.path.join(output_folder,'mttr-6-alt-h3_R2.fastq')}")

    return 1