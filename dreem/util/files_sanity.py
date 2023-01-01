
from dreem.util.fa import parse_fasta
import numpy as np
from dreem.aggregate.resources.get_attributes import read_sample_attributes
import pandas as pd

def check_library(library: pd.DataFrame, path_fasta: str):
    """Sanity check for library.csv
    
    # sections
    - Check that there is only one value for each construct for each column that's not 'construct', 'section', 'section_start', 'section_end'
    - Check that every construct is in the fasta file
    - Check that there are no duplicate rows
    - If there are multiple sections, copy paste the other attributes so that every row is covered
    - Check that the section_start and section_end are in the correct order
    - Check that the section_start and section_end are within the length of the corresponding sequence in the fasta file.
    - If there are multiple sections, check that the section_start and section_end are unique. If not, remove the duplicates.
    - If section is empty, fill it with the section_start and section_end values separated by an underscore
    - If a construct has no sections, name the section 'full'
    
    # barcodes
    - Check that every barcode is in the fasta file at the given position for the given construct
    - Check that the barcode_start is within the length of the corresponding sequence in the fasta file.
    - Check that the barcode_start is unique for each construct
    - Check that each barcode is unique for each construct
    
    Args:
        library (pandas.DataFrame): content of library.csv
        path_fasta (str): path to the fasta file
        
    Returns:
        Library (pandas.DataFrame): reframed content of library.csv
    
    """
    

    a = parse_fasta(path_fasta)
    fasta = pd.DataFrame({k.decode("utf-8") :v.decode("utf-8")  for k,v in a}, index=[0]).T.reset_index().rename(columns={"index":"construct", 0:"sequence"})

    # Check that the constructs in the library are in the fasta file
    if not set(library["construct"].unique()).issubset(set(fasta["construct"].unique())):
        raise ValueError("Some constructs are not in the fasta file")

    # Sections
    # Check that there is only one value for each construct for each column that's not 'construct', 'section', 'section_start', 'section_end'
    for idx, g in library.groupby("construct"):
        for col in g.columns:
            if col not in ["construct", "section", "section_start", "section_end"]:
                un = g[col].unique()
                if np.count_nonzero(~pd.isnull(un)) > 1:
                    raise ValueError(f"{col} has multiple values for {idx}")
    
    # Check that every construct is in the fasta file
    if not set(library["construct"].unique()).issubset(set(fasta["construct"].unique())):
        raise ValueError("Some constructs are not in the fasta file")
    
    # Check that there are no duplicate rows and remove them
    if library.duplicated().any():
        library = library.drop_duplicates()
        
    # If there are multiple sections, copy paste the other attributes so that every row is covered
    if library["section"].nunique() > 1:
        for idx, g in library.groupby("construct"):
            for col in g.columns:
                if col not in ["construct", "section", "section_start", "section_end"]:
                    library.loc[library["construct"] == idx, col] = g[col].unique()[0]
                    
    # Make sure that barcode_start, section_start and section_end are integers
    for col in ["barcode_start", "section_start", "section_end"]:
        library[col] = library[col].apply(lambda x: int(x) if not pd.isnull(x) else x)
            
    # Check that the section_start and section_end are in the correct order
    if (library["section_start"] > library["section_end"]).any():
        raise ValueError("section_start is greater than section_end")
    
    # Check that section_end is within the length of the corresponding sequence in the fasta file.
    for idx, g in library.groupby("construct"):
        sequence = fasta.loc[fasta["construct"] == idx, "sequence"].unique()[0]
        if (g["section_end"] > len(sequence)).any():
            raise ValueError(f"section_end is greater than the length of the sequence for {idx} and section {max(g['section'].unique())} (sequence length: {len(sequence)})")
        
    # If there are multiple sections, check that the section_start and section_end are unique. If not, remove the duplicates.
    if library["section"].nunique() > 1:
        if library.duplicated(subset=["construct", "section_start", "section_end"]).any():
            library = library.drop_duplicates(subset=["construct", "section_start", "section_end"])
    
    # Assert no section is named 'full'
    if "full" in library["section"].unique():
        raise ValueError("Section cannot be named 'full'")
    
    # If section, section_start, and section_end are all empty for a certain row, fill section of this row with 'full' and set the section_start and section_end to 0 and the length of the sequence
    for idx, row in library.iterrows():
        if np.count_nonzero(pd.isnull(row[['section', 'section_start', 'section_end']])) == 3:
            library.loc[idx, "section"] = "full"
            library.loc[idx, "section_start"] = 0
            library.loc[idx, "section_end"] = len(fasta.loc[fasta["construct"] == row["construct"], "sequence"].unique()[0])
    
    # If section is empty but not the section start and end, fill it with the section_start and section_end values separated by an underscore
    if library["section"].isna().any():
        library.loc[library["section"].isna(), "section"] = library.loc[library["section"].isna(), "section_start"].astype(str) + "_" + library.loc[library["section"].isna(), "section_end"].astype(str)
    
    # Barcodes

    # Check that the barcode_start is unique for each construct
    if library.groupby("construct")["barcode_start"].nunique().max() > 1:
        raise ValueError("barcode_start is not unique for each construct")
        
    # Check that the barcode_start is within the length of the corresponding sequence in the fasta file.
    for idx, g in library.groupby("construct"):
        sequence = fasta[fasta["construct"] == idx]["sequence"].unique()[0]
        if (g["barcode_start"] + len(g['barcode'])> len(sequence)).any():
            raise ValueError(f"barcode_start + lenght of barcode is greater than the length of the sequence for {idx}")

    # Check that every barcode is in the fasta file at the given position for the given construct
    for idx, g in library.groupby("construct"):
        barcode = g["barcode"].unique()[0]
        barcode_start = g["barcode_start"].unique()[0]
        if barcode != fasta.loc[fasta["construct"] == idx, "sequence"].unique()[0][barcode_start:barcode_start+len(barcode)]:
            raise ValueError(f"Barcode {barcode} is not in the fasta file at the given position for {idx}")
        
    return library

def check_samples(samples):
    """Check that the samples csv file is correct
    
    - Check that the columns are correct
    - Check that "exp_env" is either "in_vivo" or "in_vitro"
    - Check that the sample names are unique
    - Check that no mandatory columns are missing or empty
    
    Args:
        samples (pandas.DataFrame): The samples csv file
        
    Returns:
        pandas.DataFrame: The samples csv file
    """
    
    samples_attributes = read_sample_attributes()['mandatory']
    mandatory_columns = samples_attributes['all']
    in_vivo_cols = samples_attributes['in_vivo']
    in_vitro_cols = samples_attributes['in_vitro']
    
    # Check that the columns are correct
    if not set(mandatory_columns).issubset(set(samples.columns)):
        missing_col = list(set(mandatory_columns) - set(samples.columns))
        raise ValueError(f"Some mandatory columns are missing. Columns are: {missing_col}")

    if 'in_vivo' in samples['exp_env'].unique():
        if not set(in_vivo_cols).issubset(set(samples.columns)):
            missing_col = list(set(in_vivo_cols + mandatory_columns) - set(samples.columns))
            raise ValueError(f"Some mandatory columns are missing. Columns are: {missing_col}")
        
    if 'in_vitro' in samples['exp_env'].unique():
        if not set(in_vitro_cols).issubset(set(samples.columns)):
            missing_col = list(set(in_vitro_cols + mandatory_columns) - set(samples.columns))
            raise ValueError(f"Some mandatory columns are missing. Columns are: {missing_col}")
    
    # Check that "exp_env" is either "in_vivo" or "in_vitro"
    if not set(samples['exp_env'].unique()).issubset(set(['in_vivo', 'in_vitro'])):
        raise ValueError("exp_env can only be 'in_vivo' or 'in_vitro'")
    
    for _, sample in samples.iterrows():
        # Check that the sample names are unique
        if samples[samples['sample'] == sample['sample']].shape[0] > 1:
            raise ValueError(f"Sample name {sample['sample']} is not unique")
        
        # Check that no mandatory columns are missing or empty
        if sample['exp_env'] == 'in_vivo':
            if sample[in_vivo_cols].isna().any():
                raise ValueError(f"Sample {sample['sample']} has missing or empty mandatory columns")
        elif sample['exp_env'] == 'in_vitro':
            if sample[in_vitro_cols].isna().any():
                raise ValueError(f"Sample {sample['sample']} has missing or empty mandatory columns")
        else:
            raise ValueError(f"Sample {sample['sample']} has an invalid sample_type")
    
    return samples
    
    
def compare_fields(d1, d2, fields):
    assert read_dict(d1, fields) == read_dict(d2, fields), 'fields {} are not equal'.format(fields)
    
def read_dict(d, fields):
    if fields != []:
        if fields[0] not in d:
            raise ValueError('Key {} not found. Check if it\'s implemented.'.format(fields[0]))
        return read_dict(d[fields[0]], fields[1:])
    return d