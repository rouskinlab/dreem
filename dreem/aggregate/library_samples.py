import pandas as pd
import yaml
import numpy as np

from dreem.aggregate.resources.get_attributes import read_sample_attributes

def get_samples_info(df_samples, sample, verbose= False):
    if verbose: print(f"Adding samples info for {sample}")

    # Sanity check
    # Keep only the columns that are in sample_attributes.yml
    df_samples = df_samples[df_samples['sample']==sample]

    assert len(df_samples) <= 1, f"{sample} has more than one line in samples.csv"
    assert len(df_samples) == 1, f"{sample} doesn't have a corresponding line in samples.csv"
    
    df_samples = df_samples.iloc[0]
    
    exp_env = df_samples['exp_env']
    assert exp_env in ['in_vivo','in_vitro'], f"{exp_env} is not a valid value for exp_env. Should be in_vivo or in_vitro"

    # Load list of mandatory columns and check that they are in samples.csv and not empty for this sample 
    sample_attributes = read_sample_attributes()
    for mand in sample_attributes['mandatory']['all'] + sample_attributes['mandatory'][exp_env]:
        assert mand in list(df_samples.index), f"{mand} is not in samples.csv"
        assert not df_samples[mand] == np.nan, f"{mand} is empty in samples.csv for sample {sample}"
        
    return df_samples.to_dict()

def get_library_info(df_library, construct, verbose= False):
    if verbose: print(f"Adding library info for {construct}")

    # Sanity check
    df_library = df_library[df_library['construct']==construct]

    df_library.drop(columns = [c for c in df_library.columns if c in ['section_start','section_end','section','construct']], inplace=True)
    
    return df_library.iloc[0].to_dict()

