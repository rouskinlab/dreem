import pandas as pd
import numpy as np
import math 

from ..resources.get_attributes import read_sample_attributes

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

def get_library_info(df_library, reference, verbose= False):
    if verbose: print(f"Adding library info for {reference}")

    if df_library is None:
        return {},{}
    # Sanity check
    df_library = df_library[df_library['reference']==reference]
    
    section_translation = {str(ss)+'-'+str(se):s for s,ss,se in zip(df_library['section'],df_library['section_start'],df_library['section_end'])}

    df_library.drop(columns = [c for c in df_library.columns if c in ['section_start','section_end','section','reference']], inplace=True)
    
    d = df_library.iloc[0].to_dict()
    for k,v in d.copy().items():
        if type(v) is float:
            if math.isnan(v):
                del d[k]
    
    return d, section_translation
