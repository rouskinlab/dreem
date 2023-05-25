from collections import defaultdict
from logging import getLogger
import math

import numpy as np
import pandas as pd

from ..resources.get_attributes import read_sample_attributes

logger = getLogger(__name__)


def get_samples_info(df_samples, sample):
    logger.info(f"Adding samples info for {sample}")

    # Sanity check
    # Keep only the columns that are in sample_attributes.yml
    df_samples = df_samples[df_samples['sample'] == sample]

    assert len(df_samples) <= 1, f"{sample} has more than one line in samples.csv"
    assert len(df_samples) == 1, f"{sample} doesn't have a corresponding line in samples.csv"

    df_samples = df_samples.iloc[0]

    exp_env = df_samples['exp_env']
    assert exp_env in ['in_vivo',
                       'in_vitro'], f"{exp_env} is not a valid value for exp_env. Should be in_vivo or in_vitro"

    # Load list of mandatory columns and check that they are in samples.csv and not empty for this sample 
    sample_attributes = read_sample_attributes()
    for mand in sample_attributes['mandatory']['all'] + sample_attributes['mandatory'][exp_env]:
        assert mand in list(df_samples.index), f"{mand} is not in samples.csv"
        assert not df_samples[mand] == np.nan, f"{mand} is empty in samples.csv for sample {sample}"

    return df_samples.to_dict()


def get_library_info(df_library, reference):
    logger.info(f"Adding library info for {reference}")

    if df_library is None:
        return {}, {}
    # Sanity check
    df_library = df_library[df_library['reference'] == reference]

    section_names = {(int(start), int(end)): str(name)
                     for name, start, end in zip(df_library['section'],
                                                 df_library['section_start'],
                                                 df_library['section_end'],
                                                 strict=True)}

    df_library.drop(
        columns=[c for c in df_library.index if c in ['section_start',
                                                        'section_end',
                                                        'section',
                                                        'reference']],
        inplace=True)

    d = df_library.iloc[0].to_dict()
    for k, v in d.copy().items():
        if type(v) is float:
            if math.isnan(v):
                del d[k]

    return d, section_names


def get_library_sections(df_library: pd.DataFrame):
    sections: dict[tuple[str, int, int], str] = dict()
    for ref, end5, end3, sect in zip(df_library["reference"],
                                     df_library["section_start"],
                                     df_library["section_end"],
                                     df_library["section"],
                                     strict=True):
        try:
            # Convert the reference and section names to strings and the
            # 5' and 3' ends to integers.
            if not ref or (isinstance(ref, float) and math.isnan(ref)):
                raise ValueError(f"Missing ref for {sect}({end5}-{end3})")
            section = str(ref), int(end5), int(end3)
            if not sect or (isinstance(sect, float) and math.isnan(sect)):
                sname = ""
            else:
                sname = str(sect)
            # Check for duplicates.
            if (sname0 := sections.get(section)) is None:
                # The section has not already been named: name it.
                sections[section] = sname
            else:
                # The section has already been named. Check if its name
                # matches the names given previously.
                if sname0 == sname:
                    logger.warning(f"Section {section} ({sname}) was redefined")
                else:
                    logger.error(f"Section {section} named '{sname0}' was "
                                 f"redefined as '{sname}'; using first name")
        except Exception as error:
            logger.error(f"Failed to interpret section '{sect}' "
                         f"({ref}:{end5}-{end3}): {error}")
    return sections
