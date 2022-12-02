import pandas as pd
from tqdm import tqdm
import numpy as np
import os

def add_library_info(df, df_library, verbose=False):

    if verbose: 
        print(f"Adding library info")

    # check the sanity of libraries.csv
    assert 'construct' in list(df_library.columns), "construct is not in library.csv"

    df_library.construct = df_library.construct.astype(str)
    if not sum([c in df_library.columns.tolist() for c in ['section_start','section_end']]) == 2:
        df = pd.merge(df, df_library, on='construct', how='left')
    else:
        new_df = pd.DataFrame()
        assert len(df['construct'].unique()) == len(df), 'constructs are not unique'
        if verbose:
            iter_fun = lambda x: tqdm(x.groupby('construct'),  total=len(df_library['construct'].unique()),desc='constructs')
        else:
            iter_fun = lambda x: x.groupby('construct')
        for idx, g in iter_fun(df_library):
            one_full_construct = False
            for _, row in g.iterrows():
                if row['construct'] not in df['construct'].tolist():
                    print(f'construct {row["construct"]} not in df')
                    continue
                row = pd.concat([row, df[df['construct']==row['construct']][[c for c in df.columns if c not in row.index]].iloc[0]])
                if 'Unnamed: 0' in row.index:
                    row = row.drop('Unnamed: 0')
                if not 'section' in df_library.columns:
                    row['section'] = row['section_start'] + '-' + row['section_end']
                row.rename({'section': 'section'}, inplace=True)
                isnan = lambda x: np.isnan(x) if isinstance(x, float) else False
                if isnan(row['section_start']) or isnan(row['section_end']):
                    if not one_full_construct:
                        one_full_construct = True
                        row['section'] = 'full'
                        row['section_start'], row['section_end'] = 0, len(row['sequence'])
                else:
                    row['section_start'], row['section_end'] = int(row['section_start']), int(row['section_end']) # CAST TO INT
                    for c in ['sequence']+[col for col in df.columns.tolist() if (col != 'num_of_mutations' and type(df[col].iloc[0]) == list)]:
                        row[c] = row[c][int(row['section_start']):int(row['section_end'])]
                        if isnan(row['section']):
                            row['section'] = '{}-{}'.format(row['section_start'], row['section_end'])
                new_df = pd.concat([new_df, pd.DataFrame(row).T])
        df = new_df.reset_index(drop=True)
        return df
