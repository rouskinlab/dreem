#!/usr/bin/env python 

import click
from click_option_group import optgroup
import os
from tqdm import tqdm

from dreem.post_processing import rnastructure, poisson
from dreem.post_processing.sanity_check import Sanity_check

import pandas as pd
import numpy as np

@click.command()
@optgroup.group("main arguments")
@optgroup.option("-c", "--config", type=click.Path(exists=True),
                 help="reference sequences in fasta format")
@optgroup.option("--samples_info", is_flag=True, help="Print the mandatory and optional columns for samples.csv")
@optgroup.option("--library_info", is_flag=True, help="Print the mandatory and optional columns for library.csv")
@optgroup.option("--generate_templates", default=None, help="Path to generate templates for samples.csv (in_vivo and in_vitro) and library.csv")

def main(**args):
    """
    DREEM processes DMS next generation sequencing data to produce mutational
    profiles that relate to DMS modification rates written by Silvi Rouskin and the
    Rouskin lab (https://www.rouskinlab.com/)
    """
    run(args)

def reformat_config(config):
    for name, val in config['use'].items():
        config['use_'+name] = val if config['use'] else False      

    config['temp_folder'] = 'temp/'
    config['samples'] = [str(s) for s in config['samples']]
    for path in [k for k in config if k.startswith('path')]:
        config[path] = config[path]+'/' if config[path][-1]!= '/' else config[path]

    return config

def make_dirs():
    for repo in ['temp']:
        if not os.path.exists(repo):
            os.makedirs(repo)

def verbose_print(s, config):
    if config['verbose']:
        print(s)

def load_csv(config, s):
    path = config['path_to_dreem_output_files']
    df = pd.read_csv(path + s +'.csv')
    for col in [ 'mut_bases', 'info_bases','del_bases','ins_bases','cov_bases','mut_rates'] + \
                    [c for c in df.columns.tolist() if c.startswith('mod_bases')]:
        df[col] = df[col].apply(lambda x: [float(b) for b in x[1:-1].replace('\n',' ').replace(',',' ').split(' ') if b != ''])
    return df

def add_samples(config, s, df):
    df_samp = pd.read_csv(config['temp_folder']+'/samples.csv').astype({'sample':str}).set_index('sample').loc[s]
    for col in df_samp.index:
        df[col] = df_samp[col]
    return df

def add_library(config, df):
    df_lib = pd.read_csv(config['temp_folder']+'/library.csv').astype({'construct':str})
    if not sum([c in df_lib.columns.tolist() for c in ['section_start','section_end']]) == 2:
        df = pd.merge(df, df_lib, on='construct', how='left')
    else:
        new_df = pd.DataFrame()
        assert len(df['construct'].unique()) == len(df), 'constructs are not unique'
        if config['verbose']:
            iter_fun = lambda x: tqdm(x.groupby('construct'),  total=len(df_lib['construct'].unique()),desc='constructs')
        else:
            iter_fun = lambda x: x.groupby('construct')
        for idx, g in iter_fun(df_lib):
            one_full_construct = False
            for _, row in g.iterrows():
                if row['construct'] not in df['construct'].tolist():
                    print(f'construct {row["construct"]} not in df')
                    continue
                row = pd.concat([row, df[df['construct']==row['construct']][[c for c in df.columns if c not in row.index]].iloc[0]])
                if 'Unnamed: 0' in row.index:
                    row = row.drop('Unnamed: 0')
                if not 'section_name' in df_lib.columns:
                    row['section_name'] = np.nan
                row.rename({'section_name': 'section'}, inplace=True)
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

def add_rnastructure(config, df,s):
    rna_pred = {}
    df.reset_index(inplace=True)
    """
    p,q,rna = {},{},{}
    MAX_PROCESS = config['rnastructure']['max_process']
    for idx, mh in tqdm(df.iterrows(), total=len(df), desc='RNAstructure start threads', postfix=s):
        rna[idx] = rnastructure.RNAstructure(config)
        q[idx] = Queue()
        p[idx] = Process(target=rna[idx].run, args=(s,mh,q[idx]))
        p[idx].start()
        while np.array(is_alive := [1 for x in p.values() if x.is_alive()]).sum() >MAX_PROCESS:
            time.sleep(0.1)
            for x, alive in enumerate(is_alive):
                if x<idx and not alive:
                    rna_pred[x] = q[x].get()
                    p[x].kill()
                    p[x], q[x] = None, None
    """
    rna = rnastructure.RNAstructure(config)
    if config['verbose']:
        iter_fun = lambda x: tqdm(x.iterrows(), total=len(x), desc='RNAstructure prediction', postfix=s)
    else:
        iter_fun = lambda x: x.iterrows()
    for idx, mh in iter_fun(df):
        rna_pred[idx] = rna.run(s,mh,config)
    df_rna = pd.DataFrame.from_dict(rna_pred, orient='index')
    df = pd.concat([df, df_rna], axis=1)
    return df

def add_poisson(config,df,s):
    ci = {}
    df.reset_index(inplace=True)
    if config['verbose']:
        iter_fun = lambda x: tqdm(x.iterrows(), total=len(x), desc='Poisson intervals', postfix=s)
    else:
        iter_fun = lambda x: x.iterrows()
    for idx, mh in iter_fun(df):
        ci[idx] = poisson.compute_conf_interval(info_bases=mh.info_bases, mut_bases=mh.mut_bases)
    df_ci = pd.DataFrame.from_dict(ci, orient='index')
    df = pd.concat([df, df_ci], axis=1)
    return df

def remove_index_columns(df):
    return df.drop(columns=[c for c in df.columns if c in ['level_0','index','Unnamed: 0']])


def run(config):
    config = reformat_config(config)
    make_dirs()
    Sanity_check(config).run()
    verbose_print('Starting add_info',config)

    for s in tqdm(config['samples'], total=len(config['samples']),desc='samples'):
        verbose_print(f'Read csv {s}.csv',config)
        df = load_csv(config, s)

        if config['use_samples']:
            verbose_print('Add samples',config)
            df = add_samples(config, s, df)

        if config['use_library']:
            verbose_print('Add library',config)
            df = add_library(config, df)

        if config['use_rnastructure']:
            verbose_print('Add RNAstructure prediction',config)
            df = add_rnastructure(config, df,s)

        if config['use_poisson']:
            verbose_print('Add Poisson intervals',config)
            df = add_poisson(config,df,s)

        df = remove_index_columns(df)

        if config['to_JSON']:
            verbose_print(f'Dump {s} to JSON',config)
            df.to_json(config['path_output']+'{}.json'.format(s))       

        if config['to_CSV']:
            verbose_print(f'Dump {s} to CSV',config)
            df.to_csv(config['path_output']+s+'.csv', index=True)

        if config['to_pickle']:
            verbose_print(f'Dump {s} to pickle',config)
            df.to_pickle(config['path_output']+s+'.p')

        verbose_print('Done with {}!'.format(s),config)
    verbose_print(f'Done!',config)
