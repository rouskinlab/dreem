import scipy.stats
import numpy as np
from tqdm import tqdm
import pandas as pd

def compute_conf_interval(info_bases, mut_bases, alpha = 0.05):
    assert len(info_bases)==len(mut_bases), "info_bases and mut_bases must be of same length"
    ci = {}
    i, ci['poisson_min'], ci['poisson_max'], ci['poisson_low'], ci['poisson_high'] = 0, np.zeros(len(info_bases)), np.zeros(len(info_bases)), np.zeros(len(info_bases)), np.zeros(len(info_bases))
    for cov, mut in zip(info_bases, mut_bases):
        if cov == 0:
            ci['poisson_min'][i], ci['poisson_max'][i],ci['poisson_low'][i], ci['poisson_high'][i] = np.nan, np.nan, np.nan, np.nan
            continue 
        ci['poisson_min'][i], ci['poisson_max'][i] = 0.5*scipy.stats.chi2.ppf(alpha/2, df=2*mut)/cov, 0.5*scipy.stats.chi2.ppf(1-alpha/2, df=2*(mut+1))/cov
        if mut == 0:
            ci['poisson_min'][i] = 0
        ci['poisson_low'][i], ci['poisson_high'][i] = mut/cov-ci['poisson_min'][i], ci['poisson_max'][i]-mut/cov
        i = i+1
    ci.pop('poisson_min')
    ci.pop('poisson_max')
    ci['poisson_high'] = np.round(ci['poisson_high'], 4).tolist()
    ci['poisson_low'] = np.round(ci['poisson_low'], 4).tolist()
    return ci


def add_poisson_confidence_intervals(df, sample, verbose = False):
    ci = {}
    df.reset_index(inplace=True)
    if verbose:
        iter_fun = lambda x: tqdm(x.iterrows(), total=len(x), desc='Poisson intervals', postfix=sample)
    else:
        iter_fun = lambda x: x.iterrows()
    for idx, mh in iter_fun(df):
        ci[idx] = compute_conf_interval(info_bases=mh.info_bases, mut_bases=mh.mut_bases)
    df_ci = pd.DataFrame.from_dict(ci, orient='index')
    df = pd.concat([df, df_ci], axis=1)
    return df