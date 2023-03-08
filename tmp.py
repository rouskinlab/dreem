import json
import numpy as np
import os

translation = {
    'mod_bases_A': 'sub_A',
    'mod_bases_C': 'sub_C',
    'mod_bases_G': 'sub_G',
    'mod_bases_T': 'sub_T',
    'mod_bases_N': 'sub_N',
    'mut_rates': 'sub_rate',
    'num_of_mutations': 'sub_hist',
    'info_bases': 'info',
    'cov_bases': 'cov',
    'ins_bases': 'ins',
    'del_bases': 'del',
    'mut_bases': 'sub_N',
    'worst_cov_bases': 'min_cov',
    'worst_cov': 'min_cov'
}


def update_json(json_file):
    
    with open(json_file) as f:
        data = json.load(f)
        
    for ref, v in data.items():
        if type(v) != dict:
            continue
        for section, v2 in v.items():
            if type(v2) != dict:
                continue
            for attr in v2['pop_avg'].copy().keys():
                if attr in translation.keys():
                    data[ref][section]['pop_avg'][translation[attr]] = data[ref][section]['pop_avg'].pop(attr)
                if attr in ['poisson_high','poisson_low']:
                    data[ref][section]['pop_avg'].pop(attr)
            L = len(data[ref][section]['pop_avg']['cov'])
            data[ref][section]['pop_avg']['sub_hist'] = np.histogram(data[ref][section]['pop_avg']['sub_hist'], bins=range(0, L, 1))[0].tolist()
            
    with open(json_file, 'w') as f:
        json.dump(data, f, indent=2)
    
    
path = '/Users/ymdt/Downloads/Archive (2)'

for s in os.listdir(path):
    update_json(os.path.join(path, s))