import json
import pandas as pd

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
            for attr in v2['pop_avg'].keys():
                if attr in translation.keys():
                    data[ref][section]['pop_avg'][translation[attr]] = data[ref][section]['pop_avg'].pop(attr)
    
    with open(json_file, 'w') as f:
        json.dump(data, f, indent=2)
    
    
def update_feather(feather_file):
    
    df = pd.read_feather(feather_file)

    for attr in df.columns:
        if attr in translation.keys():
            if attr in ['mut_bases', 'poisson_high','poisson_low']:
                df.drop(attr, axis=1, inplace=True)
            df = df.rename(columns={attr: translation[attr]})
    df.to_feather(feather_file)
    

if __name__ == '__main__':
    update_feather('/Users/ymdt/dreem/dreem/dreem/draw/resources/my_dataset.feather')
