import numpy as np
import json

def cast_to_json_compat(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj

def cast_dict(mp):
    for k, v in mp.items():
        mp[k] = cast_to_json_compat(v)
        if type(v) is dict:
            cast_dict(v)
    return mp
    
def sort_dict(mut_profiles):
    sorting_key = lambda item: 1000*{int:0, str:0, float:0, list:1, dict:2}[type(item[1])] + ord(item[0][0])
    mut_profiles = {k:mut_profiles[k] for k in sorted(mut_profiles)}
    mut_profiles = dict(sorted(mut_profiles.items(), key=sorting_key))
    for k, v in mut_profiles.items():
        if type(v) is not dict:
            continue
        mut_profiles[k] = sort_dict(v)
    return mut_profiles


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)
    
def flatten_json(data):
    out, row = [], {}
        
    for k,v in data.copy().items():
        if type(v) != dict:
            row[k]= v
        else:
            row['reference'] = k
            for k2,v2 in v.items():
                if type(v2) != dict:
                    row[k2] = v2
                else:
                    row['section'] = k2
                    for k3,v3 in v2.items():
                        row['cluster'] = k3
                        if type(v3) != dict:
                            row[k3] = v3
                        else:
                            for k4,v4 in v3.items():
                                row[k4] = v4
                            out.append(row.copy())
    return out