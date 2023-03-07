import yaml, os

def read_sample_attributes():
    return {
        'mandatory': {
            'all': [
                'sample',
                'user',
                'date',
                'exp_env',
                'temperature_k',
                'inc_time_tot_secs',
                'DMS_conc_mM',
            ],
            'in_vitro': [
                'buffer',
            ],
            'in_vivo': [
                'cell_line',
            ],
        },
    }

