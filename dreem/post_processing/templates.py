from get_info import get_attributes
from util import Path
import os
import itertools

class TemplateGenerator(object):
    def __init__(self, path):
        if path == '.':
            path = os.path.abspath('') + '/'
        if path[0] != '/':
            self.path = os.path.abspath('') + '/' + path + '/'
        else:
            self.path = path if path[-1] == '/' else path+'/'
        if not os.path.exists(path):
            os.makedirs(path)

    def _write_cols_to_csv(self, file,all_cols):
        chain = ''
        for col in all_cols:
            chain = chain+ col+','
        with open(file,'w') as f:
            f.write(chain[:-1]+'\n'*3)
            f.close()

    def generate_template_samples(self, exp_env:str):
        attributes = get_attributes('samples.csv','yml')
        assert exp_env in ['in_vivo','in_vitro'], "exp_env must be 'in_vivo' or 'in_vitro'"
        all_cols = [attributes[a][b] for a in ['mandatory','optional'] for b in ['all',exp_env]][:-1]
        all_cols = [item for sublist in all_cols for item in sublist if item is not None]
        self._write_cols_to_csv(f"{self.path}template_samples_{exp_env}.csv",all_cols)
        print(f"{self.path}template_samples_{exp_env}.csv")

    def generate_template_library(self):
        attributes = get_attributes('library.csv','yml')
        all_cols = attributes['mandatory']['all'] + attributes['optional']['all']
        self._write_cols_to_csv(f"{self.path}template_library.csv",all_cols)
        print(f"{self.path}template_library.csv")    

    def generate_config_template(self):
        p = Path()
        with open(p.config_template,'r') as f:
            temp = f.read()
            f.close()
        with open(f"{self.path}template_config.yml",'w') as f:
            f.write(temp)
            f.close()
        print(f"{self.path}template_config.yml")

    def generate_studies_template(self):
        with open(self.path+'template_studies.yml','w') as f:
            f.write('name,description,samples,conditions,title\n\n')
            f.close()
        print(f"{self.path}template_studies.csv")

    def run(self):
        print("Generating templates...")
        self.generate_template_samples('in_vivo')
        self.generate_template_samples('in_vitro')
        self.generate_template_library()
        self.generate_config_template()
        self.generate_studies_template()
