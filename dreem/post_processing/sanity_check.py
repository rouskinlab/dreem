
import pandas as pd
import os
import yaml
from util import Path, run_command


class Sanity_check(object):
    def __init__(self, config) -> None:
        self.samples = config['samples']
        self.path_to_dreem_output_files = config['path_to_dreem_output_files']
        self.sample_file = self.path_to_dreem_output_files+'samples.csv'
        self.library_file = self.path_to_dreem_output_files+'library.csv'
        self.path_to_rnastructure = config['rnastructure']['path']
        self.verbose = config['verbose']
        self.path = Path()
        self.config = config

    def files(self):
        # check that every file is there
        for s in self.samples:
            assert os.path.exists(self.path_to_dreem_output_files+s+'.csv'), f"{self.path_to_dreem_output_files+s+'.csv'} not found"
        if self.config['use_library']:
            assert os.path.exists(self.library_file), f"{self.library_file} not found"        
        if self.config['use_samples']:
            assert os.path.exists(self.sample_file), f"{self.sample_file} not found"
        if self.config['use_rnastructure']:
            assert os.path.exists(self.path_to_rnastructure), f"{self.path_to_rnastructure} not found"

    def check_samples(self):
        if self.verbose: print(f"Checking {self.sample_file}")
        # check the sanity of samples.csv
        df = pd.read_csv(self.sample_file)
        assert len(df['sample']) == len(df['sample'].unique()), "Every line isn't unique in samples.csv"

        # check that every sample as a corresponding line of samples.csv
        df['sample'] = df['sample'].astype(str)

        for s in self.samples:
            assert s in list(df['sample']), f"{s, type(s)} doesn't have a corresponding line in samples.csv"

        with open(self.path.sample_attribute_path, 'r') as f:
            sample_attributes = yaml.safe_load(f)

        assert len(df['exp_env'].unique()) == 1, "exp_env is not unique in samples.csv"
        exp_env = df['exp_env'].unique()[0]

        # check that every column of samples.csv is in sample_attributes.yml
        assert exp_env in ['in_vivo','in_vitro'], f"{exp_env} is not a valid value for exp_env. Should be in_vivo or in_vitro"
        # Check that you have all mandatory columns
        
        for mand in sample_attributes['mandatory']['all'] + sample_attributes['mandatory'][exp_env]:
            assert mand in list(df.columns), f"{mand} is not in samples.csv"
        
        # check that every mandatory column of samples.csv is not empty for every sample
        for mand in sample_attributes['mandatory']['all'] + sample_attributes['mandatory'][exp_env]:
            for s in self.samples:
                assert df[df['sample']==s][mand].isnull().sum() == 0, f"{mand} is empty in samples.csv for sample {s}"
            
        df.to_csv(self.config['temp_folder']+'/samples.csv', index=False)
        if self.verbose: print('Checking samples.csv done\n')
        return 1

    def check_library(self):
        # check the sanity of libraries.csv
        if self.verbose: print(f"Checking library.csv")
        assert os.path.exists(self.library_file), f"{self.library_file} not found"
        df = pd.read_csv(self.library_file)
        assert 'construct' in list(df.columns), "construct is not in library.csv"
        df.to_csv(f"{self.config['temp_folder']}/library.csv", index=False)
        if self.verbose: print(f"Checking library.csv done\n")

    def make_folder(self):
       if not os.path.exists( self.config['temp_folder']):
            os.makedirs( self.config['temp_folder'])

    def run(self):
        self.make_folder()
        if self.verbose: print("Checking files")
        self.files()
        if self.config['use_library']:
            self.check_library()
        if self.config['use_samples']:
            self.check_samples()
        if self.verbose: print("Checking files done\n")
        return 1