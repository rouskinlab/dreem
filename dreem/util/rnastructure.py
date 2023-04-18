import json
import os
import sys


# Set environment variable DATAPATH to RNAstructure data tables.
data_path_key = "DATAPATH"
if os.environ.get(data_path_key) is None:
    rs_dir = "rnastructure"
    dt_dir = "data_tables"
    env_dir = os.path.dirname(os.path.dirname(sys.executable))
    for root, dirs, files in os.walk(env_dir):  # FIXME: this might not be portable
        if (os.path.basename(root) == dt_dir
                and os.path.basename(os.path.dirname(root)) == rs_dir):
            # RNAstructure data tables were found in root.
            # Set DATAPATH environment variable to root.
            os.environ[data_path_key] = root
            break
    else:
        # RNAstructure data tables were not found.
        #raise FileNotFoundError(os.path.join(rs_dir, dt_dir))
        pass


def run_command(cmd):
    import subprocess
    print(cmd)
    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    return output.decode('utf-8')
                

class RNAstructure(): 
    def __init__(   self, 
                    rnastructure_path='', 
                    temp = 'temp/rnastructure/', 
                    rnastructure_use_temp = None,  
                    rnastucture_fold_args = None,
                    rnastructure_use_dms = None,
                    rnastructure_dms_min_unpaired_value = None,
                    rnastructure_dms_max_paired_value = None,
                    rnastructure_deltag_ensemble = None,
                    rnastructure_probability = None,
                    verbose = False,
                 ) -> None:
        self.rnastructure_path = rnastructure_path
        self.rnastructure_use_temp = rnastructure_use_temp
        self.rnastucture_fold_args = rnastucture_fold_args
        self.rnastructure_use_dms = rnastructure_use_dms
        self.rnastructure_dms_min_unpaired_value = rnastructure_dms_min_unpaired_value
        self.rnastructure_dms_max_paired_value = rnastructure_dms_max_paired_value
        self.rnastructure_deltag_ensemble = rnastructure_deltag_ensemble
        self.rnastructure_probability = rnastructure_probability
        self.verbose = verbose
        
        if self.rnastructure_use_dms != None:
            assert self.rnastructure_dms_min_unpaired_value != None, "If you want to use DMS, you need to specify the min and max values for the DMS data"
            assert self.rnastructure_dms_max_paired_value != None, "If you want to use DMS, you need to specify the min and max values for the DMS data"
            
        if len(self.rnastructure_path) > 0:
            if self.rnastructure_path[-1] != '/':
                self.rnastructure_path += '/'
        self.directory = temp
        if os.path.exists(os.path.join(self.directory, 'ledger.json')):
            self.ledger = json.load(open(os.path.join(self.directory, 'ledger.json'), 'r'))
        else:
            self.ledger = {}

    def __make_files(self, temp_prefix='temp'):
        self.pfs_file = os.path.join(self.directory, temp_prefix+'.pfs')
        self.ct_file = os.path.join(self.directory, temp_prefix+'.ct')
        self.dot_file = os.path.join(self.directory, temp_prefix+'.dot')
        self.fasta_file = os.path.join(self.directory, temp_prefix+'.fasta')
        self.prob_file = os.path.join(self.directory, temp_prefix+'_prob.txt')

    def fit(self, sequence, reference='reference'):
        self.sequence = sequence
        self.__make_temp_folder()
        self.__make_files()
        self.__create_fasta_file(reference, sequence)

    def predict_ensemble_energy(self):
        cmd = f"{self.rnastructure_path}EnsembleEnergy {self.fasta_file} --sequence"
        splitted_output = run_command(cmd).split(' ')
        return float(splitted_output[splitted_output.index(f"kcal/mol\n\nEnsemble")-1])

    def predict_partition(self, temperature_k =None):
        cmd = f"{self.rnastructure_path}partition {self.fasta_file} {self.pfs_file}"
        if temperature_k != None:
            cmd += ' --temperature '+str(temperature_k)
        run_command(cmd)
        run_command(self.rnastructure_path+'ProbabilityPlot '+ self.pfs_file + ' -t '+self.prob_file)
        with open(self.prob_file,"r") as f:
            lines=f.readlines()
            out={'i':[],'j':[],'p':[]}
            for x in range(len(lines)):
                if x>1:
                    ls=lines[x].split("\t")
                    out["i"]+=[int(ls[0])]
                    out["j"]+=[int(ls[1])]
                    out["p"]+=[float(ls[2])]
        return self.__cast_pairing_prob(out)

    def predict_reference_deltaG(self):
        cmd = f"{self.rnastructure_path}Fold {self.fasta_file} {self.ct_file}"
        run_command(cmd)
        assert os.path.exists(self.ct_file), f"{self.ct_file} does not exist, check that RNAstructure works. Command: {cmd}"
        assert os.path.getsize(self.ct_file) != 0, f"{self.ct_file} is empty, check that RNAstructure works. Command: {cmd}"
        return self.__extract_deltaG_struct()

    def __make_temp_folder(self):
        isExist = os.path.exists(self.directory)
        if not isExist:
            os.makedirs(self.directory)
        return self.directory

    def __create_fasta_file(self, reference, sequence):
        # push the ref into a temp file
        temp_fasta = open(self.fasta_file, 'w')
        temp_fasta.write('>'+reference+'\n'+sequence)
        temp_fasta.close()

    # cast the temp file into a dot_bracket structure and extract the attributes
    def __extract_deltaG_struct(self):
        run_command(f"{self.rnastructure_path}ct2dot {self.ct_file} 1 {self.dot_file}")
        temp_dot = open(self.dot_file, 'r')
        first_line = temp_dot.readline().split()
        # If only dots in the structure, no deltaG 

        if len(first_line) == 4:
            _, _, deltaG, _ = first_line
            deltaG = float(deltaG)
        if len(first_line) == 1:
            deltaG, _ = 0.0, first_line[0][1:]

        sequence = temp_dot.readline()[:-1] #  Remove the \n
        structure = temp_dot.readline()[:-1] # Remove the \n
        return deltaG, structure
    
    def dump_ledger(self):
        json.dump(self.ledger, open(os.path.join(self.directory, 'ledger.json'), 'w'))

    def run(self, sequence):
        if sequence in self.ledger:
            deltaG, structure = self.ledger[sequence]['deltaG'], self.ledger[sequence]['structure']
        else:
            self.fit(sequence)
            deltaG, structure = self.predict_reference_deltaG()
            self.ledger[sequence] = {'deltaG':deltaG, 'structure':structure}
            for file in os.listdir(self.directory):
                if file.startswith('temp'):
                    os.remove(os.path.join(self.directory, file))
        if deltaG == 'void': 
            deltaG = 0.0
        else:
            deltaG = float(deltaG)
        return {'deltaG':deltaG, 'structure':structure}

if __name__ == "__main__":
    rna = RNAstructure('/Users/ymdt/src/RNAstructure/exe/', temp = 'my_temp')
    rna.fit(sequence='AAGATATTCGAAACCACTCGATCGACTAGCATCAGCTGACTAGCTAGCATGCATCAAGAATATCTT')
    print("DeltaG + structure:", rna.predict_reference_deltaG())
    print("Ens. energy:", rna.predict_ensemble_energy())
    print("One line command:", rna.run('AAGATATTCGAAACCACTCGATCGACTAGCATCAGCTGACTAGCTAGCATGCATCAAGAATATCTT'))
    rna.dump_ledger()
