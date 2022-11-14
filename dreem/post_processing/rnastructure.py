
import os
import util
import numpy as np
import pandas as pd
                
class RNAstructure(object): #TODO
    def __init__(self, config) -> None:
        self.config = config
        self.rnastructure_path = config['rnastructure']['path'] if config['rnastructure']['path'][-1] == '/' else config['rnastructure']['path']+'/'

    def make_temp_folder(self, samp):
        temp_folder = 'temp/'+ samp + '/rnastructure/'
        isExist = os.path.exists(temp_folder)
        if not isExist:
            os.makedirs(temp_folder)
        return temp_folder

    def make_files(self, temp_prefix):
        self.pfs_file = f"{temp_prefix}.pfs"
        self.ct_file = f"{temp_prefix}.ct"
        self.dot_file = f"{temp_prefix}_dot.txt"
        self.fasta_file = temp_prefix+'.fasta'
        self.prob_file = temp_prefix+'_prob.txt'

    def create_fasta_file(self, construct, sequence):
        # push the ref into a temp file
        temp_fasta = open(self.fasta_file, 'w')
        temp_fasta.write('>'+construct+'\n'+sequence)
        temp_fasta.close()

    def predict_construct(self, use_dms, dms_file, use_temperature, temperature_k=None):
        # Run RNAstructure
        suffix = ''
        if use_temperature:
            suffix += f' --temperature {temperature_k}'
        if use_dms:
            suffix = f' -dms {dms_file}'
        cmd = f"{self.rnastructure_path}Fold {self.fasta_file} {self.ct_file} -d" + suffix
        util.run_command(cmd)
        assert os.path.getsize(self.ct_file) != 0, f"{self.ct_file} is empty, check that RNAstructure works"

    # cast the temp file into a dot_bracket structure and extract the attributes
    def extract_deltaG_struct(self):
        util.run_command(f"ct2dot {self.ct_file} 1 {self.dot_file}")
        temp_dot = open(self.dot_file, 'r')
        first_line = temp_dot.readline().split()
        # If only dots in the structure, no deltaG 
        out = {}
        if len(first_line) == 4:
            _, _, deltaG, _ = first_line
            deltaG = float(deltaG)
        if len(first_line) == 1:
            deltaG, _ = 'void', first_line[0][1:]

        sequence = temp_dot.readline()[:-1] #  Remove the \n
        structure = temp_dot.readline()[:-1] # Remove the \n
        return deltaG, structure

    def generate_normalized_mut_rates(self,temp_prefix, info_bases, mut_bases):
        mut_rates = np.array(mut_bases)/np.array(info_bases) 
        mut_rates = [max(r,self.config['rnastructure']['max_paired_mut_rate']) - self.config['rnastructure']['max_paired_mut_rate'] for r in mut_rates]     
        mut_rates = np.array([min(r,self.config['rnastructure']['min_unpaired_mut_rate']) for r in mut_rates])  
        pd.DataFrame((mut_rates)/(max(mut_rates)-min(mut_rates)),index=list(range(1,1+len(mut_rates))))\
                    .to_csv(temp_prefix+'_DMS_signal.txt', header=False, sep='\t')

    def predict_ensemble_energy(self):
        cmd = f"{self.rnastructure_path}EnsembleEnergy {self.fasta_file} --DNA --sequence"
        splitted_output = util.run_command(cmd)[0].split(' ')
        return float(splitted_output[splitted_output.index(f"kcal/mol\n\nEnsemble")-1])

    def predict_mut_probability(self, use_temperature, temperature_k):

        cmd = f"{self.rnastructure_path}partition {self.fasta_file} {self.pfs_file} --DNA"
        if use_temperature:
            cmd += ' --temperature '+str(temperature_k)
        util.run_command(cmd)
        util.run_command(self.rnastructure_path+'ProbabilityPlot '+ self.pfs_file + ' -t '+self.prob_file)
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

    def __cast_pairing_prob(self, prob:dict)->list:
        """Computes the pairing probability list.

        Args:
            prob (dict): Ouput of RNAstructure.
            index (int): index of the row that you want the pairing probabilities of.

        Returns:
            list: pairing probability, under the form of a list of probabilities
        """
        # Create local dataframe to play here
        pref_list, data_type = ['i', 'j', 'p'], {'i':int, 'j':int, 'p':float}
        df_loc = pd.DataFrame(prob)
        

        df_loc = pd.concat((df_loc, df_loc.rename(columns={'i':'j','j':'i'})))
        # Group the probabilities by i and get an ordered list of the pairing probability
        g = df_loc.groupby('i')['p'].agg(lambda row: [pow(10,-float(r)) for r in row])
        g.index = g.index.astype(int)
        g = g.sort_index()
        g['sum_log_p'] = g.apply(lambda row: sum(row))
        return list(g['sum_log_p'])

    def run(self, samp, mh, queue=None):
        out = {}
        temp_folder = self.make_temp_folder(samp)
        temp_prefix = f"{temp_folder}{mh.construct}_{mh.section}"
        self.generate_normalized_mut_rates(temp_prefix, mh.info_bases, mh.mut_bases)
        for temperature, temperature_suf in {False:'', True:'_T'}.items():
            if temperature and not self.config['rnastructure']['temperature']:
                continue
            this_sequence = mh.sequence
            for dms, dms_suf in {False:'', True:'_DMS'}.items():
                if dms and min(mh.info_bases) == 0:
                    continue
                suffix = dms_suf+temperature_suf
                self.make_files(f"{temp_prefix}{suffix}")
                self.create_fasta_file(mh.construct, this_sequence)
                self.predict_construct(use_dms = dms, dms_file=temp_prefix+"_DMS_signal.txt", use_temperature=temperature, temperature_k=mh.temperature_k)
                out['deltaG'+suffix], out['structure'+suffix] = self.extract_deltaG_struct()
                if not dms:
                    out['deltaG_ens'+suffix] = self.predict_ensemble_energy()
                if not dms  and not temperature:
                    out['mut_probability'+suffix] = self.predict_mut_probability(use_temperature=temperature, temperature_k=mh.temperature_k)
        if queue != None:
            queue.put(dict(sorted(out.items())))
        else:
            return dict(sorted(out.items()))