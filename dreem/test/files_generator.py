import os, random, sys
import pandas as pd
import numpy as np
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from dreem.util.rnastructure import RNAstructure
from dreem.draw import Study
import json

test_file_folder = os.path.join(os.path.dirname((os.path.dirname(os.path.dirname(__file__)))), 'test_files')

def create_sequence(L):
    return ''.join([random.choice('ACGT') for _ in range(L)])

def invert_sequence(seq):
    return ''.join([{'A':'T','T':'A','C':'G','G':'C'}[s] for s in seq])[::-1]

def next_base(base):
    return {'A':'C','C':'G','G':'T','T':'A',0:1}[base]

def prev_base(base):
    return {'A':'T','C':'A','G':'C','T':'G',0:1}[base]

def hamming_distance(s1, s2):
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def print_fasta_line(f, id, seq):
    f.write('>{}\n{}\n'.format(id, seq))

def print_fastq_line(f, id, seq, qual):
    f.write('@{}\n{}\n+\n{}\n'.format(id, seq, qual))    
    
def sort_dict(mut_profiles):
    sorting_key = lambda item: 1000*{int:0, str:0, float:0, list:1, dict:2}[type(item[1])] + ord(item[0][0])
    mut_profiles = {k:mut_profiles[k] for k in sorted(mut_profiles)}
    mut_profiles = dict(sorted(mut_profiles.items(), key=sorting_key))
    for k, v in mut_profiles.items():
        if type(v) is not dict:
            continue
        mut_profiles[k] = sort_dict(v)
    return mut_profiles

def flatten_list_of_lists(v):
    return [item for sublist in v for item in sublist]

def empty_list_of_lists(L):
    return [[] for _ in range(L)]

def vector_of_zeros(L):
    return [0]*L

def count_pop_avg(mutations, sequence, section_start, section_end, base='ACGT'):
    # turn a list of lists into a list
    if base == 'N':
        base = 'ACGT'
    mutations = flatten_list_of_lists(mutations)
    out = [0]*len(sequence)

    for i in mutations:
        if next_base(sequence[i]) in base and section_start <= i < section_end:
            out[i] += 1
    return out
        
def count_mutations_per_read(mutations, section_start, section_end):
    # turn a list of lists into a list
    return [len([m for m in m_ if section_start <= m < section_end]) for m_ in mutations]


class Test_files_generator():
    
    def __init__(self, 
                 sample = 'my_test_sample',
                 path = os.getcwd(),
                 rnastructure_path = ''
                 ) -> None:
        
        self.sample = sample
        self.path = path
        self.rna = RNAstructure(rnastructure_path)

        # Sanity checks
        self.remove_files()
        os.makedirs(self.path, exist_ok=True)
        self.reference_count = 0
        self.sequences_stack = []
        self.mutations_stack = []
        self.reference_stack = []
        self.min_hamming_distance = 10
        
        self.samples_csv_file = os.path.join(self.path, 'samples.csv')
        self.library_file = os.path.join(self.path, 'library.csv')
        self.json_file = os.path.join(self.path, '{}.json'.format(self.sample))
        

    def iterate_reference_count(self):
        self.reference_count += 1
        return self.reference_count
        
    def remove_files(self):
        # assert that there's only json and fastq files in the folder
        if os.path.exists(self.path):
            for f in os.listdir(self.path):
                assert f.endswith('.json') or f.endswith('.fastq') or f.endswith('fasta') or f.endswith('.csv') or f == 'temp', "There are files in the folder that are not json or fastq files. Please remove them before running the test"
            os.system(' '.join(['rm', '-fr', self.path]))

    def create_read(self, mutations):
        read = self.sequence
        for m in mutations:
            read = read[:m] + next_base(read[m]) + read[m+1:]
        return read

    def assert_mutations_is_a_list_of_tuples(self, mutations):
        assert isinstance(mutations, list)
        for m in mutations:
            assert isinstance(m, tuple) or isinstance(m, list)
            for mm in m:
                assert isinstance(mm, int)
                assert mm < len(self.sequence)
                assert mm >= 0
            
    def generate_fastq_files(self, mutations):
        """Write a fastq file with the given parameters
        
        Arguments:
            self.path {str} -- where to write the fastq files
            self.mutations {dict} -- dictionary of mutations to be written in the fastq file
            
        """
        
        # Name the fastq files
        prefix = os.path.join(self.path, self.sample)
        fastq1_name, fastq2_name = prefix + '_R1.fastq', prefix + '_R2.fastq'

        # Write the fastq files
        with open(fastq1_name, 'a') as f1, open(fastq2_name, 'a') as f2:
            for i, m in enumerate(mutations):
                read = self.create_read(m)
                print_fastq_line(f1, '{}_read_{}'.format(self.reference_stack[-1], 1+i), read, 'F'*len(read))
                print_fastq_line(f2, '{}_read_{}'.format(self.reference_stack[-1], 1+i), invert_sequence(read), 'F'*len(read))
                
    def generate_fasta_files(self):
        """Write a fasta file with the given parameters
        
        Arguments:
            self.path {str} -- where to write the fasta files
            
        """
        
        # Name the fasta files
        file = os.path.join(self.path, self.sample + '.fasta')

        # Write the fasta files
        with open(file, 'a') as f:
            print_fasta_line(f, self.reference_stack[-1], self.sequence)
    
    def generate_samples_file(self):
        pd.DataFrame({
            'sample': [self.sample],
            'user': ['Napoleon'],
            'date': ['1769-08-15'],
            'exp_env': ['in_vitro'],
            'temperature_k': [310],
            'buffer': ['MgCl2'],
            'DMS_conc_mM': [2],
            'inc_time_tot_secs': [100]
        }, index=[0]).to_csv(self.samples_csv_file, index=False)
        
    def generate_library_file(self):
        df = {}
        for reference, sequence in zip(self.reference_stack, self.sequences_stack):
            df[reference] = {}
            df[reference]['reference'] = reference
            df[reference]['section_start'] = 1
            df[reference]['section_end'] = len(sequence)
            df[reference]['section'] = "{}-{}".format(1, len(sequence))
            df[reference]['some_random_attribute'] = 'some_random_value_{}'.format(reference)
        df = pd.DataFrame(df).T
        df.to_csv(self.library_file, index=False)
        
                
    def generate_json_file(self):
        
        if not os.path.isfile(self.samples_csv_file):
            self.generate_samples_file()
        samples_csv = pd.read_csv(self.samples_csv_file)
        samples_csv['temperature_k'] = samples_csv['temperature_k'].astype(int)
        
        if not os.path.isfile(self.library_file):
            self.generate_library_file()
        library = pd.read_csv(self.library_file)
        
        out = samples_csv.to_dict('records')[0]
        for idx, (reference, sequence, mutations) in enumerate(zip(self.reference_stack, self.sequences_stack, self.mutations_stack)):
            print(reference, sequence, mutations)
            n_reads = len(mutations)
            out[reference] = {}
            out[reference]['num_reads'] = n_reads
            out[reference]['num_aligned'] = n_reads
            out[reference]['skips_short_reads'] = 0
            out[reference]['skips_too_many_muts'] = 0
            out[reference]['skips_low_mapq'] = 0
            out[reference]['some_random_attribute'] = library['some_random_attribute'].values[idx]
            out[reference]['sequence'] = sequence   
            section = library['section'].values[idx]
            section_start = library['section_start'].values[idx]
            section_end = library['section_end'].values[idx]
            out[reference][section] = {}
            out[reference][section]['section_start'] = int(section_start)
            out[reference][section]['section_end'] = int(section_end)
            out[reference][section]['sequence'] = sequence[section_start-1:section_end]
            out[reference][section]['pop_avg'] = {}
            for base in ['A', 'C', 'G', 'T','N']:
                out[reference][section]['pop_avg']['sub_{}'.format(base)] = count_pop_avg(mutations, sequence, section_start-1, section_end, base)
            out[reference][section]['pop_avg']['sub_N'] = count_pop_avg(mutations, sequence, section_start-1, section_end, 'N')
            out[reference][section]['pop_avg']['sub_hist'] = np.histogram(count_mutations_per_read(mutations, section_start-1, section_end), bins=range(0, len(out[reference][section]['sequence'])))[0].tolist()
            out[reference][section]['pop_avg']['del'] = [0]*(section_end-section_start+1)
            out[reference][section]['pop_avg']['ins'] = [0]*(section_end-section_start+1)
            out[reference][section]['pop_avg']['cov'] = [n_reads]*(section_end-section_start+1)
            out[reference][section]['pop_avg']['info'] = [n_reads]*(section_end-section_start+1)
            out[reference][section]['pop_avg']['sub_rate'] = np.array( np.array(out[reference][section]['pop_avg']['sub_N'])/np.array(out[reference][section]['pop_avg']['info'])).tolist()
            out[reference][section]['pop_avg']['min_cov'] = min(out[reference][section]['pop_avg']['cov'])
            # RNAstructure
            rna_structure_prediction = self.rna.run(out[reference][section]['sequence'])
            out[reference][section]['structure'] = rna_structure_prediction['structure']
            print(rna_structure_prediction)
            print(out[reference][section]['sequence'])
            out[reference][section]['deltaG'] = float(rna_structure_prediction['deltaG'])
            out = sort_dict(out)
        with open(self.json_file, 'w') as f:
            json.dump(out, f, indent=2)
    
    def safe_sequence_generation(self, L):
        self.sequence = create_sequence(L)
        for past_sequence in self.sequences_stack:
            if len(past_sequence) != len(self.sequence):
                continue
            if hamming_distance(self.sequence, past_sequence) < self.min_hamming_distance:
                self.safe_sequence_generation(L)
        self.stack_sequence(self.sequence)
        return self.sequence
    
    def stack_sequence(self, sequence):
        self.sequences_stack.append(sequence)
        
    def stack_mutations(self, mutations):
        self.mutations_stack.append(mutations)
            
    def stack_reference(self):
        count = self.iterate_reference_count()
        self.reference_stack.append('reference_{}'.format(count))
            
    def add_reference(self, substitutions, L = 100):
        # Check and save inputs
        self.sequence = self.safe_sequence_generation(L)
        self.stack_mutations(substitutions)
        self.stack_reference()
        self.assert_mutations_is_a_list_of_tuples(substitutions)
        # Generate files
        self.generate_fastq_files(substitutions)
        self.generate_fasta_files()
        
if __name__ == '__main__':
    
    print('Generating test files to {}'.format(test_file_folder))
    
    t = Test_files_generator(
        path=test_file_folder, 
        sample='my_test_sample',
        rnastructure_path='/Users/ymdt/src/RNAstructure/exe',
        )
    t.generate_samples_file()
    t.add_reference(substitutions = [[4]]+[[9]]+[[14]]+[[19]]*4+[[24]]+[[29]]+[[34, 39]], L=50)
    t.add_reference(substitutions = [[3]]+[[8]]+[[13]]+[[18]]*4+[[23]]+[[28]]+[[33, 38]], L=50) # KEEP A SINGLE VALUE FOR L
    t.generate_library_file()
    t.generate_json_file()

    # Plot the data 
    s = Study(data = json.load(open(t.json_file, 'r')))
    for sa, re, se in zip(s.df['sample'].values, s.df['reference'].values, s.df['section'].values):
        s.mutation_fraction(sample=sa, reference=re, section=se)['fig'].show()
