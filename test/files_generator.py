import os, random, sys
import pandas as pd
import numpy as np
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from dreem.util.rnastructure import RNAstructure
from dreem.draw import Study
import json

test_file_folder = os.path.join(os.path.dirname((os.path.dirname(__file__))), 'test_files')

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



class BaseTestGenerator():
    def __init__(self, 
                 sample = 'my_test_sample',
                 path = os.getcwd(),
                 rnastructure_path = '',
                 demultiplexed = False,
                 ) -> None:
        
        self.sample = sample
        self.path = os.path.join(path, sample)
        self.rna = RNAstructure(rnastructure_path, temp=os.path.join(self.path, 'temp'))
        self.demultiplexed = demultiplexed

        # Sanity checks
        #self.remove_files()
        os.makedirs(self.path, exist_ok=True)
        self.reference_count = 0
        self.sequences_stack = []
        self.mutations_stack = []
        self.reference_stack = []
        self.min_hamming_distance = 10
        
        self.samples_csv_file = os.path.join(self.path, 'samples.csv')
        self.library_file = os.path.join(self.path, 'library.csv')
        self.json_file = os.path.join(self.path, '{}.json'.format(self.sample))
        self.fasta = os.path.join(self.path, self.sample + '.fasta')

    def make_fastq_names(self, prefix):
        return [os.path.join(self.path, '{}_R{}.fastq'.format(prefix, i)) for i in [1,2]]
        
    def remove_files(self):
        # assert that there's only json and fastq files in the folder
        if os.path.exists(self.path):
            for f in os.listdir(self.path):
                assert f.endswith('.json') or f.endswith('.fastq') or f.endswith('fasta') or f.endswith('.csv') or f == 'temp', "There are files in the folder that are not json or fastq files. Please remove them before running the test"
            os.system(' '.join(['rm', '-fr', self.path]))
            
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

class FakeData(BaseTestGenerator):

    def iterate_reference_count(self):
        self.reference_count += 1
        return self.reference_count
        

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
        fastq1_name, fastq2_name = self.make_fastq_names(self.reference_stack[-1] if self.demultiplexed else self.sample)
        fastq1_name, fastq2_name = os.path.join(self.path, fastq1_name), os.path.join(self.path, fastq2_name)
        
        # Write the fastq files
        with open(fastq1_name, 'a') as f1, open(fastq2_name, 'a') as f2:
            for i, m in enumerate(mutations):
                read = self.create_read(m)
                print_fastq_line(f1, '{}_read_{}:1:1:1:1:{}:{} 1:N:0:1'.format(self.reference_stack[-1], 1+i, self.reference_stack[-1][-1], i+1), read, 'F'*len(read))
                print_fastq_line(f2, '{}_read_{}'.format(self.reference_stack[-1], 1+i), invert_sequence(read), 'F'*len(read))
                
    def generate_fasta_files(self):
        """Write a fasta file with the given parameters
        
        Arguments:
            self.path {str} -- where to write the fasta files
            
        """
    
        # Write the fasta files
        with open(self.fasta, 'a') as f:
            print_fasta_line(f, self.reference_stack[-1], self.sequence)
        
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
            n_reads = len(mutations)
            out[reference] = {}
            out[reference]['num_reads'] = n_reads
            out[reference]['num_aligned'] = n_reads
            out[reference]['skips_short_reads'] = 0
            out[reference]['skips_too_many_muts'] = 0
            out[reference]['skips_low_mapq'] = 0
            out[reference]['some_random_attribute'] = library['some_random_attribute'].values[idx]
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
            while out[reference][section]['pop_avg']['sub_hist'][-1] == 0:
                out[reference][section]['pop_avg']['sub_hist'] = out[reference][section]['pop_avg']['sub_hist'][:-1]

            out[reference][section]['pop_avg']['del'] = [0]*(section_end-section_start+1)
            out[reference][section]['pop_avg']['ins'] = [0]*(section_end-section_start+1)
            out[reference][section]['pop_avg']['cov'] = [n_reads]*(section_end-section_start+1)
            out[reference][section]['pop_avg']['info'] = [n_reads]*(section_end-section_start+1)
            out[reference][section]['pop_avg']['sub_rate'] = np.array( np.array(out[reference][section]['pop_avg']['sub_N'])/np.array(out[reference][section]['pop_avg']['info'])).tolist()
            out[reference][section]['pop_avg']['min_cov'] = min(out[reference][section]['pop_avg']['cov'])
            # RNAstructure
            rna_structure_prediction = self.rna.run(out[reference][section]['sequence'])
            out[reference][section]['structure'] = rna_structure_prediction['structure']
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
        
        
class RealData(BaseTestGenerator):
    
    def open_fastq(self, reference):
        fn1, fn2 = self.make_fastq_names(reference if self.demultiplexed else self.sample)
        def make_fastq_path(fn):
            if self.demultiplexed:
                return os.path.join(self.path, self.sample, fn)
            else:
                return os.path.join(self.path, fn)
        f1, f2 = make_fastq_path(fn1), make_fastq_path(fn2)
        return open(f1, 'w' if self.demultiplexed else 'a'), open(f2, 'w' if self.demultiplexed else 'a')
    
    
    def add_library_and_fasta(self, ref, seq):
        # Add to fasta
        print_fasta_line(open(self.fasta,'a'), ref, seq)
        
        # Add to library
        df = pd.read_csv(self.library_file)
        row_full = {
            'reference': ref,
            'section': 'roi',
            'section_start': 70,
            'section_end': 81,
            'some_other_column': 'some_other_value'
        }   
        row_ms2 = {
            'reference': ref,
            'section': 'ms2',
            'section_start': 20,
            'section_end': 41,
            'some_other_column': 'some_other_value'
        }
        df = pd.concat([df, pd.DataFrame([row_full, row_ms2])])
        df.to_csv(self.library_file, index=False)
    
    def add_reference_3(self):
        self.add_library_and_fasta('reference_3', 'TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTCTTCCTTCGGGTCTCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTCTTCCTCTTTCTCTTTCTCTTTCTCTTCTCTTCTCTCTCTCTTCGTGAACGATTCTCGCTACTCGTTCCTTTCGA')
        f1, f2 = self.open_fastq('reference_3')
        f1.write('''@reference_3_read_1:1:1:1:1:3:1 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTCTTCCTTCGGGTCTCCACAGTCGAAAGACTGTGTCCTCTCTTCCTTTTTCTCTCCTCTTCTCTTTCTCTTTCCTTTCTTCTCTCTCCTTCGTGAACGATTCTCGCTAC
+
CCC;CCC;-CCCCCCCC-CCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCC;CCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCC;
@reference_3_read_2:1:1:1:1:3:2 1:N:0:1
TCGAAAGGAACGAGTAGCGATAATCGTTCACGAAAGAGAGAGAGAGAAAAGAAAGAGAAAGAGAAAGAGGAAGAGAAAAAGGAAGAGAGAGACACAGTCTTTCGACTGTGGAGACCCGAAGAAAGAGCATATGGGTGACCTCATATGCGGT
+
CCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@reference_3_read_3:1:1:1:1:3:3 1:N:0:1
TCGAAAGGAACGAGTAGCGAGAATCGTCCACGAAAGAGAGAGAGAGAAGAGAAAGAGAAAGAGAAAGAGGAAGAGAAAAAGGAAGAGACAGACACAGTCTTACGACTGTGGAGACCCGAAGGAAGAGCATATGGGTTATCCTCATATGCGG
+
CCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCC;C;CCCCCCCCCCCCCCCCCCCCCCC;CCC;CCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-C-CCCCCC;CCCC-CCCCCCCC-CCCCCCCCCCC
@reference_3_read_4:1:1:1:1:3:4 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTCTTCCTTCGGGTCTCCACAGTCGAAAGACTGTGTCTCTCTCTTCTTTTTCTCTCTCTCTCTCTTTCTCTTTCTCTTTCTCTTCTCTTCTCTCTCTCTTCGTGAACGAT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCC;CCC;CCCCCCCCCC-C-C;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@reference_3_read_5:1:1:1:1:3:5 1:N:0:1
TCGAAAGGAACGAGTAGCGAGAATCGTTCACAAAGGAGAGAGAAGAAAGGAAAGAGAAAGAAAGAGGAGAGAAAAAGGAAGAGAGGACACAGTCTTTCGACTGTGGAGACCCGAAGGAAGAGCATATGGGGGATCCTCATATGCGGTATGT
+
CCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCC;C
@reference_3_read_6:1:1:1:1:3:6 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGTCACCCATATGCTCTTCTTCGGGTCTCCACAGTCGAAAGACTGTGTCTCCTCTTCCTTTTTATCTTCCTCTTTCTTTTCCTTTTCTTTTCTTCCTCTCTCTTCGTGAATGATTCTCGCTACT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CC;CCCCCCCCCCCC;CCCCCCCCCCCCCCCCC;CCC-CCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@reference_3_read_7:1:1:1:1:3:7 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTCTTCCTTCGGGTCTCCACAGTCGAAAGACTGTGTCCTCTCTTCCTTTTTCTCTTCTCTTCTCTTTCTATTTCCTTTCTTCTCTCTCCTTAGTGAAAGATTCTCGCTAC
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCC;C;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCC;
@reference_3_read_8:1:1:1:1:3:8 1:N:0:1
TCGAAAGGAACGAGTAGCGAGAATCGTTCACGAAGAGAGTGGAAGAGAAAAAAGGAAAAGAAAGAGTAAGAGAAAAAGGAAGAGGAGACACAGTCATTCGACTGTGGAGACCCGAAGAAGAGCATATGGGTGAGCCTCATATGCGGTATGT
+
CCCCCC;CCCCCCCCCCCCCC-CCCC;-;C-CCCCCCC---CCC;---CCC;CC-CC-CC-CCCCC-CC-CC-CC;CCCCCCC--CCCCCC-;-CCCCC;CCCCCC-CCC;CC-C;CCCC;C;CCCCC-CCC;C;CCCCCC;CCC;;-CCC
@reference_3_read_9:1:1:1:1:3:9 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTTTTCTTCGGGTCTCCACAGTCGAAGACTGTGTCTCCTCTTCCTTTTTCTCTTCCTCTTTCTTTTCCTTTTCTTCTCTTCCTCTCTCTTCGTGAACGATTCTCGCTACT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;C;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCC;CCC;CCCCCC;CCCCCCCCCCCCCCCCC;CCCCCCCC;C-C;C;CC
@reference_3_read_10:1:1:1:1:3:10 1:N:0:1
TCGAAAGGAACGAGTAGCGAGAATCGTTCACGAAGAAAGAGAGAAGAGAAGAGAAAGAGAAAGAGAAAGAGGAAGAGAAAAAAGAAGAGAGACACACAGTCTTTCGACTGTGGAGACCCGAAGGAAGAGCAAGTGGGTGATCCTCATATGC
+
;CCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCC;CCCC;-CC--CCCCCCCCCCCCCCCCCCCCCCCCC
''')
            
        f2.write('''@reference_3_read_1:1:1:1:1:3:1 1:N:0:1
TCGAAAGGAACGAGTAGCGAGAATCGTTCACGAAGGAGAGAGAAGAAAGGAAAGAGAAAGAGAAGAGGAGAGAAAAAGGAAGAGAGGACACAGTCTTTCGACTGTGGAGACCCGAAGGAAGAGCATATGGGTGATCCTCATATGCGGTATG
+
CCCCCCCCCC;C;CCCCCCCC-CCC;CCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;C;CCCCCCCCCCCC;CCCCCCCCCCCCCCCC-CCCC;CCCCCC-CCCC;CCCCCCCCCCCCCC-CC;CCCCCC
@reference_3_read_2:1:1:1:1:3:2 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGTCACCCATATGCTCTTTCTTCGGGTCTCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTCTTCCTCTTTCTCTTTCTCTTTCTTTTCTCTCTCTCTCTTTCGTGAACGATTATCGC
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;;CCCCCCCCCC
@reference_3_read_3:1:1:1:1:3:3 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGATAACCCATATGCTCTTCCTTCGGGTCTCCACAGTCGTAAGACTGTGTCTGTCTCTTCCTTTTTCTCTTCCTCTTTCTCTTTCTCTTTCTCTTCTCTCTCTCTCTTTCGTGGACGATTCTCG
+
CCCCCCCC;CCCCC-CCCCCCCCCCCCCCCC;CCCCC;CCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCC;C;CCCCCCCCCC;CCCCC;CCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCC;CCCCCCC-C-CCCC
@reference_3_read_4:1:1:1:1:3:4 1:N:0:1
TCGAAAGGAACGAGTAGCGAGAATCGTTCACGAAGAGAGAGAGAAGAGAAGAGAAAGAGAAAGAGAAAGAGAGAGAGAGAGAAAAAGAAGAGAGAGACACAGTCTTTCGACTGTGGAGACCCGAAGGAAGAGCATATGGGTGATCCTCATA
+
CCCCCCCCCCCCCCCCCCCCC-CCCCCCC-CCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCC;CCCCCC;CCCCCCC;CCCCCCCCC;CCCCCC-CCCCCCC-CCCCC
@reference_3_read_5:1:1:1:1:3:5 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGATCCCCCATATGCTCTTCCTTCGGGTCTCCACAGTCGAAAGACTGTGTCCTCTCTTCCTTTTTCTCTCCTCTTTCTTTCTCTTTCCTTTCTTCTCTCTCCTTTGTGAACGATTCTCGCTACT
+
CCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@reference_3_read_6:1:1:1:1:3:6 1:N:0:1
TCGAAAGGAACGAGTAGCGAGAATCATTCACGAAGAGAGAGGAAGAAAAGAAAAGGAAAAGAAAGAGGAAGATAAAAAGGAAGAGGAGACACAGTCTTTCGACTGTGGAGACCCGAAGAAGAGCATATGGGTGACCTCATATGCGGTATGT
+
CCCCCCCCCCCC;CCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;C;CCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@reference_3_read_7:1:1:1:1:3:7 1:N:0:1
TCGAAAGGAACGAGTAGCGAGAATCTTTCACTAAGGAGAGAGAAGAAAGGAAATAGAAAGAGAAGAGAAGAGAAAAAGGAAGAGAGGACACAGTCTTTCGACTGTGGAGACCCGAAGGAAGAGCATATGGGTGATCCTCATATGCGGTATG
+
CCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCC;CCCCCCCCCCCCCCCCCC;CCC;CCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CC-CCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@reference_3_read_8:1:1:1:1:3:8 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGCTCACCCATATGCTCTTCTTCGGGTCTCCACAGTCGAATGACTGTGTCTCCTCTTCCTTTTTCTCTTCCTCTTTCTTTTCCTTTTTTCACTTCCTCTCTCATCGTGAACGATTCTCGCTACT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCC-CCCCC--CCCCCCCCC-CC;;;-;CC-CC-C--;;CCCC-CC;CC;CCCC;CCC-CCCCC-CCCCC-CCCC----CC--CCCCC-C--C-;CCC;CCCCC--C---
@reference_3_read_9:1:1:1:1:3:9 1:N:0:1
TCGAAAGGAACGAGTAGCGAGAATCGTTCACGAAGAGAGAGGAAGAGAAGAAAAGGAAAAGAAAGAGGAAGAGAAAAAGGAAGAGGAGACACAGTCTTCGACTGTGGAGACCCGAAGAAAAGCATATGGGTGATCCTCATATGCGGTATGT
+
CCCCCCCCCCCCCCCCCCCCC-CCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCC;CCCC;CCCCCCCCC;CC;C;CCCCCCCCCCCCCCCCCCCCCCC-C-CCCCCC;CCC
@reference_3_read_10:1:1:1:1:3:10 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGATCACCCACTTGCTCTTCCTTCGGGTCTCCACAGTCGAAAGACTGTGTGTCTCTCTTCTTTTTTCTCTTCCTCTTTCTCTTTCTCTTTCTCTTCTCTTCTCTCTTTCTTCGTGAACGATTCT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCC;CCCCCCCCCCCCCCCC;CCCCCCCCCCCCCC;CCCC-CCCCCCCCC;CCCCCCCCC;CCCC-C;CCC-;CC;C;CCC
''')

        self.update_json(
            ref='reference_3', 
            data = {
                "num_aligned": 10,
                "some_other_column": "some_other_value",
                "ms2": {
                "deltaG": -6.2,
                "section_end": 41,
                "section_start": 20,
                "sequence": "GCATATGAGGATCACCCATATG",
                "structure": ".((((((.((....))))))))",
                "pop_avg": {
                    "min_cov": 10,
                    "cov": [  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10],
                    "del": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0],
                    "info": [  10,  10,  10,  8,  10,  9,  10,  10,  8,  9,  7,  10,  10,  10,  10,  9,  9,  10,  8,  9,  10,  10],
                    "ins": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0],
                    "sub_A": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_C": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_G": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_N": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_T": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_hist": [  8,  2],
                    "sub_rate": [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.1,  0.1,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0]
                }
                },
                "roi": {
                "deltaG": -0.4,
                "section_end": 81,
                "section_start": 70,
                "sequence": "GACTGTGTCTCT",
                "structure": "(((...)))...",
                "pop_avg": {
                    "min_cov": 10,
                    "cov": [  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10],
                    "del": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  2],
                    "info": [  10,  7,  10,  9,  10,  10,  9,  10,  9,  6,  9,  7],
                    "ins": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_A": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_C": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_G": [  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0],
                    "sub_N": [  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0],
                    "sub_T": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_hist": [  8,  2],
                    "sub_rate": [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.1111111111111111,  0.0,  0.1111111111111111,  0.0]
                }}})

    def add_reference_4(self):
        
        self.add_library_and_fasta('reference_4', 'TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTCTGTCCCATGACTGGGATATCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTCTTCCTCTTTCTCTTTCTCTTTCTCTTCTCTTCTCTCTCATCAAGTACACCGCTACTCGTTCCTTTCGA')
        
        f1, f2 = self.open_fastq('reference_4')
        
        f1.write('''@reference_4_read_1:1:1:1:1:4:1 1:N:0:1
TCGAAAGGAACGAGTAGCGGTGTACTTGAAGAGAGAAGAGAGAGAGAGAGAAAGAGAAAGGGAAGAGAAAAAGGAAGGAGGACACAGTCTTTCGACTGTGGATATCCCAATCATGGGACAGGCATATGGGTGATCTCATATGCGGTATGTT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@reference_4_read_2:1:1:1:1:4:2 1:N:0:1
TTAACCGGCCAACATACCGCATATGAGGTTCACCCATATGCTCTGTCCCATGACTGGGATATCCACAGTCGAAAGACTGTGTCTCTCTTTCCTTTTTCTCTTCCTCTTTCTCTTTCTCTTTCTCTTCTCTTTCTCTCATCAAGTACACCGC
+
CCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCC;;CCCCCCCCCCCCCC;CC;;C
@reference_4_read_3:1:1:1:1:4:3 1:N:0:1
TCGAAAGGAACGAGTAGCGGTGTACTTGATGAGAGAGAAGAGAGAGAAAGAGAAAGAGAAAGGGAAGAGAAAAAGGAAGGAGGACACAGTCTATCGACTGTGGATATCCCAGTCATGGGACAGGCATATGGGTGATCTCATATGCGGTATG
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@reference_4_read_4:1:1:1:1:4:4 1:N:0:1
TCGAAAGGAACGAGTAGCGGTGTACTTGATGAGAGAGAAGAGAAGAGAAAGAGAAAGAGAAAGAGGAAGAGAAAAGGAAGAGATAGACACAGTCTTTCGACTGTGGATATCCCAGTCATGGGACAAAGCATATGGGTGATCCTCATATGCG
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCC;CCCCCCCC;CCCCCCCCCC
@reference_4_read_5:1:1:1:1:4:5 1:N:0:1
TTAAACCGGCCAACATACCGAATATGAGATCACCCATATGCCTGTCCCATGACTGGGATATCCACAGTCGAAAGACTGTGTCCTCCTTCCTTTTTCTCTTCCCTTTCTCTTTCTCTTTCTCTCTCTTCTCTCTCATCAAGTACAGCGCTAC
+
CCCCCCCCCCCCCCCCC;CCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@reference_4_read_6:1:1:1:1:4:6 1:N:0:1
TCGAAAGGAACGAGTAGCGGTGTACTTGATGAGAGAGAAGAGAGAGAAAGAGAAAGAGAAAGGGAAGAGAAAAAGGAAGGAGGACACAGTCTTTCGACTGTGGATATCCCAGTCATGGGACAGGCATATGGGTGATCTCATATGCGGTATG
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCC;CCC;CCCCCCCCCCCCCCCCCCCCCC;CCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;;CCCCCCCCCCCCCCCCCCCCC
@reference_4_read_7:1:1:1:1:4:7 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTCGGTCCCATGACTGGGATATCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCCTTCCTTTTTCTCTTTCTCTTTCTCTTCCTTCTCTTTCATCAAGTACACCC
+
;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCC-CCCC-CCCCC;CC;CCCCCCCCC;CCCCCC;CC;CCCCCCC;CC;CC;CCC-CCC;CCC-CCC;;;CCCCCCCCCCCCCCC;CCC-CCCCCC-
@reference_4_read_8:1:1:1:1:4:8 1:N:0:1
TCGAAAGGAACGAGTAGCGGTGTACTTGATGATAGAGAAGAGAAGAGAAAGAGAAAGAGAAAGAGGAACAGAAAAAGGAAGAGAGAGACACAGTCTTTCGACTGTGGATATCCCAGTCATGGGACAGAGCATATGGGTGATCCTCATATGC
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCC;CCCCCCCCCCCCCCCCCCCCC;CCCC-;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCC;CCCCCCCCCCCCCCCCCCC
@reference_4_read_9:1:1:1:1:4:9 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGATCACCCATATGCCTGTCCCATGAATGGGATATCCACAGTCGAAAGACTGTGTCCTCCTTCCTTTTTTCTTCCCTTTCTCTTTATCTTTCTCTCTCTTCTCTCTCATCAAGTACACCGCTACT
+
;CCCCCCCCCCCCCCCCCCCCCCC;CCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCC;CCCCC;CCCCCCCCCCCCC-CCCC;CCCCCC;;CCCCC
@reference_4_read_10:1:1:1:1:4:10 1:N:0:1
TCGAAAGGAACGAGTAGCGGTGTACATGATGAGAGAGAAGAGAAAAGAAAGAGAAAGAGAAAGAGGAAGAGAAAAAGGAAGAGAGAGACACAGTCTTTCGACTGTGGATATCCCAGTCATGGGACAGAGCATATGGGTGATCCTCATATGC
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
''')
        f2.write('''@reference_4_read_1:1:1:1:1:4:1 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGATCACCCATATGCCTGTCCCATGATTGGGATATCCACAGTCGAAAGACTGTGTCCTCCTTCCTTTTTCTCTTCCCTTTCTCTTTCTCTCTCTCTCTCTTCTCTCTTCAAGTACACCGCTACTC
+
CCCCCCCCCCCCCCCCCCCCCC;CCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCC;CCC
@reference_4_read_2:1:1:1:1:4:2 1:N:0:1
TCGAAAGGAACGAGTAGCGGTGTACTTGATGAGAGAAAGAGAAGAGAAAGAGAAAGAGAAAGAGGAAGAGAAAAAGGAAAGAGAGACACAGTCTTTCGACTGTGGATATCCCAGTCATGGGACAGAGCATATGGGTGAACCTCATATGCGG
+
C;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CC-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@reference_4_read_3:1:1:1:1:4:3 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGATCACCCATATGCCTGTCCCATGACTGGGATATCCACAGTCGATAGACTGTGTCCTCCTTCCTTTTTCTCTTCCCTTTCTCTTTCTCTTTCTCTCTCTTCTCTCTCATCAAGTACACCGCTAC
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCC;C;-
@reference_4_read_4:1:1:1:1:4:4 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTTTGTCCCATGACTGGGATATCCACAGTCGAAAGACTGTGTCTATCTCTTCCTTTTCTCTTCCTCTTTCTCTTTCTCTTTCTCTTCTCTTCTCTCTCATCAAGTACACC
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCC
@reference_4_read_5:1:1:1:1:4:5 1:N:0:1
TCGAAAGGAACGAGTAGCGCTGTACTTGATGAGAGAGAAGAGAGAGAAAGAGAAAGAGAAAGGGAAGAGAAAAAGGAAGGAGGACACAGTCTTTCGACTGTGGATATCCCAGTCATGGGACAGGCATATGGGTGATCTCATATTCGGTATG
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCC-C-CCCCCCCCCCCCCC
@reference_4_read_6:1:1:1:1:4:6 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGATCACCCATATGCCTGTCCCATGACTGGGATATCCACAGTCGAAAGACTGTGTCCTCCTTCCTTTTTCTCTTCCCTTTCTCTTTCTCTTTCTCTCTCTTCTCTCTCATCAAGTACACCGCTAC
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;-CCCC;CCCCCCC
@reference_4_read_7:1:1:1:1:4:7 1:N:0:1
TCGAAAGGAACGAGTAGCGGTGTACTTGATGAAAGAGAAGGAAGAGAAAGAGAAAGAGAAAAAGGAAGGAAAAAGGAAGAGAGAGACACAGTCTTTCGACTGTGGATATCCCAGTCATGGGACCGAGCATATGGGTGATCCTCATATGCGG
+
-C-CCCCCCCC;CCCC-CC;CCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCC;CCCCCCCC;CCCCCCCCCCCCCCC;CCCCCCCCCCCCCCC;CCCCCCCCCCCCCC-
@reference_4_read_8:1:1:1:1:4:8 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTCTGTCCCATGACTGGGATATCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTGTTCCTCTTTCTCTTTCTCTTTCTCTTCTCTTCTCTATCATCAAGTACAC
+
-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-C;CCCCC;CCCCCC;CCCCCCCCCCCCC;CCCCC;;C;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCC-C--CCCCC-CCCCCC;CCCCCC
@reference_4_read_9:1:1:1:1:4:9 1:N:0:1
TCGAAAGGAACGAGTAGCGGTGTACTTGATGAGAGAGAAGAGAGAGAAAGATAAAGAGAAAGGGAAGAAAAAAGGAAGGAGGACACAGTCTTTCGACTGTGGATATCCCATTCATGGGACAGGCATATGGGTGATCTCATATGCGGTATGT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCC;CCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCC;CCCC-CCCCCCCCCCCCCCCCCC;CCCCCCCCC;CCCCCCCCCCCC
@reference_4_read_10:1:1:1:1:4:10 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTCTGTCCCATGACTGGGATATCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTCTTCCTCTTTCTCTTTCTCTTTCTTTTCTCTTCTCTCTCATCATGTACAC
+
CCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCC;CCC-CCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCCCC-CCCC;CCCC-;CCCC
''')
            
        self.update_json(
            ref='reference_4', 
            data = {"num_aligned": 10,
                    "some_other_column": "some_other_value",
                    "ms2": {
                    "deltaG": -6.2,
                    "section_end": 41,
                    "section_start": 20,
                    "sequence": "GCATATGAGGATCACCCATATG",
                    "structure": ".((((((.((....))))))))",
                    "pop_avg": {
                        "min_cov": 10,
                        "cov": [  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10],
                        "del": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                        "info": [  10,  10,  10,  10,  10,  9,  10,  10,  5,  5,  10,  9,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10],
                        "ins": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                        "sub_A": [  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                        "sub_C": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                        "sub_G": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                        "sub_N": [  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                        "sub_T": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                        "sub_hist": [  8,  2],
                        "sub_rate": [  0.0,  0.1,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.1,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0]
                    }
                    },
                    "roi": {
                    "deltaG": 0.0,
                    "section_end": 81,
                    "section_start": 70,
                    "sequence": "TCGAAAGACTGT",
                    "structure": "............",
                    "pop_avg": {
                        "min_cov": 10,
                        "cov": [  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10,  10],
                        "del": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                        "info": [  10,  10,  10,  10,  10,  10,  10,  9,  10,  10,  10,  10],
                        "ins": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                        "sub_A": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                        "sub_C": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                        "sub_G": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                        "sub_N": [  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0],
                        "sub_T": [  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0],
                        "sub_hist": [  9,  1],
                        "sub_rate": [  0.0,  0.0,  0.0,  0.0,  0.1,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0]
                    }}})
            
    def add_reference_5(self):
        
        self.add_library_and_fasta('reference_5', 'TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTCGTCCCATGACTGGGATTCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTCTTCCTCTTTCTCTTTCTCTTTCTCTTCTCTTCTCTCTCTTATCTCAGTCTGCGCTACTCGTTCCTTTCGA')
        
        f1, f2 = self.open_fastq('reference_5')
        
        f1.write('''@reference_5_read_1:1:1:1:1:5:1 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGATTACCCATATGCTCGTCCCATGACTGGGATTCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTCTTCCTCTTTATCTTTCTCTTTTCTTCTCTTCTCTCTCTTATCTCAGTCTGC
+
CCCCCCCCCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@reference_5_read_2:1:1:1:1:5:2 1:N:0:1
TCGAAAGGAACGAGTAGCGCAGATTGAGATAAGAGAGAGAAGAGAAGAGAAAGAGAAAAAGAAAGAGGAAGAAAAAGGAAGAAAGAGACACAGTCTTTCGACTGTGGAATCCCAGTCATGGGACGAGCATATGGGTGATCCTCATATGCGG
+
CCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@reference_5_read_3:1:1:1:1:5:3 1:N:0:1
TCGAAAGGAACGAGTAGCGCAGATTGAGATAAGAGAGAGAAGAGAAGAGAAAGAGAAAGAGAAAGAGGAAGAGAAAAAGGAAGAGAGAGACACAGTCTTTCGACTGTGGAATCCCAGTCATGGGACAAGCATATGGGTGATCCTCATATGC
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCC;CCCCC
@reference_5_read_4:1:1:1:1:5:4 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGTCACCCATATGCTTGACCCATGACTGGGATCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTCTTCCTCTTTCTCTTTCTCTTTCTCTTCTTTCTCTCTCTTATCAAAGTCTGCGC
+
-CC;CCCC;;CC--CCCCC;-CCCCCCCC;-CCCC---C-C-CC-CCCCC-C;C;CCCCC;;CCCCC;-CC-CCC;C;CCCCCCCCC-;-CCCCC;CCCC;C;CC;C-CC;-C;C-;CC-C;C;CCCCCCC---CCC-----C-CCC;CCC
@reference_5_read_5:1:1:1:1:5:5 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTGTCCCATGACTGGGATTCCACAGCGAAAGCTGTGTCCTCTCTTCCTTTTTTCTTCTCTTTCTTTCTTCTTCTTCTCTCTCTTATCTCAGTCTGCGCTACTCGTTCCTT
+
CCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@reference_5_read_6:1:1:1:1:5:6 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTCGTCCCATGACTGGGATTCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTCTTCCTCTTTCTCTTTCTCTTTTCTTCTCTTCTCTCTCTTATCTCAGTCTGC
+
CCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@reference_5_read_7:1:1:1:1:5:7 1:N:0:1
CGAAAGGAACGAGTAGCGCAGACTGAGATAAGAGAGAGAAGAGGAAAAAGAGAAAGAGAAAAGGAAGAGAAAAAGGAAGAGAGAAACAGTCTTTCGACTGGGAATCCCAGTCATGGGACAAGCATATGGGTGATCCTCATATGCGGTATGT
+
CCCCCCCCCCCCCCCCCCCC;-CCCCCCCCCCCCCCCC-CCC-CCCCCCCCCCCCCC;CCCCC-CCCCCCCCCCCCCCCCCCC;CCCC;CCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCC;CCCCCCCC;;CCCCCCCCCCC;CCC
@reference_5_read_8:1:1:1:1:5:8 1:N:0:1
TCGAAAGGAACGAGTAGCGCAGACAGAGATAAGAGAGAGAAGAGAAGAGAAAGAGAAAGAGAGAGAGGAAGAGAAAAAGGAAGAGAGAGACACAGTCTTTCGACTGTGGAATCCCAGTCATGGGACGATCATATGGGTGATCCTCATATGC
+
CCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCC;CCCCCCCCCCC;CCCCCCCCCCC;CCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@reference_5_read_9:1:1:1:1:5:9 1:N:0:1
TTAAACCGGCCAACATACCGGATATGAGGATCACCCATATGCTCGTCCCATGACTGGTATTCCACAGACGAAAGACTGTGGCTCTCTCTTCCGTTTTCTGTTGCTCTTTCTCTTTTTATTTCTTTTCTCTTCTCTCTCTTATCTCTGTCTG
+
C;CCC-CC;CCCCCC-;CCC-CCCCCCCCCC-CCCCCCCC;CC;;CC-CCC-;C--;--;--;CC-C-C;CCCCCCC;C---CC;----C;--C;-;-C-;--C--C-;-;-CCC-C;-CC-C-C-CC;-;;;CCC;CC;C;-C-;-CC;C
@reference_5_read_10:1:1:1:1:5:10 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTTGTCCCTTGACTGGGATTCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTCTTCCTCTTTATCTTTCTCTTTCTCTTCTCTTCTCTCTCTTATCTAAGTCTG
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC''')
            
        f2.write('''@reference_5_read_1:1:1:1:1:5:1 1:N:0:1
TCGAAAGGAACGAGTAGCGCAGACTGAGATAAGAGAGAGAAGAGAAGAAAAGAGAAAGATAAAGAGGAAGAGAAAAAGGAAGAGAGAGACACAGTCTTTCGACTGTGGAATCCCAGTCATGGGACGAGCATATGGGTAATCCTCATATGCG
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCC;CC;CCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-CCCCCC;CC
@reference_5_read_2:1:1:1:1:5:2 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTCGTCCCATGACTGGGATTCCACAGTCGAAAGACTGTGTCTCTTTCTTCCTTTTTCTTCCTCTTTCTTTTTCTCTTTCTCTTCTCTTCTCTCTCTTATCTCAATCTGCG
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC-;;CCCC
@reference_5_read_3:1:1:1:1:5:3 1:N:0:1
ATAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTTGTCCCATGACTGGGATTCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTCTTCCTCTTTCTCTTTCTCTTTCTCTTCTCTTCTCTCTCTTATCTCAATCTG
+
-CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCC;C;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCC;;CCCCCC-CCCCCCCC
@reference_5_read_4:1:1:1:1:5:4 1:N:0:1
GCGCAAGGAACGAGTAGCGCAGACTGAGATAAGAGAGAGAAAGAAGAGAAAGAGAAAGAGAAAGAGGAAGAGAAAAAGGAAGAGAGAGACACAGTCTTTCGACTGTGGATCCCAGTCATGGGACAAGCATATGGGTGACCTCATATGCGGT
+
--C--;CC-CCCCCCCCCCCC-;C-CCCCCC-CCC-CCCC;C-;CC;C-C;--C-C;-;CCCCCC-C;CC-CCC;CCCCCCC-C;-C-C-C;CCCCCC;CC-CCC;-CCCCCCCCCCCCCCCCC;CCCC-C-;C-C;CCCC--C;;CCCC-
@reference_5_read_5:1:1:1:1:5:5 1:N:0:1
TCGAAAGGAACGAGTAGCGCAGACTGAGATAAGAGAGAGAAGAAGAAGAAAGAAAGAGAAGAAAAAAGGAAGAGAGGACACAGCTTTCGCTGTGGAATCCCAGTCATGGGACAGCATATGGGTGATCCTCATATGCGGTATGTTGGCCGGT
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
@reference_5_read_6:1:1:1:1:5:6 1:N:0:1
TCGAAAGGAACGAGTAGCGCAGACTGAGATAAGAGAGAGAAGAGAAGAAAAGAGAAAGAGAAAGAGGAAGAGAAAAAGGAAGAGAGAGACACAGTCTTTCGACTGTNGAATCCCAGTCATGGGACGAGCATATGGGTGATCCTCATATGCG
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCC;CCCCC;CCCCCCCCCCCCCCCCC#CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCC
@reference_5_read_7:1:1:1:1:5:7 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTTGTCCCATGACTGGGATTCCCAGTCGAAAGACTGTTTCTCTCTTCCTTTTTCTCTTCCTTTTCTCTTTCTCTTTTTCCTCTTCTCTCTCTTATCTCAGTCTGCGCTAC
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;C;CCCCCCCCCCC;C;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCC-CCCCCCCCCCCCCCC-CCCCC-CCCC
@reference_5_read_8:1:1:1:1:5:8 1:N:0:1
TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGATCGTCCCATGACTGGGATTCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTCTTCCTCTCTCTCTTTCTCTTTCTCTTCTCTTCTCTCTCTTATCTCTGTCTG
+
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;C;;CCC;C-CCCC;CCCCC;CCCCCC-CCCCCCCCCC;CCCC;CCCCCCCCCCCCCCCCCCCCCCCCCC-CCC-
@reference_5_read_9:1:1:1:1:5:9 1:N:0:1
TCGAAAGGAACGAGTAGCGCAGACAGAGTTAAGAGAGAGTAGAGAAAAGAAAGAAAAAGAGTAAGAGGAAGAGATAAAGGTAGAGACAGACACTGTCTTTCGACTGTGGAATCCCAGTCAAGTGACGAGCATATGGGAGATCATCATATCT
+
C-CCC-CC;CCCC-C;CCCC;C;-CCCC--CC----CC--CCC-;CCC--C-C-C;C-CC--C--C-C---CCC-CCC-C;-;-C---;CCCC-C-CC;-;;C;CC--;C-CC-CC-CCC-C-;--C;--CCC--;CCCCC----CCC---
@reference_5_read_10:1:1:1:1:5:10 1:N:0:1
CGAAAGGAACGAGTAGCGCAGACTTAGATAAGAGAGAGAAGAGAAGAGAAAGAGAAAGATAAAGAGGAAGAGAAAAAGGAAGAGAGAGACACAGTCTTTCGACTGTGGAATCCCAGTCAAGGGACAAGCATATGGGTGATCCTCATATGCG
+
CCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC;CCC;CCCCC;CCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC''')
        
        self.update_json(
            ref='reference_5', 
            data = {
                "num_aligned": 9,
                "some_other_column": "some_other_value",
                "ms2": {
                "deltaG": -6.2,
                "section_end": 41,
                "section_start": 20,
                "sequence": "GCATATGAGGATCACCCATATG",
                "structure": ".((((((.((....))))))))",
                "pop_avg": {
                    "min_cov": 9,
                    "cov": [  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9],
                    "del": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "info": [  9,  7,  9,  9,  9,  8,  8,  9,  8,  9,  8,  9,  7,  9,  8,  9,  9,  8,  8,  8,  9,  8],
                    "ins": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_A": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_C": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_G": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_N": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_T": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_hist": [  8,  1],
                    "sub_rate": [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.14285714285714285,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0]
                }
                },
                "roi": {
                "deltaG": 0.0,
                "section_end": 81,
                "section_start": 70,
                "sequence": "GAAAGACTGTGT",
                "structure": "............",
                "pop_avg": {
                    "min_cov": 9,
                    "cov": [  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9],
                    "del": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0],
                    "info": [  9,  8,  9,  9,  8,  8,  9,  9,  9,  9,  6,  8],
                    "ins": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_A": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_C": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_G": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_N": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_T": [  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0],
                    "sub_hist": [  9],
                    "sub_rate": [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0]
                }}})
            
        
    def update_json(self, ref, data):
        """Update the json file with the new data.
        """
        my_json = json.load(open(self.json_file,'r'))
        if ref in my_json:
            print('WARNING: the reference {} is already in the json file, the data will be overwritten.'.format(ref))
        my_json[ref] = data
        json.dump(my_json, open(self.json_file,'w'))
        


class TestFilesGenerator(RealData, FakeData):
    """The two inheraed classes are done a be widly, this is just a class to centralize the methods to run.
    """
    def run(self):
        """Generate all the test files.
        """
        # all
        self.generate_samples_file()
        
        # fake data specific
        self.add_reference(substitutions = [[4]]+[[9]]+[[14]]+[[19]]*4+[[24]]+[[29]]+[[34, 39]], L=50) # these are 0-indexed
        self.add_reference(substitutions = [[3]]+[[8]]+[[13]]+[[18]]*4+[[23]]+[[28]]+[[33, 38]], L=50) # KEEP A SINGLE VALUE FOR L
        self.generate_library_file()
        self.generate_json_file()
        
        # real data specific
        self.add_reference_3() # also add to the library file
        self.add_reference_4()
        self.add_reference_5()
        
if __name__ == '__main__':
    for sample in ['my_cli_sample','my_python_sample']:
        print('Generating test files to {}'.format(test_file_folder))
        TestFilesGenerator(
            path=test_file_folder, 
            sample=sample,
            rnastructure_path='/Users/ymdt/src/RNAstructure/exe',
            demultiplexed=False,
            ).run()

    # Plot the data 
#    s = Study(data = json.load(open(t.json_file, 'r')))
#    for sa, re, se in zip(s.df['sample'].values, s.df['reference'].values, s.df['section'].values):
#        s.mutation_fraction(sample=sa, reference=re, section=se)['fig'].show()

