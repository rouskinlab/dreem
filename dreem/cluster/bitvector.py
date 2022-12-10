import numpy as np

import numpy as np
import pyarrow as pa
import pyarrow.orc as po
import pyarrow.compute as pc

class BitVector:
    """Container object. Contains the name of the construct, the sequence, the bitvector, the read names and the read count.
    """
    
    def __init__(self,path) -> None:
        preprocessing = self.preprocessing(path)

        self.name = path.split('/')[-1][:-(len('.orc'))]
        self.sequence = preprocessing[0]
        self.bv = preprocessing[1]
        self.read_hist = preprocessing[2]
        self.read_names = preprocessing[3]
        self.report = preprocessing[4]
        self.publish_preprocessing_report(path=path[:-len('.orc')]+'_preprocessing_report.txt')
        
    def preprocessing(self, path, low_mut_rate = 0.015, use_G_U = False, too_many_mutations = 5, max_mut_close_by = 3):
        """Preprocess the bitvector.
        
        - Remove the bases G and U
        - Count and remove duplicates
        - Remove the bases with a mutation rate below low_mut_rate
        - save the read names and the read count
        
        Parameters:
        -----------
        
        path: str
            Path to the directory containing the bitvector.
            
        low_mut_rate: float
            Mutation rate below which the bases are removed.
        
        use_G_U: bool
            If True, keep the bases G and U.
            
            
        Output:
        -------
        
        bv: array (N x D)
            Preprocessed bitvector.
            
        read_hist: array (N)
            Count of reads per bitvector.
            
        read_names: list of str
            Read names.
            
        report: str
            Report of the preprocessing.
            #TODO define
            
        """
        
        report = {}
        
        bv = po.read_table(path)
        report['total_number_of_reads'] = bv.shape[1]
        
        # Take the read names
        read_names = np.array(bv.column('id'), dtype = str)
        bv = bv.drop(['id'])
        
        # Remove the non-informative bases types
        sequence = ''.join([c[0] for c in bv.column_names])
        if not use_G_U:
            bv = bv.drop([c for c in bv.column_names if c[0] in ['G','T']])
        
        # Remove the low-mutation-rate bases
        mut_rates = np.sum(bv, axis = 1)/bv.shape[0]
        paired_bases = [unpaired for unpaired, mut_rate in zip(bv.column_names, mut_rates) if mut_rate < low_mut_rate]
        bv = bv.drop(paired_bases)

        # Remove the bit vectors with too many mutations
        temp_len = bv.shape[1]
        bv = bv.filter(np.sum(bv, axis = 0) < too_many_mutations)
        report['too_many_mutations'] = temp_len - bv.shape[1]
        
        #Remove the bit vectors with two consecutive mutations
        idx_remove_consecutive_mutations = []
        for i0, c0 in enumerate(bv.column_names[:-1]):
            for c1 in bv.column_names[i0+1:i0+max_mut_close_by+1]:
                if int(c1[1:]) - int(c0[1:]) <= max_mut_close_by:
                    break
                idx_remove_consecutive_mutations += np.nonzero(np.logical_xor(bv[c0], bv[c1]))[0].tolist()
        mask = np.ones(bv.shape[0], dtype=bool)
        mask[idx_remove_consecutive_mutations] = False
        report['mutations_close_by'] = bv.shape[0] - np.sum(mask)
        valid_reads = np.arange(bv.shape[0])[mask]
        bv = pc.take(bv, valid_reads)
        
        # Remove too few informative bits #TODO
        report['too_few_informative_bits'] = '#TODO'
        
        # Remove the bit vectors with no information around mutations #TODO
        report['no_info_around_mutations'] = '#TODO'  

        # Remove the duplicates and count the reads
        _, read_idx, read_hist = np.unique(bv, axis = 1, return_counts = True, return_index=True)
        read_names = [read_names[i] for i in read_idx]      
        report['number_of_unique_reads'] = bv.shape[1]
        report['number_of_used_reads'] = np.sum(read_hist)
       
        return sequence, bv, read_hist, read_names, report
    

   
   
    def publish_preprocessing_report(self, path):
       """Publish the report in a text file.
       
       Parameters:
       -----------
       
       path: str
           Path to the directory where the report will be saved.
           
       """
       
       report = "Number of bit vectors used: "+ str(self.report['number_of_used_reads']) \
       + "Number of unique bit vectors used: " + str(self.report['number_of_unique_reads']) \
       + "Number of bit vectors discarded: "+ str(self.report['total_number_of_reads'] - self.report['number_of_used_reads']) \
       + " Bit vectors removed because of too many mutations: " + str(self.report['too_many_mutations']) \
       + " Bit vectors removed because of too few informative bits: " + str(self.report['too_few_informative_bits']) \
       + " Bit vectors removed because of mutations close by: " + str(self.report['mutations_close_by']) \
       + " Bit vectors removed because of no info around mutations:  " + str(self.report['no_info_around_mutations']) 

       with open(path, 'w') as f:
           f.write(report)
           
           
    def associate_reads_with_likelihoods(self, likelihood_per_read, path):
       """Associates the reads with their likelihood, using the attributes read_names and read_hist. 
       Publish the reads in a json file.
       
       Parameters:
       -----------
       
        likelihood_per_read: array (N x K)
            Likelihood of each read.
        
        path: str
            Path to the directory where the report will be saved.
            
        Output:
        -------
        
        reads: dict
            Dictionary associating the read name with the likelihood.
        
       """
       
       reads = dict(placeholder = 'PLACEHOLDER #TODO')
       return reads
   
def drop_duplicates(table: pa.Table, column_name: str) -> pa.Table:
    unique_values = pc.unique(table[column_name])
    unique_indices = [pc.index(table[column_name], value).as_py() for value in unique_values]
    mask = np.full((len(table)), False)
    mask[unique_indices] = True
    return table.filter(mask=mask)