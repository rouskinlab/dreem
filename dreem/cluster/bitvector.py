
import numpy as np
import pyarrow as pa
import pyarrow.orc as po
import pyarrow.compute as pc

from dreem.util.util import *

class BitVector:
    """Container object. Contains the name of the construct, the sequence, the bitvector, the read names and the read count.
    """
    
    def __init__(self,path, **args) -> None:
        preprocessing = self.preprocessing(path)

        self.name = path.split('/')[-1][:-(len('.orc'))]
        self.sequence = preprocessing[0]
        self.bv = preprocessing[1]
        self.read_index = preprocessing[2]
        self.read_inverse = preprocessing[3]
        self.read_hist = preprocessing[4]
        self.read_names = preprocessing[5]
        self.report = preprocessing[6]
        self.base_to_keep = preprocessing[7]
        self.publish_preprocessing_report(path=path[:-len('.orc')]+'_preprocessing_report.txt')
        
    #TODO optimize this 
    def preprocessing(self, path, low_mut_rate = 0.015, use_G_U = False, max_mut_close_by = 4):
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
        
        bv: numpy array (N x D)
            Preprocessed bitvector.
            
        read_hist: numpy array (N)
            Count of reads per bitvector.
            
        read_names: list of str
            Read names.
            
        report: str
            Report of the preprocessing.
            #TODO define
            
        """
        
        report = {}
        
        bv = po.read_table(path)
        report['total_number_of_reads'] = bv.shape[0]     
        
        # Take the read names
        read_names = np.array(bv.column('id'), dtype = str)
        bv = bv.drop(['id'])   
        
        ## PER BASE REMOVALS
        # Remove the non-informative bases types
        sequence = ''.join([c[0] for c in bv.column_names])  
        report['sequence'] = sequence
        
        temp_n_cols = bv.shape[1]
        bases_to_drop = set()
        if not use_G_U:
            # bv = bv.drop([c for c in bv.column_names if c[0] in ['G','T']])
            bases_to_drop |= set([c for c in bv.column_names if c[0] in ['G','T']])
        report['removed_G_U'] = temp_n_cols - bv.shape[1]

        # Remove the low-mutation-rate bases
        bin_bv = mutations_bin_arr(bv)
        mut_rates = np.sum(bin_bv, axis = 1)/bin_bv.shape[1]
        paired_bases = [unpaired for unpaired, mut_rate in zip(bv.column_names, mut_rates) if mut_rate < low_mut_rate]
        bases_to_drop |= set(paired_bases)
        bases_to_keep = set(bv.column_names)-bases_to_drop
        bases_to_keep = [int(i[1:]) for i in bases_to_keep]; bases_to_keep.sort()

        temp_n_cols = bv.shape[1]
        bv = bv.drop(bases_to_drop)
        report['too_low_mutation_rate'] = temp_n_cols - bv.shape[1]
        report['min_mutation_rate'] = low_mut_rate


        ## PER READ REMOVALS
        # Remove the bit vectors with too many mutations
        temp_n_reads = bv.shape[0]
        bin_bv = mutations_bin_arr(bv)
        n_muts_per_reads = np.sum(bin_bv, axis = 0)
        idx = np.nonzero( n_muts_per_reads < np.mean(n_muts_per_reads) + 3*np.std(n_muts_per_reads) )[0]
        bv, read_names = pc.take(bv, idx), read_names[idx]
        report['too_many_mutations'] = temp_n_reads - bv.shape[0]
        
        # Remove the bitvectors with deletions
        temp_n_reads = bv.shape[0]
        dels = deletions_bin_array(bv)
        no_deletion_in_the_read = np.nonzero(~np.sum(dels, axis=0, dtype=bool))[0]
        # MAKE BV BINARY
        for i, c in enumerate(bv.column_names):
            bv = bv.set_column(i, c, [mutations_bin_arr(bv.column(c))])
        bv, read_names = pc.take(bv, no_deletion_in_the_read), read_names[no_deletion_in_the_read]
        report['no_info_around_mutations'] = temp_n_reads - bv.shape[0]

        #Remove the bit vectors with 'max_mut_close_by' consecutive mutations
        idx_remove_consecutive_mutations = []
        for i0, c0 in enumerate(bv.column_names[:-1]):
            for c1 in bv.column_names[i0+1:i0+max_mut_close_by+1]:
                if int(c1[1:]) - int(c0[1:]) <= max_mut_close_by:
                    idx_remove_consecutive_mutations += np.nonzero(np.logical_and(bv[c0], bv[c1]))[0].tolist()
        mask = np.ones(bv.shape[0], dtype=bool)
        mask[idx_remove_consecutive_mutations] = False
        temp_n_reads = bv.shape[0]
        bv, read_names = pc.take(bv, np.arange(bv.shape[0])[mask]), read_names.take(np.arange(bv.shape[0])[mask])

        # Turn bv into a np array
        bv = np.array(bv, dtype = np.uint8).T        
        report['mutations_close_by'] = temp_n_reads - bv.shape[0] 
        
            
        # What's this #TODO
        report['too_few_informative_bits'] = '#TODO'
        
        # Remove the duplicates and count the reads
        bv, read_idx, read_inverse, read_hist = np.unique(bv, axis = 0, return_index=True, return_inverse=True, return_counts = True)
        report['number_of_unique_reads'] = read_hist.shape[0]
        report['number_of_used_reads'] = np.sum(read_hist)
        report['bases_used'] = bv.shape[1]
        
        
        # Sanity check
        assert len(report['sequence']) == report['bases_used'] + report['too_low_mutation_rate'] + report['removed_G_U']
        assert report['total_number_of_reads'] == report['number_of_used_reads'] + report['too_many_mutations'] + report['no_info_around_mutations'] + report['mutations_close_by']
        return sequence, bv, read_idx, read_inverse, read_hist, read_names, report, bases_to_keep
    

   
   
    def publish_preprocessing_report(self, path):
        """Publish the report in a text file.

        Parameters:
        -----------

        path: str
            Path to the directory where the report will be saved.
            
        """

        report = "Sequence: "+ self.report['sequence'] \
        + "\nLength of the sequence: " + str(len(self.report['sequence'])) \
        + "\nResidues used: " + str(self.report['bases_used']) \
        + "\nResidues removed because they are Gs and Us: " + str(self.report['removed_G_U']) \
        + "\nResidues removed because they have mutation rate below {} : {}".format(self.report['min_mutation_rate'], self.report['too_low_mutation_rate']) \
        + "\nNumber of bit vectors in total: "+ str(self.report['total_number_of_reads'])\
        + "\nNumber of bit vectors used: "+ str(self.report['number_of_used_reads']) \
        + "\nNumber of unique bit vectors used: " + str(self.report['number_of_unique_reads']) \
        + "\nNumber of bit vectors discarded: "+ str(self.report['total_number_of_reads'] - self.report['number_of_used_reads']) \
        + "\nBit vectors removed because of too many mutations: " + str(self.report['too_many_mutations']) \
        + "\nBit vectors removed because of too few informative bits: " + str(self.report['too_few_informative_bits']) \
        + "\nBit vectors removed because of mutations close by: " + str(self.report['mutations_close_by']) \
        + "\nBit vectors removed because of no info around mutations:  " + str(self.report['no_info_around_mutations']) 
        with open(path, 'w') as f:
            print(report)
            f.write(report)
           
           
    def associate_reads_with_likelihoods(self, likelihood_per_read):
        """Associates the reads with their likelihood, using the attributes read_names and read_hist. 
        Publish the reads in a json file.

        Parameters:
        -----------

        likelihood_per_read: array (N x K)
            Likelihood of each read.

        Output:
        -------

        reads: dict
            Dictionary associating the read name with the likelihood.

        """

        reads = {}
        for name, idx in zip(self.read_names, self.read_inverse):
            reads[name] = {}
            for k in range(likelihood_per_read.shape[1]):
                reads[name]['K'+str(k+1)] = likelihood_per_read[idx,k]

        return reads


def mutations_bin_arr(bv):
    return query_muts(np.array(bv, dtype=np.uint8), SUB_N[0] | DELET[0] | INS_5[0] | INS_3[0], sum_up=False)


def deletions_bin_array(bv):
    return query_muts(np.array(bv, dtype=np.uint8), DELET[0], sum_up=False)