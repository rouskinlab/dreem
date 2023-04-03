from ..util.seq import *
from ..util.util import *
from ..vector.profile import VectorReader


class BitVector:
    """Container object. Contains the name of the reference, the sequence, the bitvector, the read names and the read count.
    """
    
    def __init__(self, path, signal_thresh, include_gu, min_reads, include_del) -> None:
        preprocessing = self.preprocessing(path, signal_thresh=signal_thresh, include_gu=include_gu, min_reads=min_reads, include_del=include_del)

        self.name = path.split('/')[-1][:-(len('.json'))]
        self.sequence = preprocessing[0]
        self.bv = preprocessing[1]
        self.read_index = preprocessing[2]
        self.read_inverse = preprocessing[3]
        self.read_hist = preprocessing[4]
        self.read_names = preprocessing[5]
        self.report = preprocessing[6]
        self.base_to_keep = preprocessing[7]
        self.publish_preprocessing_report(path=path[:-len('.json')]+'_preprocessing_report.txt')
        
    #TODO optimize this 
    def preprocessing(self, path, signal_thresh = 0.005, include_gu = False, min_reads=1000, include_del=False, max_mut_close_by = 4):
        """Preprocess the bitvector.
        
        - Remove the bases G and U
        - Count and remove duplicates
        - Remove the bases with a Mutation fraction below signal_thresh
        - save the read names and the read count
        
        Parameters:
        -----------
        
        path: str
            Path to the directory containing the bitvector.
            
        signal_thresh: float
            Mutation fraction below which the bases are removed.
        
        include_gu: bool
            If True, keep the bases G and U.

        min_reads: int
            Minimum number of reads to keep a base.

        include_del: bool
            If True, keep the bases with a deletion.

        max_mut_close_by: int
            Maximum number of mutations close by to keep a base.

        
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

        reader = VectorReader.load(path)
        bv = reader.get_all_vectors()
        bv.columns = reader.columns

        report['total_number_of_reads'] = bv.shape[0]     
        
        # Take the read names
        read_names = np.array(bv.index.tolist(), dtype = str)
        
        ## PER BASE REMOVALS
        # Remove the non-informative bases types
        sequence = reader.seq.decode()
        report['sequence'] = sequence
        
        # Remove the G and U bases
        temp_n_cols = bv.shape[1]
        bases_to_drop = []
        if not include_gu:
            bases_to_drop = [c for c in bv.columns if c[0] in ['G','T']]
            bv = bv.drop(columns=bases_to_drop)
        report['removed_G_U'] = len(bases_to_drop)

        # Remove the low-mutation-rate bases
        bin_bv = mutations_bin_arr(bv)
        sub_rate = np.sum(bin_bv, axis = 0)/bin_bv.shape[0]
        bases_to_drop = [unpaired for unpaired, mut_rate in zip(bv.columns, sub_rate) if mut_rate < signal_thresh]
        bases_to_keep = set(bv.columns)-set(bases_to_drop)
        bases_to_keep = [int(i[1:]) for i in bases_to_keep]; bases_to_keep.sort()

        temp_n_cols = bv.shape[1]
        bv = bv.drop(columns=bases_to_drop)
        report['too_low_mutation_rate'] = temp_n_cols - bv.shape[1]
        report['min_mutation_rate'] = signal_thresh


        ## PER READ REMOVALS
        # Remove the bit mut_vectors with too many mutations
        temp_n_reads = bv.shape[0]
        bin_bv = mutations_bin_arr(bv)
        n_muts_per_reads = np.sum(bin_bv, axis = 1)
        idx = np.nonzero( n_muts_per_reads < np.mean(n_muts_per_reads) + 3*np.std(n_muts_per_reads) )[0]
        bv, read_names = bv.take(idx), read_names[idx]
        report['too_many_mutations'] = temp_n_reads - bv.shape[0]
        
        # Remove the bitvectors with deletions
        if not include_del:
            temp_n_reads = bv.shape[0]
            dels = deletions_bin_array(bv)
            no_deletion_in_the_read = np.nonzero(~np.sum(dels, axis=1, dtype=bool))[0]
            bv, read_names = bv.take(no_deletion_in_the_read), read_names[no_deletion_in_the_read]
            report['no_info_around_mutations'] = temp_n_reads - bv.shape[0]

        # MAKE BV BINARY
        for i, c in enumerate(bv.columns):
            bv[c] = mutations_bin_arr(bv[c])

        #Remove the bit mut_vectors with 'max_mut_close_by' consecutive mutations
        idx_remove_consecutive_mutations = []
        for i0, c0 in enumerate(bv.columns[:-1]):
            for c1 in bv.columns[i0+1:i0+max_mut_close_by+1]:
                if int(c1[1:]) - int(c0[1:]) <= max_mut_close_by:
                    idx_remove_consecutive_mutations += np.nonzero(np.logical_and(bv[c0], bv[c1]).values)[0].tolist()
        mask = np.ones(bv.shape[0], dtype=bool)
        mask[idx_remove_consecutive_mutations] = False
        temp_n_reads = bv.shape[0]
        bv, read_names = bv.take(np.arange(bv.shape[0])[mask]), read_names[mask]
        report['mutations_close_by'] = temp_n_reads - bv.shape[0] 

        # Turn bv into a np array
        bv = np.array(bv, dtype = np.uint8)    
        
        # What's this #TODO
        report['too_few_informative_bits'] = '#TODO'
        
        # Remove the duplicates and count the reads
        bv, read_idx, read_inverse, read_hist = np.unique(bv, axis = 0, return_index=True, return_inverse=True, return_counts = True)
        report['number_of_unique_reads'] = read_hist.shape[0]
        report['number_of_used_reads'] = np.sum(read_hist)
        report['bases_used'] = bv.shape[1]
        
        
        # Sanity check
        assert report['number_of_used_reads'] >= min_reads, "Too few reads after preprocessing"
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
        + "\nResidues removed because they have Mutation fraction below {} : {}".format(self.report['min_mutation_rate'], self.report['too_low_mutation_rate']) \
        + "\nNumber of bit mut_vectors in total: "+ str(self.report['total_number_of_reads'])\
        + "\nNumber of bit mut_vectors used: "+ str(self.report['number_of_used_reads']) \
        + "\nNumber of unique bit mut_vectors used: " + str(self.report['number_of_unique_reads']) \
        + "\nNumber of bit mut_vectors discarded: "+ str(self.report['total_number_of_reads'] - self.report['number_of_used_reads']) \
        + "\nBit mut_vectors removed because of too many mutations: " + str(self.report['too_many_mutations']) \
        + "\nBit mut_vectors removed because of too few informative bits: " + str(self.report['too_few_informative_bits']) \
        + "\nBit mut_vectors removed because of mutations close by: " + str(self.report['mutations_close_by']) \
        + "\nBit mut_vectors removed because of no info around mutations:  " + str(self.report['no_info_around_mutations'])
        
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
    return query_muts(np.array(bv, dtype=np.uint8), SUB_N | DELET | INS_5 | INS_3, sum_up=False)


def deletions_bin_array(bv):
    return query_muts(np.array(bv, dtype=np.uint8), DELET, sum_up=False)
