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
    def preprocessing(self, path, signal_thresh = 0.005, include_gu = False, min_reads=1000, include_del=False, min_mut_dist = 4):
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
        bv = reader.get_all_vectors(numeric=True)

        report['total_number_of_reads'] = bv.shape[0]
        
        ## PER BASE REMOVALS
        # Remove the non-informative bases types
        sequence = reader.seq.decode()
        report['sequence'] = sequence
        
        # Remove the G and U bases
        temp_n_cols = bv.shape[1]
        bases_to_drop = []
        if not include_gu:
            bases_to_drop = [pos for pos, base in zip(bv.columns,
                                                      reader.seq,
                                                      strict=True)
                             if base != A_INT and base != C_INT]
            bv.drop(columns=bases_to_drop, inplace=True)
        report['removed_G_U'] = len(bases_to_drop)

        ## PER READ REMOVALS
        # Remove the bit vectors with too many mutations
        temp_n_reads = bv.shape[0]
        n_muts_per_read = reader.query_vectors(bv, SUB_N | INDEL, subsets=True).sum(axis=1)
        n_muts_per_read_limit = n_muts_per_read.mean() + 3 * n_muts_per_read.std()
        bv.drop(index=bv.index[n_muts_per_read > n_muts_per_read_limit], inplace=True)
        report['too_many_mutations'] = temp_n_reads - bv.shape[0]
        
        # Remove the low-mutation-rate bases
        match_count = (reader.query_vectors(bv, MATCH).sum(axis=0) +
                       reader.query_vectors(bv, MATCH | INS_5).sum(axis=0))
        mut_query = SUB_N | INDEL if include_del else SUB_N
        mut_count = reader.query_vectors(bv, mut_query, subsets=True).sum(axis=0)
        mut_rate = mut_count / (match_count + mut_count)
        temp_n_cols = bv.shape[1]
        bv.drop(columns=bv.columns[mut_rate < signal_thresh], inplace=True)
        report['too_low_mutation_rate'] = temp_n_cols - bv.shape[1]
        report['min_mutation_rate'] = signal_thresh

        # MAKE BV BINARY
        mut_query = SUB_N | INDEL if include_del else SUB_N
        bv = reader.query_vectors(bv, mut_query, subsets=True)

        # Remove the bit mut_vectors with 'max_mut_close_by' consecutive mutations
        muts_too_close = pd.Series(np.zeros(bv.shape[0], dtype=bool),
                                   index=bv.index)
        for pos5 in bv.columns:
            pos3 = pos5 + min_mut_dist - 1
            muts_too_close |= bv.loc[:, pos5: pos3].sum(axis=1) > 1
        temp_n_reads = bv.shape[0]
        bv.drop(index=muts_too_close[muts_too_close].index, inplace=True)
        report['mutations_close_by'] = temp_n_reads - bv.shape[0]

        # Turn bv into a np array
        bases_to_keep = bv.columns - reader.end5  # convert to 0 indexing
        read_names = bv.index.to_list()
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
        assert report['total_number_of_reads'] == report['number_of_used_reads'] + report['too_many_mutations'] + report['mutations_close_by']
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
        + "\nBit mut_vectors removed because of mutations close by: " + str(self.report['mutations_close_by'])
        #+ "\nBit mut_vectors removed because of no info around mutations:  " + str(self.report['no_info_around_mutations']
        
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
        for name, idx in zip(self.read_names, self.read_inverse, strict=True):
            reads[name] = {}
            for k in range(likelihood_per_read.shape[1]):
                reads[name]['K'+str(k+1)] = likelihood_per_read[idx,k]

        return reads


def mutations_bin_arr(bv):
    return query_muts(np.array(bv, dtype=np.uint8), SUB_N | DELET | INS_5 | INS_3, sum_up=False)


def deletions_bin_array(bv):
    return query_muts(np.array(bv, dtype=np.uint8), DELET, sum_up=False)
