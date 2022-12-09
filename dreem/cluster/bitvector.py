import numpy as np

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
                
        
    def preprocessing(self, path, low_mut_rate = 0.005):
        """Preprocess the bitvector.
        
        - Remove the bases G and U
        - Count and remove duplicates
        - Remove the bases with a mutation rate below low_mut_rate
        - save the read names and the read count
        
        Parameters:
        -----------
        
        path: str
            Path to the directory containing the bitvector.
            
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
        
        sequence = ''
        bv = np.zeros((0,0))
        read_hist = []
        read_names = []
        report = dict(
            place_holder = 'PLACEHOLDER #TODO'
        )
       
        
        return sequence, bv, read_hist, read_names, report
   
   
    def publish_preprocessing_report(self, path):
       """Publish the report in a text file.
       
       Parameters:
       -----------
       
       path: str
           Path to the directory where the report will be saved.
           
       """
       with open(path + self.name + '.txt', 'w') as f:
           f.write('PLACEHOLDER #TODO')
           
           
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
       