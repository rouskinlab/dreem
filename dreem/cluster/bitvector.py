    

class BitVector:
    """Container object. Contains the name of the construct, the sequence, the bitvector, the read names and the read count.
    """
    
    def __init__(self,path) -> None:
        preprocessing = self.preprocessing(path)

        self.name = path.split('/')[-1][:-(len('.orc'))]
        self.sequence = preprocessing[0]
        self.bv = preprocessing[1]
        self.reads_hist = preprocessing[2]
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
            
        reads_hist: array (N)
            Count of reads per bitvector.
            
        read_names: list of str
            Read names.
            
        report: str
            Report of the preprocessing.
            #TODO define
            
        """
        
       
        
       return sequence, bv, reads_hist, read_names, report
   
   
    def publish_report(self, path):
       """Publish the report in a text file.
       
       Parameters:
       -----------
       
       path: str
           Path to the directory where the report will be saved.
           
       """