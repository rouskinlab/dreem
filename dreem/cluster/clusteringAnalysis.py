from EMclustering import EMclustering

class ClusteringAnalysis:
    """Launches the clustering algorithm many times to iterate over K_max the number of clusters and N_runs the number of runs.
    
    Parameters from args:
    ---------------------
        
    construct: object
        Construct object. Conatins the processed bitvector and the count of reads per bitvector.
    
    K_max: int
        Number of clusters. The algorithm will iterate from 1 to K_max.
        
    N_runs: int
        Number of runs. The algorithm will be launched N_runs times.
        
    N_runs_save: int
        Number of runs to save. The algorithm will save the results of the N_run_save best runs.
        
    Key methods:
    ------------
    
    run: Run the clustering algorithm many times to iterate over K_max the number of clusters and N_runs the number of runs.
    
        
    """
    
    def __init__(self, construct, K_max, N_runs, N_runs_save):
        self.construct = construct
        self.K_max = K_max
        self.N_runs = N_runs
        self.N_run_save = N_runs_save

    def run(self):
        '''Run the clustering algorithm many times to iterate over K_max the number of clusters and N_runs the number of runs.
        
        Output:
        -------
        
        The output is a dictionary containing a list of tensors and their clustering score, each tensor being the result of a run of the clustering algorithm (the mutation profile).
        The dimension of the tensor is (N_runs_save, K, D). K is the number of clusters and D is the dimension of the bitvector.
        
        Example: 
        --------
        
        If you want a maximum of K_max=3 clusters and you want to save the results of the N_runs_save=2 best runs, the output will be a dictionary containing 3 lists of tensors.
        
        out = {
            K1: 
                score: clustering score
                data: tensor of dim [2, 1, D]
            K2:
                score: clustering score
                data: tensor of dim [2, 2, D]
            K3:
                score: clustering score
                data: tensor of dim [2, 3, D]
            }
        
        '''
        
        