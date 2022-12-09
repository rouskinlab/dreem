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
        
    Key methods:
    ------------
    
    run: Run the clustering algorithm many times to iterate over K_max the number of clusters and N_runs the number of runs.
    
        
    """
    
    def __init__(self, construct, K_max, N_runs):
        self.construct = construct
        self.K_max = K_max
        self.N_runs = N_runs

    def run(self):
        '''Run the clustering algorithm many times to iterate over K_max the number of clusters and N_runs the number of runs.
        
        Output:
        -------

        The output is a dictionary containing the results of EMclustering for each K and each run.
        The first level of the dictionary is the number of clusters.
        The second level of the dictionary is a list containing the results of EMclustering for each run.
        The second levels of the dictionary are sorted by performance (log-likelihood).
        
        Example:
        --------
        
        results = {
            K1: [EMclustering1, EMclustering2, ...],
            K2: [EMclustering1, EMclustering2, ...],
            ...
        }
        
        The EMclustering output dictionary contains the following attributes:
            - mu: array (K x D)
                Mean of each cluster.
            - pi: array (K)
                Probability of each cluster.
            - log_likelihood: float
                Log-likelihood of the model.
        
        '''
        
        results = [{'mu': [], 'pi': [], 'log_likelihood': []} for _ in range(self.K_max)]
        return results
        
    
