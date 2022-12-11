from EMclustering import EMclustering
from bitvector import BitVector

class ClusteringAnalysis:
    """Launches the clustering algorithm many times to iterate over K_max the number of clusters and N_runs the number of runs.
    
    Parameters from args:
    ---------------------
        
    bitvector: object
        bitvector object. Conatins the processed bitvector and the count of reads per bitvector.
    
    K_max: int
        Number of clusters. The algorithm will iterate from 1 to K_max.
        
    N_runs: int
        Number of runs. The algorithm will be launched N_runs times.
        
    Key methods:
    ------------
    
    run: Run the clustering algorithm many times to iterate over K_max the number of clusters and N_runs the number of runs.
    plot: Plot the results of the clustering analysis. #TODO
        
    """
    
    def __init__(self, bitvector:BitVector, max_clusters:int, N_runs:int, clustering_args):
        self.bitvector = bitvector
        self.K_max = max_clusters
        self.N_runs = N_runs
        self.clustering_args = clustering_args
        # TODO add all

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
        results = {'K1': EMclustering(self.bitvector.bv, 1, self.bitvector.read_hist, **self.clustering_args)}
        for k in range(2,self.K_max+1):
            results['K'+str(k)] = []
            em = EMclustering(self.bitvector.bv, k, self.bitvector.read_hist, **self.clustering_args)
            for _ in range(self.N_runs):
                results['K'+str(k)].append(em.run())
                
        #results = {'K'+str(k):[{'mu': [], 'pi': [], 'log_likelihood': []} for _ in range(self.N_runs)] for k in range(self.K_max)}
        return results
        
    
