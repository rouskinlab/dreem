from ..cluster.EMclustering import EMclustering
from ..cluster.bitvector import BitVector
import numpy as np
import multiprocessing

class ClusteringAnalysis:
    """Launches the clustering algorithm many times to iterate over K_max the number of clusters and num_runs the number of runs.
    
    Parameters from args:
    ---------------------
        
    bitvector: object
        bitvector object. Conatins the processed bitvector and the count of reads per bitvector.
    
    max_clusters: int
        Number of clusters. The algorithm will iterate from 1 to max_clusters.
        
    num_runs: int
        Number of runs. The algorithm will be launched num_runs times.
        
    Key methods:
    ------------
    
    run: Run the clustering algorithm many times to iterate over K_max the number of clusters and num_runs the number of runs.
        
    """
    
    def __init__(self, bitvector:BitVector, max_clusters:int, num_runs:int, clustering_args):
        self.bitvector = bitvector
        self.K_max = max_clusters
        self.num_runs = num_runs
        self.clustering_args = clustering_args
        # TODO add all

    def run(self):
        '''Run the clustering algorithm many times to iterate over K_max the number of clusters and num_runs the number of runs.
        
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

        # global dT # !! For test_input !!
        em = EMclustering(self.bitvector.bv, 1, self.bitvector.read_hist, self.bitvector.base_to_keep, self.bitvector.sequence,
                                **self.clustering_args)
        results = {1: [em.run()]}
        for k in range(2, self.K_max + 1):
            em = EMclustering(self.bitvector.bv, k, self.bitvector.read_hist, self.bitvector.base_to_keep, self.bitvector.sequence,
                                **self.clustering_args)

            pool = multiprocessing.Pool(processes=self.clustering_args["max_procs"])
            results[k] = sorted(pool.starmap(em.run, [() for _ in range(self.num_runs)]), key=lambda res: res['log_likelihood'], reverse=True)
            pool.close()
            pool.join()
        
        return results
        

if __name__ == '__main__':
    pass
    # if False:

    #     import matplotlib.pyplot as plt # !! For test_input !!
    #     import time, copy # !! For test_input !!
    #     import tracemalloc # !! For test_input !!
        
    #     dT = 0 # !! For test_input !!
        
    #     # dummy class just for test_input
    #     class testBV_class:
    #         def __init__(self, bit_vector, read_hist) -> None:
    #             self.bv = bit_Vector
    #             self.read_hist = read_hist
                

    #     exp_path = "/Users/Alberic/Desktop/Pro/RouskinLab/projects/DREEM/"
    #     bit_Vector_total = np.load(exp_path+"bit_vector.npy")
    #     read_hist_total =  np.load(exp_path+"read_hist.npy")

    #     with open(exp_path+"data_EM_analysis_v2.txt", 'a') as f:
    #         f.write("N_reads n_cpu dT mean_memory peak_memory \n")
    #         f.close()
        
    #     for N_partial in np.geomspace(10000, 120000, 5).astype(np.int64):

    #         N_partial = min(N_partial, bit_Vector_total.shape[0])

    #         for max_procs in range(1, 11):

    #             print("Starting experiment with {} reads and {} cpus".format(N_partial, max_procs))
                
    #             tracemalloc.start()

    #             # Get part of the total reads
    #             bit_Vector = copy.copy(bit_Vector_total[:N_partial])
    #             read_hist = copy.copy(read_hist_total[:N_partial])
    #             BV_class = testBV_class(bit_vector=bit_Vector, read_hist=read_hist)
                
    #             # Run EM and record stats
    #             clustering_args = dict(
    #                 min_iter = 10,
    #                 signal_thresh = 0.5, 
    #                 info_thresh = 0.5, 
    #                 include_g_u = True, 
    #                 include_del = True, 
    #                 min_reads = 10,
    #                 convergence_cutoff = 0.5,
    #                 num_runs = 100,
    #                 max_procs = max_procs,
    #                 verbose = True
    #             )

    #             clustering = ClusteringAnalysis(BV_class, 2, 10, clustering_args)
    #             clustering.run()
                
    #             # Log experiment stats
    #             memory_stats = [ mem/1e6 for mem in tracemalloc.get_traced_memory() ]
                
    #             print("Mean and max memory [MB]:", memory_stats)
    #             print("Total time [s]", "{:.1f}s".format(dT)) # !! For test_input !!
                
    #             log_data = np.array((N_partial,
    #                                     n_cpu,
    #                                     dT,
    #                                     memory_stats[0], 
    #                                     memory_stats[1] ))
    
    #             with open(exp_path+"data_EM_analysis_v2.txt", 'a') as f:
    #                 np.savetxt(f, log_data, newline=" ")
    #                 f.write("\n")
    #                 f.close()

    #             # Stop memory tracing and clean memory
    #             tracemalloc.stop()
    #             bit_Vector = None
    #             read_hist = None
    #             time.sleep(5)
                
    #             # Comparing output with previous software -> different test
    #             # mu_reference = np.load("/Users/Alberic/Desktop/Pro/RouskinLab/projects/DREEM/result.npy")
                # print("Reference matched:",(mu_reference == result["mu"]).all())
