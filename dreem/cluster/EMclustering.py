import numpy as np
import scipy.stats
from scipy.optimize import newton_krylov

from ..util.cli import *

## ------------- Utility functions ------------- ##

def calc_BIC(N, PARAMS_LEN, K, log_like):
    """
    """
    return np.log(N) * PARAMS_LEN * K - (2 * log_like)

def calc_denom(i, mu, denom_probs, s2_probs):
    """
    """
    mu = np.clip(mu, 0.0, 1.0)
    if i in denom_probs:  # Already encountered
        return (denom_probs[i], s2_probs)
    elif i >= len(mu):  # Base case
        return (1, s2_probs)
    else:  # Make the calc
        s1 = calc_denom(i + 1, mu, denom_probs, s2_probs)[0]
        s2 = (1.0 - mu[i + 1: i + 4]).prod() * \
            calc_denom(i + 4, mu, denom_probs, s2_probs)[0]
        denom_probs[i] = ((1 - mu[i]) * s1) + (mu[i] * s2)
        s2_probs[i] = s2
        return (denom_probs[i], s2_probs)

def mu_der(mu_k, x_bar_k):
    """
    """
    mu_k_rev = mu_k[::-1]
    denom_k = calc_denom(0, mu_k, {}, {})
    denom_k_rev = calc_denom(0, mu_k_rev, {}, {})
    upd_mu = [(mu_k[i] * denom_k[1][i] * denom_k_rev[1][len(mu_k) - i - 1] /
              denom_k[0]) - x_bar_k[i]
              for i in range(len(mu_k))]
    return np.array(upd_mu)

    

class EMclustering:
    """This class runs the EM clustering algorithm.
    
    Parameters from args:
    ---------------------
    
    bv: array (N x D)
        Bitvector of the reference. 
    
    K: int
        Number of clusters.
        
    read_hist: array (N)
        Number of reads per bitvector.

    bases_to_keep: array (D)
        Array of booleans indicating which bases to keep.
    
    sequence: str
        Reference sequence.

    max_procs: int
        Maximum number of processes to use.

    min_iter: int
        Minimum number of iterations.

    max_iter: int
        Maximum number of iterations.
        
    convergence_cutoff: float
        Convergence threshold. When the difference between the log-likelihood of two consecutive iterations is lower than convergence_cutoff, the algorithm stops.
        
    Key methods:
    ------------
    
    run: Run the EM clustering algorithm. 
    maximize: Maximization step of the EM algorithm.
    expectation: Expectation step of the EM algorithm.
    
    """
    
    def __init__(self, bv, K, read_hist, bases_to_keep, sequence, max_procs: int, max_clusters: int, signal_thresh: float, include_gu: bool, include_del: bool, min_reads: int, min_iter: int, max_iter:int, convergence_cutoff: float, num_runs: int, verbose: bool = True):
        self.bv = bv
        self.sparse_mu = np.zeros((K, len(sequence)))
        self.bases_to_keep = bases_to_keep
        self.K = K
        self.read_hist = read_hist
        self.N = bv.shape[0]
        self.D = bv.shape[1]
        
        self.min_iter = min_iter
        self.max_iter = max_iter

        self.convergence_cutoff = convergence_cutoff

        self.verbose = verbose
    
    def expectation(self, mu, pi):
        """
        Run the Expectation step of the EM algorithm - calc log like
        Args:
            X (EM_Class object): Contains the list of bit mut_vectors
            K (int): Number of clusters
            mu (list): DMS reactivities in each cluster
            pi (list): Proportion of each cluster
            calc_inds (list): Indices to split the matrix at
            max_procs (int): Number of processes to use
            max_clusters (int): Maximum number of clusters
            signal_thresh (float): Threshold for signal
            info_thresh (float): Threshold for information
            include_g_u (bool): Include G and U in the calculation
            include_del (bool): Include deletions in the calculation
            min_reads (int): Minimum number of reads for a cluster
            min_iter (int): Minimum number of iterations
            convergence_cutoff (float): Convergence threshold
            num_runs (int): Number of runs
            verbose (bool): Print verbose output
            
        Returns:
            resps (list): Responsibity of each cluster
            log_like (float): Log likelihood of observing the data
            denom (float): Corrected denom of EM algorithm
        """
        N, D = self.bv.shape[0], self.bv.shape[1]
        log_pi = np.log(pi)
        log_pmf = np.zeros((N, D, self.K))
        for k in range(self.K):
            self.sparse_mu[k][self.bases_to_keep] = mu[k]
        denom = [calc_denom(0, self.sparse_mu[k], {}, {})[0] for k in range(self.K)]
        

        # Compute probability mass function for all reads, for each cluster
        for k in range(self.K):
            log_pmf[:, :, k] = scipy.stats.bernoulli.logpmf(self.bv, mu[k]) 

        log_pmf = np.sum(log_pmf, axis=1)  # Sum of log - like taking product

        # Substract log of denominator - like dividing by it
        for k in range(self.K):
            log_pmf[:, k] -= np.log(denom[k])

        log_resps_numer = log_pi + log_pmf
        log_resps_denom = scipy.special.logsumexp(log_resps_numer, axis=1)
        log_resps = (log_resps_numer.T - log_resps_denom).T
        resps = np.exp(log_resps)

        log_like = np.dot(log_resps_denom, self.read_hist)
        
        return (resps, log_like, denom)


    def maximisation(self, mu, resps, denom):
        """
        Run the Maximization step of the EM algorithm - update mu and pi
        Args:
            X (EM_Class object): Contains the list of bit mut_vectors
            K (int): Number of clusters
            mu (list): DMS reactivities in each cluster
            resps (list): Responsibity of each cluster
            denom (float): Corrected denom of EM algorithm
        Returns:
            mu (list): DMS reactivities in each cluster
            obs_pi (list): Observed proportion of each cluster
            real_pi (list): Corrected proportion of each cluster
        """
        D = self.bv.shape[1]
        mu, obs_pi, real_pi = np.zeros((self.K, D)), np.zeros(self.K), np.zeros(self.K)
        for k in range(self.K):
            N_k = np.sum(resps[:, k] * self.read_hist)
            x_bar_k = np.sum((resps[:, k] * self.read_hist *
                            self.bv.T).T, axis=0) / N_k

            sparse_x_bar_k = np.zeros(self.sparse_mu.shape[1])
            sparse_x_bar_k[self.bases_to_keep] = x_bar_k
            self.sparse_mu[k][self.bases_to_keep] = mu[k]

            upd_mu = newton_krylov(lambda mu_k: mu_der(mu_k, sparse_x_bar_k), self.sparse_mu[k])
            upd_mu = upd_mu[self.bases_to_keep]

            # Check if mu is less than 0. This should not happen.
            zeros = np.where(upd_mu < 0.0)[0]
            if len(zeros) > 0:
                print('Mu is 0 or less:', upd_mu[upd_mu <= 0.0])
                upd_mu = np.maximum(upd_mu, np.zeros(len(upd_mu)))

            mu[k] = upd_mu  # Mu with denom correction
            # mu[k] = x_bar_k  # Mu without denom correction
            obs_pi[k] = N_k / np.sum(self.read_hist)
        real_pi = [obs_pi[k] / denom[k] for k in range(self.K)]
        real_pi = real_pi / np.sum(real_pi)
        
        return (mu, obs_pi, real_pi)


    def run(self):

        """Run the EM clustering algorithm.
        
        Output:
        -------
        A dictionary containing the following attributes:
            mu: array (K x D)
                Mean of each cluster.
            pi: array (K)
                Probability of each cluster.
            log_likelihood: float
                Log-likelihood of the model.
        
        """
        """
        Run the EM algorithm on DMSMaPseq data

        Returns:
            log_like_list (list): Log likelihoods
            final_mu (list): DMS reactivities in each cluster
            final_obs_pi (list): Obs proportion of each cluster
            final_real_pi (list): Corrected proportion of each cluster
            resps (list): Responsibity of each cluster
            BIC (float): BIC value for this K
        """

        N, D = self.bv.shape[0], self.bv.shape[1]

        # ---------------------- Iterations start ---------------------------- #

        # Initialize DMS modification rate for each base in each cluster
        # by sampling from a beta distribution
        BETA_A = 1.5  # Beta dist shape parameter
        BETA_B = 20  # Beta dist shape parameter
        mu = np.asarray([scipy.stats.beta.rvs(BETA_A, BETA_B, size=D)
                        for k in range(self.K)])
        # Initialize cluster probabilties with uniform distribution
        obs_pi = np.asarray([1.0 / self.K] * self.K)

        converged = False
        iter = 1
        log_like_list, mu_list, obs_pi_list, real_pi_list = [], [], [], []

        while not converged:  # Each iteration of the EM algorithm

            if self.verbose:
                print('Iteration:', iter, '|', np.min(mu), np.max(mu))

            # Expectation step
            (resps, log_like, denom) = self.expectation(mu, obs_pi)

            # Maximization step
            (mu, obs_pi, real_pi) = self.maximisation(mu, resps, denom)

            log_like_list.append(log_like)
            mu_list.append(mu)
            obs_pi_list.append(obs_pi)
            real_pi_list.append(real_pi)

            # Check if log like is decreasing - why does this happen?
            if iter > 1 and log_like < log_like_list[-2]:
                prev_loglike = log_like_list[-2]
                print('The log like decreased from {:.9f} to {:.9f}'.format(prev_loglike, log_like))
                if iter >= self.min_iter:
                    converged = True
            else:
                if iter >= self.min_iter:  # At least min iterations has run
                    prev_loglike = log_like_list[-2]
                    diff = log_like - prev_loglike
                    if diff <= self.convergence_cutoff:  # Converged
                        converged = True
                        print('Log like converged after {:d} iterations'.format(iter))

            if iter == self.max_iter:
                converged = True
                print(f'Max iterations of {self.max_iter} reached')
            iter += 1            

        final_mu, final_obs_pi, final_real_pi = mu_list[-1], obs_pi_list[-1], real_pi_list[-1]

        # ------------------------ Iterations end ---------------------------- #

        # BIC = calc_BIC(N, D, self.K, log_like_list[-1]) ## !! Technically it should be the number of non G/T bases, not D !!
        return {'mu': final_mu, 'pi': final_real_pi, 'log_likelihood': log_like_list[-1]}


## ----- Testing if the above code gives same results as original code ----- ##

if __name__ == '__main__' and False:

    exp_path = "/Users/Alberic/Desktop/Pro/RouskinLab/projects/DREEM/"
    bit_Vector_total = np.load(exp_path+"bit_vector.npy")
    read_hist_total =  np.load(exp_path+"read_hist.npy")

    # with open(exp_path+"data_EM_analysis.txt", 'a') as f:
    #     f.write("N_reads max_procs dT mean_memory peak_memory \n")
    #     f.close()
    
    # for N_partial in np.geomspace(10000, 120000, 5).astype(np.int64):
    for N_partial in [10000]:

        N_partial = min(N_partial, bit_Vector_total.shape[0])

        # for max_procs in range(1, 11):
        for max_procs in [2]:

            print("Starting experiment with {} reads and {} cpus".format(N_partial, max_procs))
            
            tracemalloc.start()

            # Get part of the total reads
            bit_Vector = copy.copy(bit_Vector_total[:N_partial])
            read_hist = copy.copy(read_hist_total[:N_partial])
            
            # Run EM and record stats
            EM = EMclustering(bit_Vector, K=2, read_hist=read_hist, min_iter=10, max_procs=max_procs,
                            max_clusters=3, signal_thresh=0.5, info_thresh=0.5, include_g_u=True, include_del=True,
                            min_reads=10,convergence_cutoff=0.5,num_runs=100, verbose=True )
            
            result = EM.run()

            # Log experiment stats
            # memory_stats = [ mem/1e6 for mem in tracemalloc.get_traced_memory() ]
            # log_data = np.array((N_partial,
            #                         max_procs,
            #                         result["dT"],
            #                         memory_stats[0], 
            #                         memory_stats[1] ))
 
            # # with open(exp_path+"data_EM_analysis.txt", 'a') as f:
            #     np.savetxt(f, log_data, newline=" ")
            #     f.write("\n")
            #     f.close()

            # Stop memory tracing and clean memory
            tracemalloc.stop()
            bit_Vector = None
            read_hist = None
            # time.sleep(5)
            
            # Comparing output with previous software -> different test
            mu_reference = np.load("/Users/Alberic/Desktop/Pro/RouskinLab/projects/DREEM/result.npy")
            print("Reference matched:",np.allclose(mu_reference, result["mu"], atol=1e-9))

            for k in range(2):
                plt.subplot(1, 2, k+1)
                plt.plot(result["mu"][k], alpha=0.5)
                plt.plot(mu_reference[k], alpha=0.5)
