import numpy as np
import scipy.stats
from scipy.optimize import newton_krylov

from multiprocessing.dummy import Pool as ThreadPool

import matplotlib.pyplot as plt # !! For testing !!
import time # !! For testing !!

## ------------- Utility functions ------------- ##

def calc_BIC(N, PARAMS_LEN, K, log_like):
    """
    """
    return np.log(N) * PARAMS_LEN * K - (2 * log_like)

def calc_matrixIndices(cpus, N, K):
    """
    Compute the indexes to slice the list of ready in one chunk per processor
    """
    calcsPerCPU = max(round(N * K / cpus), 1)
    inds, start = [], 0
    while start < N:
        coord = (start, start + calcsPerCPU - 1)
        inds.append(coord)
        start = start + calcsPerCPU
    return inds

def calc_denom(i, mu, denom_probs, s2_probs):
    """
    """
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

def logpmf_parallel(bit_vector, mu, ind, k):
    """
    """
    start, end = ind[0], ind[1]
    return scipy.stats.bernoulli.logpmf(bit_vector[start:end + 1], mu[k])


def logpmf_parallel_denom(log_pmf, denom, ind, k):
    """
    """
    start, end = ind[0], ind[1]
    return log_pmf[start:end + 1, k] - np.log(denom[k][0])

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
        Bitvector of the construct. 
    
    K: int
        Number of clusters.
        
    read_hist: array (N)
        Number of reads per bitvector.
        
    min_iter: int
        Minimum number of iterations.
        
    convergence_eps: float
        Convergence threshold. When the difference between the log-likelihood of two consecutive iterations is lower than convergence_eps, the algorithm stops.
        
    Key methods:
    ------------
    
    run: Run the EM clustering algorithm. 
    maximize: Maximization step of the EM algorithm.
    expectation: Expectation step of the EM algorithm.
    
    """
    
    def __init__(self, bv, K, read_hist, n_cpus, max_clusters:int, signal_thresh:float, info_thresh:float, include_g_u:bool, include_del:bool, min_reads:int, min_iter:int, convergence_cutoff:float, num_runs:int, verbose:bool):
        self.bv = bv
        self.K = K
        self.read_hist = read_hist
        self.N = bv.shape[0]
        self.D = bv.shape[1]
        
        self.min_iter = min_iter
        self.cpus = n_cpus

        self.convergence_eps = convergence_cutoff
    
    def expectation(self, mu, pi, calc_inds):
        """
        Run the Expectation step of the EM algorithm - calc log like
        Args:
            X (EM_Class object): Contains the list of bit vectors
            K (int): Number of clusters
            mu (list): DMS reactivities in each cluster
            pi (list): Proportion of each cluster
            calc_inds (list): Indices to split the matrix at
        Returns:
            resps (list): Responsibity of each cluster
            log_like (float): Log likelihood of observing the data
            denom (float): Corrected denom of EM algorithm
        """
        N, D = self.bv.shape[0], self.bv.shape[1]
        log_pi = np.log(pi)
        log_pmf = np.zeros((N, D, self.K))
        denom = [calc_denom(0, mu[k], {}, {}) for k in range(self.K)]

        input_array1 = [[self.bv, mu, ind, k] for ind in calc_inds for k in range(self.K)]
        pool1 = ThreadPool(self.cpus)
        logpmf_results1 = pool1.starmap(logpmf_parallel,
                                        input_array1)
        pool1.close()
        pool1.join()
        for i in range(len(logpmf_results1)):
            ind, k = input_array1[i][2], input_array1[i][3]
            start, end = ind[0], ind[1]
            log_pmf[start:end + 1, :, k] = logpmf_results1[i]

        log_pmf = np.sum(log_pmf, axis=1)  # Sum of log - like taking product

        input_array2 = [[log_pmf, denom, ind, k] for ind in calc_inds
                        for k in range(self.K)]
        pool2 = ThreadPool(self.cpus)
        logpmf_results2 = pool2.starmap(logpmf_parallel_denom,
                                        input_array2)
        pool2.close()
        pool2.join()
        for i in range(len(logpmf_results2)):
            ind, k = input_array2[i][2], input_array2[i][3]
            start, end = ind[0], ind[1]
            log_pmf[start:end + 1, k] = logpmf_results2[i]

        log_resps_numer = np.add(log_pi, log_pmf)
        log_resps_denom = scipy.special.logsumexp(log_resps_numer, axis=1)
        log_resps = np.subtract(log_resps_numer.T, log_resps_denom).T
        resps = np.exp(log_resps)

        log_like = np.dot(log_resps_denom, self.read_hist)
        
        return (resps, log_like, denom)


    def maximisation(self, mu, resps, denom):
        """
        Run the Maximization step of the EM algorithm - update mu and pi
        Args:
            X (EM_Class object): Contains the list of bit vectors
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
            upd_mu = newton_krylov(lambda mu_k: mu_der(mu_k, x_bar_k),
                                mu[k])

            # Check if mu is less than 0. This should not happen.
            zeros = np.where(upd_mu < 0.0)[0]
            if len(zeros) > 0:
                print('Mu is 0 or less:', upd_mu[upd_mu <= 0.0])
                upd_mu = np.maximum(upd_mu, np.zeros(len(upd_mu)))

            mu[k] = upd_mu  # Mu with denom correction
            # mu[k] = x_bar_k  # Mu without denom correction
            obs_pi[k] = N_k / np.sum(self.read_hist)
        real_pi = [obs_pi[k] / denom[k][0] for k in range(self.K)]
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

        # Start and end coordinates for each thread
        calc_inds = calc_matrixIndices(self.cpus, N, self.K)

        # ---------------------- Iterations start ---------------------------- #

        # Initialize DMS modification rate for each base in each cluster
        # by sampling from a beta distribution
        BETA_A = 1.5  # Beta dist shape parameter
        BETA_B = 20  # Beta dist shape parameter
        # np.random.seed(seed=42000) # !! For testing !!
        mu = np.asarray([scipy.stats.beta.rvs(BETA_A, BETA_B, size=D)
                        for k in range(self.K)])

        # Initialize cluster probabilties with uniform distribution
        obs_pi = np.asarray([1.0 / self.K] * self.K)

        converged = False
        iter = 1
        log_like_list, mu_list, obs_pi_list, real_pi_list = [], [], [], []

        dt = []
        while not converged:  # Each iteration of the EM algorithm
            print('Iteration:', iter)
            time_now = time.time()

            # Expectation step
            (resps, log_like, denom) = self.expectation(mu, obs_pi, calc_inds)

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
            else:
                if iter >= self.min_iter:  # At least min iterations has run
                    prev_loglike = log_like_list[-2]
                    diff = log_like - prev_loglike
                    if diff <= self.convergence_eps:  # Converged
                        converged = True
                        print('Log like converged after {:d} iterations'.format(iter))
            iter += 1
            dt.append(time.time()-time_now)
        final_mu, final_obs_pi, final_real_pi = mu_list[-1], obs_pi_list[-1], real_pi_list[-1]

        print("Average dT", np.mean(dt))
        plt.hist(dt, bins=30)
        plt.show()


        # ------------------------ Iterations end ---------------------------- #

        # BIC = calc_BIC(N, D, self.K, log_like_list[-1]) ## !! Technically it should be the number of non G/T bases, not D !!

        return {'mu': final_mu, 'pi': final_real_pi, 'log_likelihood': log_like_list[-1]}


## ----- Testing if the above code gives same results as original code ----- ##

if False:
    bit_Vector = np.zeros((3,3)) # np.load("/Users/Alberic/Desktop/Pro/RouskinLab/projects/DREEM/bit_vector.npy")
    read_hist =  np.zeros((3,3)) # np.load("/Users/Alberic/Desktop/Pro/RouskinLab/projects/DREEM/read_hist.npy")
    EM = EMclustering(bit_Vector, 2, read_hist, min_iter=10)

    result = EM.run()


    mu_reference = np.load("/Users/Alberic/Desktop/Pro/RouskinLab/projects/DREEM/result.npy")

    print("Reference matched:",(mu_reference == result["mu"]).all())



    for k in range(2):
        plt.subplot(1, 2, k+1)
        plt.plot(result["mu"][k])

    plt.show()