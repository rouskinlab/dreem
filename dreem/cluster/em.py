from logging import getLogger

import numpy as np
from scipy.optimize import newton_krylov
from scipy.special import logsumexp
from scipy.stats import bernoulli, beta

from .bitvector import BitVector

logger = getLogger(__name__)


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
        return denom_probs[i], s2_probs
    elif i >= len(mu):  # Base case
        return 1, s2_probs
    else:  # Make the calc
        s1 = calc_denom(i + 1, mu, denom_probs, s2_probs)[0]
        s2 = ((1.0 - mu[i + 1: i + 4]).prod() *
              calc_denom(i + 4, mu, denom_probs, s2_probs)[0])
        denom_probs[i] = ((1 - mu[i]) * s1) + (mu[i] * s2)
        s2_probs[i] = s2
        return denom_probs[i], s2_probs


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
    """
    This class runs the EM clustering algorithm.
    
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

    def __init__(self, bv: BitVector,
                 n_clusters: int,
                 min_iter: int,
                 max_iter: int,
                 convergence_cutoff: float):
        self.bv = bv
        self.bases_to_keep = bv.positions - bv.end5  # 0-indexed
        self.nclust = n_clusters
        self.nvec = bv.n_reads_used
        self.npos = bv.n_pos_used
        self.seqlen = len(bv.seq)
        self.min_iter = min_iter
        self.max_iter = max_iter
        self.convergence_cutoff = convergence_cutoff
        self.sparse_mu = np.zeros((n_clusters, self.seqlen), dtype=float)

    def expectation(self, mu, pi):
        """
        Run the Expectation step of the EM algorithm - calc log like
        Args:
            mu (list): DMS reactivities in each cluster
            pi (list): Proportion of each cluster
            
        Returns:
            resps (list): Responsibity of each cluster
            log_like (float): Log likelihood of observing the data
            denom (float): Corrected denom of EM algorithm
        """
        log_pi = np.log(pi)
        log_pmf = np.zeros((self.nvec, self.npos, self.nclust), dtype=float)
        for k in range(self.nclust):
            self.sparse_mu[k][self.bases_to_keep] = mu[k]
        denom = [calc_denom(0, self.sparse_mu[k], {}, {})[0] for k in range(self.nclust)]

        # Compute probability mass function for all reads, for each cluster
        for k in range(self.nclust):
            log_pmf[:, :, k] = bernoulli.logpmf(self.bv.muts, mu[k])

        log_pmf = np.sum(log_pmf, axis=1)  # Sum of log - like taking product

        # Substract log of denominator - like dividing by it
        for k in range(self.nclust):
            log_pmf[:, k] -= np.log(denom[k])

        log_resps_numer = log_pi + log_pmf
        log_resps_denom = logsumexp(log_resps_numer, axis=1)
        log_resps = (log_resps_numer.T - log_resps_denom).T
        resps = np.exp(log_resps)

        log_like = np.dot(log_resps_denom, self.bv.read_hist)

        return resps, log_like, denom

    def maximisation(self, mu, resps, denom):
        """
        Run the Maximization step of the EM algorithm - update mu and pi
        Args:
            mu (list): DMS reactivities in each cluster
            resps (list): Responsibity of each cluster
            denom (float): Corrected denom of EM algorithm
        Returns:
            mu (list): DMS reactivities in each cluster
            obs_pi (list): Observed proportion of each cluster
            real_pi (list): Corrected proportion of each cluster
        """
        mu, obs_pi, real_pi = np.zeros((self.nclust, self.npos)), np.zeros(self.nclust), np.zeros(self.nclust)
        for k in range(self.nclust):
            N_k = np.sum(resps[:, k] * self.bv.read_hist)
            x_bar_k = np.sum((resps[:, k] * self.bv.read_hist *
                              self.bv.muts.T).T, axis=0) / N_k

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
            obs_pi[k] = N_k / np.sum(self.bv.read_hist)
        real_pi = [obs_pi[k] / denom[k] for k in range(self.nclust)]
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

        # ---------------------- Iterations start ---------------------------- #

        # Initialize DMS modification rate for each base in each cluster
        # by sampling from a beta distribution
        BETA_A = 1.5  # Beta dist shape parameter
        BETA_B = 20  # Beta dist shape parameter
        mu = beta.rvs(BETA_A, BETA_B, size=(self.nclust, self.npos))
        # Initialize cluster probabilties with uniform distribution
        obs_pi = np.full(self.nclust, 1.0 / self.nclust)

        converged = False
        iter = 1
        log_like_list, mu_list, obs_pi_list, real_pi_list = [], [], [], []

        while not converged:  # Each iteration of the EM algorithm

            logger.debug(f"Iteration {iter}")

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
                logger.warning('The log like decreased from {:.9f} to {:.9f}'.format(prev_loglike, log_like))
                if iter >= self.min_iter:
                    converged = True
            else:
                if iter >= self.min_iter:  # At least min iterations has run
                    prev_loglike = log_like_list[-2]
                    diff = log_like - prev_loglike
                    if diff <= self.convergence_cutoff:  # Converged
                        converged = True
                        logger.info('Log like converged after {:d} iterations'.format(iter))

            if iter == self.max_iter:
                converged = True
                logger.warning(f'Max iterations of {self.max_iter} reached')
            iter += 1

        final_mu, final_obs_pi, final_real_pi = mu_list[-1], obs_pi_list[-1], real_pi_list[-1]

        # ------------------------ Iterations end ---------------------------- #

        # BIC = calc_BIC(N, D, self.K, log_like_list[-1]) ## !! Technically it should be the number of non G/T bases, not D !!
        return {'mu': final_mu, 'pi': final_real_pi, 'log_likelihood': log_like_list[-1]}
