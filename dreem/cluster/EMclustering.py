import numpy as np

class EMclustering:
    """This class runs the EM clustering algorithm.
    
    Parameters from args:
    ---------------------
    
    bv: array (N x D)
        Bitvector of the construct. 
    
    K: int
        Number of clusters.
        
    reads_hist: array (N)
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
    
    def __init__(self, bv, K, reads_hist, min_iter=100, convergence_eps=0.5) -> None:
        self.bv = bv
        self.K = K
        self.reads_hist = reads_hist
        self.N = bv.shape[0]
        self.D = bv.shape[1]
    
    def maximize(self, gamma):
        """Maximization step of the EM algorithm.
        
        Parameters:
        -----------
        
        gamma: array (N x K)
            Posterior probability of each bitvector to belong to each cluster.
            
        Output:
        -------
        
        mu: array (K x D)
            Mean of each cluster.
            
        pi: array (K)
            Probability of each cluster.
            
        """
        
        mu = np.zeros((self.K, self.D))
        pi = np.zeros(self.K)
        
        for k in range(self.K):
            # TODO: write the likelihood function
            mu[k] = np.sum(self.bv * gamma[:, k].reshape(-1, 1), axis=0) / np.sum(gamma[:, k])
            pi[k] = np.sum(gamma[:, k]) / self.N
            
        return mu, pi
    

    def expectation(self, mu, pi):
        """Expectation step of the EM algorithm.
        
        Parameters:
        -----------
        
        mu: array (K x D)
            Mean of each cluster.
            
        pi: array (K)
            Probability of each cluster.
            
        Output:
        -------
        
        gamma: array (N x K)
            Posterior probability of each bitvector to belong to each cluster.
            
        """
        
        gamma = np.zeros((self.N, self.K))
        
        for k in range(self.K):
            gamma[:, k] = pi[k] * np.prod(mu[k] ** self.bv * (1 - mu[k]) ** (1 - self.bv), axis=1)
            
        gamma = gamma / np.sum(gamma, axis=1).reshape(-1, 1)
        
        return gamma

    def run(self):
        """Run the EM clustering algorithm.
        """
            
        # Initialization
        mu = np.random.rand(self.K, self.D)
        pi = np.random.rand(self.K)
        pi = pi / np.sum(pi)
        
        # EM algorithm
        i = 0
        while True:
            gamma = self.expectation(mu, pi)
            mu, pi = self.maximize(gamma)