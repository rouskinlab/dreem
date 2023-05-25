from functools import cached_property
from logging import getLogger

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import dirichlet

from ..call.load import BitVecLoader
from dreem.core.mu import unbias, denom
from dreem.core.bit import UniqMutBits

logger = getLogger(__name__)


def calc_bic(n_params: int, n_data: int, log_like: float, factor: float = 10.):
    """
    Compute the Bayesian Information Criterion (BIC) of a model.
    Typically, the model with the smallest BIC is preferred.

    Parameters
    ----------
    n_params: int
        Number of parameters that the model estimates
    n_data: int
        Number of data points in the sample from which the parameters
        were estimated
    log_like: float
        Natural logarithm of the likelihood of observing the data given
        the parameters
    factor: float = 10.0
        In order for the BIC approximation to be valid, the sample size
        must be much larger than the number of estimated parameters.
        Issue a warning if the sample size is less than factor times the
        number of parameters, but still compute and return the BIC.

    Returns
    -------
    float
        The Bayesian Information Criterion (BIC)
    """
    if n_data < factor * n_params:
        logger.warning(f"The BIC approximation is valid only when the sample "
                       f"size (n = {n_data}) is much larger than the number "
                       f"of parameters being estimated (p = {n_params}). "
                       f"This model does not meet this criterion, so the BIC "
                       f"may not indicate the model's complexity accurately.")
    return n_params * np.log(n_data) - 2. * log_like


class EmClustering(object):
    """ Run expectation-maximization to cluster the given bit vectors
    into the specified number of clusters."""

    def __init__(self,
                 loader: BitVecLoader,
                 muts: UniqMutBits,
                 n_clusters: int,
                 i_run: int, *,
                 min_iter: int,
                 max_iter: int,
                 conv_thresh: float):
        """
        Parameters
        ----------
        loader: BitVecLoader
            Loader of the filtered bit vectors
        muts: UniqMutBits
            Container of unique bit vectors of mutations
        n_clusters: int
            Number of clusters into which to cluster the bit vectors;
            must be a positive integer
        i_run: int
            Run number; must be a positive integer
        min_iter: int
            Minimum number of iterations for clustering. Must be a
            positive integer no greater than max_iter.
        max_iter: int
            Maximum number of iterations for clustering. Must be a
            positive integer no less than min_iter.
        conv_thresh: float
            Stop the algorithm when the difference in log likelihood
            between two successive iterations becomes smaller than the
            convergence threshold (and at least min_iter iterations have
            run). Must be a positive real number (ideally close to 0).
        """
        # Filter loader
        self.loader = loader
        # Unique bit vectors of mutations
        self.muts = muts
        # Number of clusters
        if not n_clusters >= 1:
            raise ValueError(f"n_clusters must be ≥ 1, but got {n_clusters}")
        self.ncls = n_clusters
        # Run number
        if not i_run >= 1:
            raise ValueError(f"i_run must be ≥ 1, but got {i_run}")
        self.i_run = i_run
        # Minimum number of iterations of EM
        if not min_iter >= 1:
            raise ValueError(f"min_iter must be ≥ 1, but got {min_iter}")
        self.min_iter = min_iter
        # Maximum number of iterations of EM
        if not max_iter >= min_iter:
            raise ValueError(f"max_iter must be ≥ min_iter ({min_iter}), "
                             f"but got {max_iter}")
        self.max_iter = max_iter
        # Cutoff for convergence of EM
        if not conv_thresh >= 0.:
            raise ValueError(f"conv_thresh must be ≥ 0, but got {conv_thresh}")
        self.conv_thresh = conv_thresh
        # Number of reads observed in each cluster (no bias correction)
        self.nreads = np.empty(self.ncls, dtype=float)
        # Log of the denominator of each cluster
        self.log_denom = np.empty(self.ncls, dtype=float)
        # Mutation rates of used positions (row) in each cluster (col)
        self.mus = np.empty((self.loader.pos_kept.size, self.ncls), dtype=float)
        # Positions of the section that will be used for clustering
        # (0-indexed from the beginning of the section)
        self.sparse_pos = self.loader.pos_kept - self.loader.end5
        # Mutation rates of all positions, including those not used for
        # clustering (row), in each cluster (col). The rate for every
        # unused position always remains zero.
        self.sparse_mus = np.zeros((self.loader.n_pos_init, self.ncls),
                                   dtype=float)
        # Likelihood of each vector (col) coming from each cluster (row)
        self.resps = np.empty((self.ncls, self.muts.n_uniq), dtype=float)
        # Trajectory of log likelihood values.
        self.log_likes: list[float] = list()
        # Number of iterations.
        self.iter = 0
        # Whether the algorithm has converged.
        self.converged = False

    def update_sparse_mus(self):
        """ Update the sparse mutation rates using the current mutation
        rates. """
        # Copy mutation rates to the rows of sparse_mu that correspond
        # to used positions (unused rows remain zero).
        self.sparse_mus[self.sparse_pos] = self.mus
        return self.sparse_mus

    @cached_property
    def cluster_nums(self):
        """ Return an array of the cluster numbers, starting from 1. """
        return np.arange(1, self.ncls + 1, dtype=int)

    @property
    def prop_obs(self):
        """ Return the observed proportion of each cluster, without
        correcting for the drop-out bias. """
        return self.nreads / np.sum(self.nreads)

    @property
    def prop_real(self):
        """ Calculate the real proportion of each cluster, corrected for
        the drop-out bias. """
        # Correct the bias in the proportions of reads in each cluster
        # by re-weighting them by the reciprocal of the denominators.
        weighted_prop_obs = self.prop_obs / np.exp(self.log_denom)
        return weighted_prop_obs / np.sum(weighted_prop_obs)

    @property
    def log_like(self):
        """ Return the current log likelihood, which is the last item in
        the trajectory of log likelihood values. """
        try:
            return self.log_likes[-1]
        except IndexError:
            # No log likelihood values have been computed.
            return np.nan

    @property
    def log_like_prev(self):
        """ Return the previous log likelihood, which is the penultimate
        item in the trajectory of log likelihood values. """
        try:
            return self.log_likes[-2]
        except IndexError:
            # Fewer than two log likelihood values have been computed.
            return np.nan

    @property
    def delta_log_like(self):
        """ Compute the change in log likelihood from the previous to
        the current iteration. """
        return self.log_like - self.log_like_prev

    @property
    def bic(self):
        """ Compute this model's Bayesian Information Criterion. """
        # The parameters estimated by the model are the mutation rates
        # for each position in each cluster and the proportion of each
        # cluster in the population. By contrast, the cluster membership
        # values are latent variables because each describes one item in
        # the sample (a bit vector), not a parameter of the population.
        # The number of data points is the number of unique bit vectors.
        return calc_bic(self.mus.size + self.prop_real.size,
                        self.muts.n_uniq,
                        self.log_like)

    def _max_step(self):
        """ Run the Maximization step of the EM algorithm. """
        # Calculate the number of reads in each cluster by summing the
        # count-weighted likelihood that each bit vector came from the
        # cluster.
        self.nreads = self.resps @ self.muts.counts
        # logger.debug(f"NREADS:\n{self.nreads}")
        # Copy the sparse mutation rates from the previous iteration.
        sparse_mus_prev = self.sparse_mus.copy()
        # Compute the observed mutation rate at each position (j).
        for j, muts_j in enumerate(self.muts.indexes):
            # Calculate the number of mutations at each position in each
            # cluster by summing the count-weighted likelihood that each
            # bit vector with a mutation at (j) came from the cluster,
            # then divide by the count-weighted sum of the number of
            # reads in the cluster to find the observed mutation rate.
            self.mus[j] = (self.resps[:, muts_j]
                           @ self.muts.counts[muts_j]) / self.nreads
        # Solve for the real mutation rates that are expected to yield
        # the observed mutation rates after considering read drop-out.
        # Constrain the mutation rates to [min_mu, max_mu].
        self.mus = unbias(self.update_sparse_mus(),
                          self.loader.min_mut_gap,
                          sparse_mus_prev)[self.sparse_pos]

    def _exp_step(self):
        """ Run the Expectation step of the EM algorithm. """
        # SUB-STEP E1: Update the denominator of each cluster.
        # Update the denominator of each cluster using the sparse mus.
        self.log_denom = np.log(denom(self.update_sparse_mus(),
                                      self.loader.min_mut_gap))

        # SUB-STEP E2: Compute for each cluster and observed bit vector
        #              the joint probability that another random bit
        #              vector would both come from the same cluster and
        #              have exactly the same set of bits.
        # Compute the logs of the mutation and non-mutation rates.
        with np.errstate(divide="ignore"):
            # Suppress warnings about taking the log of zero, which is a
            # valid mutation rate. Thus, here we CAN compute log of 0,
            # just like Chuck Norris.
            log_mus = np.log(self.mus)
        log_nos = np.log(1. - self.mus)
        # To save memory at no cost, instead of allocating a new array
        # for the joint probabilities, write them into self.resps, which
        # has the correct shape and will anyway not be needed before the
        # next time its values are updated.
        # Loop over each cluster (k).
        for k in range(self.ncls):
            # Compute the probability that a bit vector has no mutations
            # given that it comes from cluster (k), which is the sum of
            # all not-mutated log probabilities (sum(log_nos[:, k]))
            # minus the log denominator of the cluster (log_denom[k]).
            log_prob_no_muts_given_k = np.sum(log_nos[:, k]) - self.log_denom[k]
            # Initialize the log probability for all bit vectors in
            # cluster (k) to the log probability that an observed bit
            # vector comes from cluster (k) (log(prob_obs[k])) and has
            # no mutations given that it came from cluster (k).
            log_prob_no_muts_and_k = (np.log(self.prop_obs[k])
                                      + log_prob_no_muts_given_k)
            self.resps[k].fill(log_prob_no_muts_and_k)
            # Loop over each position (j); each iteration adds the log
            # PMF for one additional position in each bit vector to the
            # accumulating total log probability of each bit vector.
            # Using accumulation with this loop also uses less memory
            # than would holding the PMF for every position in an array
            # and then summing over the position axis.
            for j in range(self.loader.pos_kept.size):
                # Compute the probability that each bit vector would
                # have the bit observed at position (j). The probability
                # is modeled using a Bernoulli distribution, with PMF:
                # log_mus[k, j] if muts[j, i] else log_nos[k, j].
                # This PMF could be computed explicitly, such as with
                # scipy.stats.bernoulli.logpmf(muts[j], mus[k, j]).
                # But few of the bits are mutated, so a more efficient
                # (read: much, MUCH faster) way is to assume initially
                # that no bit is mutated -- that is, to initialize the
                # probability of every bit vector with the probability
                # that the bit vector has no mutations, as was done by
                # prop_obs[k] + (np.sum(log_nos[k]) - log_denom[k]).
                # Then, for only the few bits that are mutated, the
                # probability is adjusted by adding the log of the
                # probability that the base is mutated minus the log of
                # the probability that the base is not mutated.
                adj_for_mut_j_given_k = log_mus[j, k] - log_nos[j, k]
                self.resps[k][self.muts.indexes[j]] += adj_for_mut_j_given_k
        # logger.debug(f"Computed joint log probs of {self}:\n{self.resps}")

        # SUB-STEP E3: Compute the marginal probability of observing
        #              each bit vector, and use them to compute the
        #              posterior probabilities that each bit vector came
        #              from each cluster, and the total log likelihood.
        # For each observed bit vector, the probability that a random
        # generator would yield the same series of bits (regardless of
        # which cluster the random generator selects) is the sum over
        # all clusters of the joint probability of selecting the cluster
        # and generating the bit vector from the cluster.
        marg_log_prob = logsumexp(self.resps, axis=0)
        # logger.debug(f"Computed marginal log probs of {self}:\n{marg_log_prob}")
        # Calculate the posterior probability that each bit vector came
        # from each cluster by dividing the joint probability (observing
        # the bit vector and coming from the cluster) by the marginal
        # probability (observing the bit vector in any cluster).
        self.resps = np.exp(self.resps - marg_log_prob)
        # logger.debug(f"Computed responsibilities of {self}:\n{self.resps}")
        # Calculate the log likelihood of observing all the bit vectors
        # by summing the log probability over all bit vectors, weighted
        # by the number of times each bit vector occurs. Cast to a float
        # explicitly to verify that the product is a scalar.
        log_like = float(np.vdot(marg_log_prob, self.muts.counts))
        # logger.debug(f"Computed log likelihood of {self}:\n{log_like}")
        self.log_likes.append(log_like)

    def run(self):
        """ Run the EM clustering algorithm. """
        # Erase the trajectory of log likelihood values (if any).
        self.log_likes.clear()
        # Choose the concentration parameters using a standard uniform
        # distribution so that the reads assigned to each cluster can
        # vary widely among the clusters (helping to make the clusters
        # distinct from the beginning) and among different runs of the
        # algorithm (helping to start the runs in very different states
        # so that they can explore much of the parameter space).
        # Use the half-open interval (0, 1] because the concentration
        # parameter of a Dirichlet distribution can be 1 but not 0.
        conc_params = 1. - np.random.default_rng().random(self.ncls)
        # Initialize cluster membership with a Dirichlet distribution.
        self.resps = dirichlet.rvs(conc_params, self.muts.n_uniq).T
        # Run EM until the log likelihood converges or the number of
        # iterations reaches max_iter, whichever happens first.
        self.converged = False
        for self.iter in range(1, self.max_iter + 1):
            # Update the mutation rates and cluster proportions.
            self._max_step()
            # Update the cluster membership and log likelihood.
            self._exp_step()
            if not np.isfinite(self.log_like):
                raise ValueError(f"{self}, iteration {self.iter} returned a "
                                 f"non-finite log likelihood: {self.log_like}")
            logger.debug(f"{self}, iteration {self.iter}: "
                         f"log likelihood = {self.log_like}")
            # Check for convergence.
            if self.delta_log_like < 0.:
                # The log likelihood should not decrease.
                logger.warning(f"{self}, iteration {self.iter} returned a "
                               f"smaller log likelihood ({self.log_like}) than "
                               f"the previous iteration ({self.log_like_prev})")
            elif (self.delta_log_like < self.conv_thresh
                  and self.iter >= self.min_iter):
                # Converge if the increase in log likelihood is
                # smaller than the convergence cutoff and at least
                # the minimum number of iterations have been run.
                self.converged = True
                logger.info(f"{self} converged on iteration {self.iter}: "
                            f"last log likelihood = {self.log_like}")
                break
        else:
            # The log likelihood did not converge within the maximum
            # number of iterations.
            logger.warning(f"{self} failed to converge within {self.max_iter} "
                           f"iterations: last log likelihood = {self.log_like}")
        # Return this instance so that any code that runs multiple
        # EmClustering objects can create and run them in one line.
        return self

    @staticmethod
    def props_index():
        return pd.Index(["Real", "Observed"])

    def output_props(self):
        """ Return a DataFrame of the real and observed proportion of
        each cluster. """
        return pd.DataFrame(np.vstack([self.prop_real, self.prop_obs]),
                            index=self.props_index(),
                            columns=self.cluster_nums)

    def output_mus(self):
        """ Return a DataFrame of the mutation rate at each position
        for each cluster. """
        return pd.DataFrame(self.mus,
                            index=self.loader.index_kept,
                            columns=self.cluster_nums)

    def output_resps(self):
        """ Return the responsibilities of the reads. """
        return pd.DataFrame(self.resps.T[self.muts.inverse],
                            index=self.loader.get_read_names(),
                            columns=self.cluster_nums)

    def __str__(self):
        return f"Clustering {self.loader}, k = {self.ncls}, run {self.i_run}"
