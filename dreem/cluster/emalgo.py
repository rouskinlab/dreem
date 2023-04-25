from functools import cached_property
from logging import getLogger

import numpy as np
import pandas as pd
from scipy.optimize import newton_krylov
from scipy.special import logsumexp
from scipy.stats import bernoulli, dirichlet

from .bvec import BitVector

logger = getLogger(__name__)


def calc_bic(n_params: int, n_data: int, log_like: float, factor: float = 10.0):
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
        Issue a warning if the samples size is less than factor times
        the number of parameters, but still compute the BIC.

    Returns
    -------
    float
        The Bayesian Information Criterion (BIC)
    """
    if n_data < factor * n_params:
        logger.warning(f"The BIC approximation is valid only when the sample "
                       f"size (n = {n_data}) is much larger than the number "
                       f"of parameters being estimated (p = {n_params}). "
                       f"This model does not meet this criterion.")
    return n_params * np.log(n_data) - 2.0 * log_like


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


class EmClustering(object):
    def __init__(self,
                 bvec: BitVector,
                 n_clusters: int,
                 replicate: int, *,
                 min_iter: int,
                 max_iter: int,
                 conv_thresh: float,
                 epsilon: float = 1e-6):
        """
        Run expectation-maximization to cluster the given bit vectors
        into the specified number of clusters.

        Parameters
        ----------
        bvec: BitVector
            Bit vector container
        n_clusters: int
            Number of clusters into which to cluster the bit vectors;
            must be a positive integer
        replicate: int
            Replicate number; must be a positive integer
        min_iter: int
            Minimum number of iterations for clustering; must be a
            positive integer no greater than max_iter
        max_iter: int
            Maximum number of iterations for clustering; must be a
            positive integer no less than min_iter
        conv_thresh: float
            Stop the algorithm when the difference in log likelihood
            between two successive iterations becomes smaller than the
            convergence threshold (and at least min_iter iterations have
            run); must be a (small) positive real number
        epsilon: float = 10^-6
            Keep all mutation rates within [epsilon, 1 - epsilon]; if
            any stray outside this range, then issue a warning and move
            them back into the range; must be a (small) positive real
            number no greater than 0.5
        """
        # Container of bit vectors
        self.bvec = bvec
        # Number of clusters
        if not n_clusters >= 1:
            raise ValueError(f"n_clusters must be ≥ 1, but got {n_clusters}")
        self.ncls = n_clusters
        # Replicate
        if not replicate >= 1:
            raise ValueError(f"replicate must be ≥ 1, but got {replicate}")
        self.rep = replicate
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
        if not conv_thresh >= 0.0:
            raise ValueError(f"conv_thresh must be ≥ 0, but got {conv_thresh}")
        self.conv_thresh = conv_thresh
        # Minimum and maximum permitted mutation rates
        if not 0.0 < epsilon <= 0.5:
            raise ValueError(f"epsilon must be in (0, 0.5], but got {epsilon}")
        self.min_mu = epsilon
        self.max_mu = 1.0 - epsilon
        # Number of reads observed in each cluster (no bias correction)
        self.nreads = np.empty(self.ncls, dtype=float)
        # Denominator of each cluster
        self.denom = np.empty(self.ncls, dtype=float)
        # Mutation rate of each position (col) in each cluster (row)
        self.mus = np.empty((self.ncls, self.bvec.n_pos_use), dtype=float)
        # Likelihood of each vector (col) coming from each cluster (row)
        self.resps = np.empty((self.ncls, self.bvec.n_bvec_use), dtype=float)
        # Trajectory of log likelihood values.
        self.log_likes: list[float] = list()
        # Number of iterations.
        self.iter = 0
        # Whether the algorithm has converged.
        self.converged = False

    @cached_property
    def cluster_nums(self):
        """ Return an array of the cluster numbers, starting from 1. """
        return np.arange(1, self.ncls + 1, dtype=int)

    @cached_property
    def sparse_pos(self):
        """ Positions of the section that will be used for clustering
        (0-indexed from the beginning of the section). """
        return self.bvec.positions - self.bvec.end5

    @property
    def prop_obs(self):
        """ Calculate the observed proportion of each cluster, without
        correcting for the drop-out bias. """
        return self.nreads / self.bvec.n_read_use

    @property
    def prop_real(self):
        """ Calculate the real proportion of each cluster, corrected for
        the drop-out bias. """
        # Correct the bias in the proportions of reads in each cluster
        # by re-weighting them by the reciprocal of the denominators.
        weighted = self.prop_obs / self.denom
        return weighted / np.sum(weighted)

    @property
    def sparse_mus(self):
        """ Mutation rate of each position (col) in each cluster (row),
        including all positions that are not used for clustering. """
        sparse_mus = np.zeros((self.ncls, len(self.bvec.seq)), dtype=float)
        # Copy all mutation rates to the columns of sparse_mu that
        # correspond to used positions (unused columns remain zero).
        sparse_mus[:, self.sparse_pos] = self.mus
        return sparse_mus

    @property
    def log_like(self):
        """ Return the current log likelihood, which is the last item in
        the trajectory of log likelihood values. """
        try:
            return self.log_likes[-1]
        except IndexError:
            # No log likelihood values have been computed.
            return float("nan")

    @property
    def delta_log_like(self):
        """ Compute the change in log likelihood from the previous to
        the current iteration. """
        try:
            return self.log_like - self.log_likes[-2]
        except IndexError:
            # Fewer than two log likelihood values have been computed.
            return float("nan")

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
                        self.bvec.n_bvec_use,
                        self.log_like)

    def _max_step(self):
        """ Run the Maximization step of the EM algorithm. """
        # SUB-STEP M1: Calculate mutation rates without bias correction.
        # Loop over each cluster (k).
        for k, resps_k in enumerate(self.resps):
            # Calculate the number of reads in the cluster, which is the
            # sum over all bit vectors of the likelihood that the vector
            # came from the cluster, weighted by the number of times the
            # bit vector occurs.
            self.nreads[k] = np.vdot(resps_k, self.bvec.counts)
            # Loop over each position (j).
            for j, muts_j in enumerate(self.bvec.bvec):
                # Calculate the number of mutations at the position by
                # summing the likelihoods of the bit vectors weighted by
                # the count of each bit vector, similar to nreads[k] but
                # for only bit vectors with mutations at this position.
                n_muts_kj = np.vdot(resps_k[muts_j], self.bvec.counts[muts_j])
                # Compute the mutation rate without bias correction.
                self.mus[k, j] = n_muts_kj / self.nreads[k]
        logger.debug(f"Computed uncorrected mus of {self}:\n{self.mus}")

        # SUB-STEP M2: Correct the bias in mutation rates.
        # Loop over every cluster (k).
        for k, sparse_mus_k in enumerate(self.sparse_mus):
            # Solve for the corrected mutation rates for the cluster
            # using Newton's method with the Krylov approximation.
            # Compute the derivative of the error using mu_der, which
            # compares the current solution for mus[k] (x) with the
            # uncorrected mutation rates (sparse_mus_k). The initial
            # guess for the solution is also sparse_mus_k, since it
            # is usually close to the optimal solution.
            sol = newton_krylov(lambda x: mu_der(x, sparse_mus_k), sparse_mus_k)
            # Copy the mutation rates of the positions that are used
            # from the sparse solution (sol) to the mutation rates for
            # the current cluster.
            self.mus[k] = sol[self.sparse_pos]
        logger.debug(f"Computed corrected mus of {self}:\n{self.mus}")

        # SUB-STEP M3: Verify that mutation rates are within limits.
        if n_low := np.sum(np.less(self.mus, self.min_mu)):
            # If any mutation rates have become too small, issue a
            # warning and reset those mutation rates to the minimum.
            logger.warning(f"Mutation rates of {n_low} positions have become "
                           f"< {self.min_mu}; resetting to {self.min_mu}")
            self.mus[:, :] = np.maximum(self.mus, self.min_mu)
        if n_high := np.sum(np.greater(self.mus, self.max_mu)):
            # If any mutation rates have become too large, issue a
            # warning and reset those mutation rates to the maximum.
            logger.warning(f"Mutation rates of {n_high} positions have become "
                           f"> {self.max_mu}; resetting to {self.max_mu}")
            self.mus = np.minimum(self.mus, self.max_mu)

    def _exp_step(self):
        """ Run the Expectation step of the EM algorithm. """
        logger.debug(f"Began expectation step {self.iter} of {self}")
        # SUB-STEP E1: Update the denominator of each cluster.
        for k, sparse_mus_k in enumerate(self.sparse_mus):
            # Calculate the denominator using the sparse mutation rates.
            self.denom[k] = calc_denom(0, sparse_mus_k, {}, {})[0]
        logger.debug(f"Computed denominators of {self}:\n{self.denom}")

        # SUB-STEP E2: Compute for each cluster and observed bit vector
        #              the joint probability that another random bit
        #              vector would both come from the same cluster and
        #              have exactly the same set of bits.
        # For each cluster, the prior probability that a randomly chosen
        # bit vector comes from the cluster, given that the bit vector
        # is observed (i.e. present in the data set, so has no pair of
        # mutations too close), is the cluster's real proportion.
        log_priors = np.log(self.prop_obs)
        logger.debug(f"Computed log priors of {self}:\n{log_priors}")
        # For each cluster, the probability of generating a specific bit
        # vector would be the product of the Bernoulli probability mass
        # function (PMF) over all the positions in the bit vector if no
        # bit vectors dropped out. Because of drop-out, for each bit
        # vector with mutations too close, the probability of observing
        # the bit vector is zero. So that the probabilities of the bit
        # vectors remaining sum to unity, they need to be inflated by
        # dividing them by the cluster's denominator (which is ≤ 1).
        # A computationally efficient method to divide all the Bernoulli
        # PMF values by the denominator is simply to divide the prior
        # probability for each cluster by the denominator, since the
        # denominator is a property of the cluster, not the bit vector.
        log_priors -= np.log(self.denom)
        logger.debug(f"Adjusted log priors of {self}:\n{log_priors}")
        # To save memory at no cost, instead of allocating a new array
        # for the joint probabilities, write them into self.resps, which
        # has the correct shape and will anyway not be needed before the
        # next time its values are updated.
        # Loop over each cluster (k).
        for resps_k, mus_k, log_prior_k in zip(self.resps, self.mus, log_priors,
                                               strict=True):
            # Initialize the probability for all bit vectors in cluster
            # (k) to the prior probability that was previously computed.
            resps_k.fill(log_prior_k)
            # Loop over each position (j); each iteration adds the log
            # PMF for one additional position in each bit vector to the
            # accumulating total log probability of each bit vector.
            # Using accumulation with this loop also uses less memory
            # than would holding the PMF for every position in an array
            # and then summing over the position axis.
            for muts_j, mu_kj in zip(self.bvec.bvec, mus_k, strict=True):
                # Compute the Bernoulli probability mass function for
                # each vector: the log probability that bit vector (i)
                # would come from cluster (k) and have the bit observed
                # at position (j). Recall that these log probabilities
                # have already been normalized by subtracting the log of
                # the cluster's denominator.
                resps_k += bernoulli.logpmf(muts_j, mu_kj)
        logger.debug(f"Computed joint log probs of {self}:\n{self.resps}")

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
        logger.debug(f"Computed marginal log probs of {self}:\n{marg_log_prob}")
        # Calculate the posterior probability that each bit vector came
        # from each cluster by dividing the joint probability (observing
        # the bit vector and coming from the cluster) by the marginal
        # probability (observing the bit vector in any cluster).
        self.resps = np.exp(self.resps - marg_log_prob)
        logger.debug(f"Computed responsibilities of {self}:\n{self.resps}")
        # Calculate the log likelihood of observing all the bit vectors
        # by summing the log probability over all bit vectors, weighted
        # by the number of times each bit vector occurs. Cast to a float
        # explicitly to verify that the product is a scalar.
        log_like = float(np.vdot(marg_log_prob, self.bvec.counts))
        logger.debug(f"Computed log likelihood of {self}:\n{log_like}")
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
        conc_params = 1.0 - np.random.default_rng().random(self.ncls)
        # Initialize cluster membership with a Dirichlet distribution.
        self.resps[:, :] = dirichlet.rvs(conc_params, self.bvec.n_bvec_use).T
        # Run EM until the log likelihood converges or the number of
        # iterations reaches max_iter, whichever happens first.
        self.converged = False
        for self.iter in range(1, self.max_iter + 1):
            logger.debug(f"Began iteration {self.iter} of {self}")
            # Update the mutation rates and cluster proportions.
            self._max_step()
            # Update the cluster membership and log likelihood.
            self._exp_step()
            logger.info(f"Ended iteration {self.iter} of {self} "
                        f"(log likelihood = {self.log_like})")
            if not np.isfinite(self.log_like):
                raise ValueError(f"Log likelihood became {self.log_like} on "
                                 f"iteration {self.iter} of {self}")
            # Check for convergence.
            if self.delta_log_like < 0.0:
                # The log likelihood should not decrease.
                logger.warning("Log likelihood decreased from "
                               f"{self.log_likes[-2]} to {self.log_like} on "
                               f"iteration {self.iter}")
            elif (self.delta_log_like < self.conv_thresh
                  and self.iter >= self.min_iter):
                # Converge if the increase in log likelihood is
                # smaller than the convergence cutoff and at least
                # the minimum number of iterations have been run.
                self.converged = True
                logger.info(f"Log likelihood of {self} converged to "
                            f"{self.log_like} after {self.iter} iterations")
                break
        else:
            # The log likelihood did not converge within the maximum
            # number of iterations.
            logger.warning(f"Log likelihood of {self} failed to converge "
                           f"within {self.max_iter} iterations")
        # Return this instance so that any code that runs multiple
        # EmClustering objects can create and run them in one line.
        return self

    @staticmethod
    def props_index():
        return ["Real", "Observed"]

    def output_props(self):
        """ Return a DataFrame of the real and observed proportion of
        each cluster. """
        return pd.DataFrame(np.vstack([self.prop_real, self.prop_obs]),
                            index=self.props_index(),
                            columns=self.cluster_nums)

    def output_mus(self):
        """ Return a DataFrame of the mutation rate at each position
        for each cluster. """
        return pd.DataFrame(self.mus.T,
                            index=self.bvec.positions,
                            columns=self.cluster_nums)

    def output_resps(self):
        """ Return the responsibilities of the reads. """
        return pd.DataFrame(self.resps.T[self.bvec.bvec_to_read],
                            index=self.bvec.read_names,
                            columns=self.cluster_nums)

    @cached_property
    def description(self):
        """ Cache the string of this EM clustering object to speed up
        the very frequent logging messages that use this string. """
        return f"EM Clustering of {self.bvec}, replicate {self.rep}"

    def __str__(self):
        return self.description
