from functools import cached_property
from logging import getLogger

import numpy as np
import pandas as pd
from scipy.optimize import newton_krylov
from scipy.special import logsumexp
from scipy.stats import dirichlet

from .bvec import BitVector

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
    return n_params * np.log(n_data) - 2. * log_like


def calc_log_denom(log_mus: np.ndarray, log_nos: np.ndarray, min_mut_dist: int):
    """ Calculate the probability that a bit vector generated randomly
    from the given mutation rates would not have any mutations closer
    than ```min_mut_dist```. Note that ```mus``` is transposed relative
    to all other uses, so that its shape is (positions x clusters).
    Transposition makes the indexing easier because this function uses
    the positions as the primary axis and clusters as secondary.

    Parameters
    ----------
    log_mus: ndarray
        A (positions x clusters) array of the logs of the probabilities
        that a base at each position in each cluster is mutated.
    log_nos: ndarray
        A (positions x clusters) array of the logs of the probabilities
        that a base at each position in each cluster is not mutated.
    min_mut_dist: int
        Minimum permitted distance between consecutive mutations.
    """
    if log_mus.shape != log_nos.shape:
        raise ValueError(f"Shapes of log_mus {log_mus.shape} "
                         f"and log_nos {log_nos.shape} differ")
    # Determine the number of positions and clusters.
    npos, ncls = log_mus.shape
    # Probability that the bit vector has no two mutations that are too
    # close given that no two mutations before base (j) are too close.
    log_prob_j = np.empty_like(log_mus, dtype=float)
    # Probability that the bit vector has no two mutations that are too
    # close given that no two mutations before base (j) are too close
    # and that base (j) is mutated.
    log_prob_j_if_j_mut = np.empty_like(log_mus, dtype=float)
    # Set the last row of each array to 100% probability: log(P) = 0.
    log_prob_j[-1].fill(0.)
    log_prob_j_if_j_mut[-1].fill(0.)
    # Loop from the penultimate (npos - 2) to the first position (j).
    for j in range(npos - 2, -1, -1):
        # Find the probability that, given that base (j) is mutated,
        # none of the next (min_mut_dist - 1) bases are mutated,
        # i.e. that no other mutation is too close to base (j).
        log_prob_j_if_j_mut[j] = np.sum(log_nos[j + 1: j + min_mut_dist],
                                        axis=0)
        # Multiply it by the probability that no two mutations after the
        # next (min_mut_dist - 1) bases are too close, which is equal to
        # log_prob_j[j + min_mut_dist].
        if j + min_mut_dist < npos:
            log_prob_j_if_j_mut[j] += log_prob_j[j + min_mut_dist]
        # The probability that no two mutations are too close between
        # position (j) and the 3' end of the sequence is the probability
        # that base (j) is mutated (mus[j]) times the probability that
        # no two mutations afterward are too close given that base (j)
        # is mutated (log_prob_j_if_j_mut[j]), plus the probability that
        # base (j) is not mutated (nos[j]) times the probability that no
        # two mutations afterward are too close given that base (j) is
        # not mutated (log_prob_j[j + 1]).
        log_prob_j[j] = np.log(np.exp(log_mus[j]
                                      + log_prob_j_if_j_mut[j])
                               + np.exp(log_nos[j]
                                        + log_prob_j[j + 1]))
    # The denominators are in the first row of probabilities.
    log_denom = log_prob_j[0]
    return log_denom, log_prob_j_if_j_mut


def dmus(curr_mus: np.ndarray, init_mus: np.ndarray, min_mut_dist: int):
    """ Compute the derivative of the mus correction function.

    Parameters
    ----------
    curr_mus: ndarray
        A (positions x clusters) array of the current guesses of each
        cluster's (corrected) mutation rates.
    init_mus: ndarray
        A (positions x clusters) array of the initial mutation rates of
        each cluster, from the .
    min_mut_dist: int
        Minimum permitted distance between consecutive mutations.
    """
    # Compute the logs of the current mutation rates.
    log_curr_mus = np.log(curr_mus)
    log_curr_nos = np.log(1. - curr_mus)
    # Compute the denominator of each cluster in forward and reverse.
    log_denom_fwd, log_prob_mut_fwd = calc_log_denom(log_curr_mus,
                                                     log_curr_nos,
                                                     min_mut_dist)
    log_denom_rev, log_prob_mut_rev = calc_log_denom(log_curr_mus[::-1],
                                                     log_curr_nos[::-1],
                                                     min_mut_dist)
    if not np.allclose(log_denom_fwd, log_denom_rev):
        raise ValueError(f"Denominators in forward ({log_denom_fwd})"
                         f"and reverse ({log_denom_rev}) directions differ")
    # Compute the derivative.
    return np.exp(log_curr_mus + log_prob_mut_fwd
                  + log_prob_mut_rev[::-1] - log_denom_fwd) - init_mus


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
        if not conv_thresh >= 0.:
            raise ValueError(f"conv_thresh must be ≥ 0, but got {conv_thresh}")
        self.conv_thresh = conv_thresh
        # Minimum and maximum permitted mutation rates
        if not 0. < epsilon <= 0.5:
            raise ValueError(f"epsilon must be in (0, 0.5], but got {epsilon}")
        self.min_mu = epsilon
        self.max_mu = 1. - epsilon
        # Number of reads observed in each cluster (no bias correction)
        self.nreads = np.empty(self.ncls, dtype=float)
        # Denominator of each cluster
        self.log_denom = np.empty(self.ncls, dtype=float)
        # Mutation rate of each position (row) in each cluster (col)
        self.mus = np.empty((self.bvec.n_pos_use, self.ncls), dtype=float)
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
    def sparse_mus(self):
        """ Mutation rate of each position (row) in each cluster (col),
        including all positions that are not used for clustering. """
        sparse_mus = np.zeros((len(self.bvec.seq), self.ncls), dtype=float)
        # Copy mutation rates to the rows of sparse_mu that correspond
        # to used positions (unused rows remain zero).
        sparse_mus[self.sparse_pos] = self.mus
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
        # Calculate the number of reads in each cluster by summing the
        # count-weighted likelihood that each bit vector came from the
        # cluster.
        self.nreads = np.dot(self.resps, self.bvec.counts)
        # Loop over each position (j).
        for j, muts_j in enumerate(self.bvec.muts):
            # Calculate the number of mutations at each position in each
            # cluster by summing the count-weighted likelihood that each
            # bit vector with a mutation at (j) came from the cluster.
            self.mus[j] = np.dot(self.resps[:, muts_j],
                                 self.bvec.counts[muts_j]) / self.nreads
        logger.debug(f"Computed uncorrected mus of {self}:\n{self.mus}")
        # Solve for the corrected mutation rates for the cluster using
        # Newton's method with the Krylov approximation. Compute the
        # derivative of the error using dmus, which compares the current
        # solution for mus[k] (x) with the uncorrected mutation rates
        # (sparse_mus). Also use sparse_mus as the initial guess for the
        # solution, since it is usually close to the corrected solution.
        sparse_mus = self.sparse_mus
        sparse_mus = newton_krylov(lambda mus: dmus(mus,
                                                    sparse_mus,
                                                    self.bvec.min_mut_dist),
                                   sparse_mus)
        self.mus[:, :] = sparse_mus[self.sparse_pos]
        logger.debug(f"Computed corrected mus of {self}:\n{self.mus}")
        # Verify that mutation rates are within limits.
        self.mus[:, :] = np.clip(self.mus, self.min_mu, self.max_mu)

    def _exp_step(self):
        """ Run the Expectation step of the EM algorithm. """
        # SUB-STEP E1: Update the denominator of each cluster.
        # Compute the logs of the sparse mutation and no-mutation rates.
        sparse_mus = self.sparse_mus
        log_sparse_mus = np.log(sparse_mus)
        log_sparse_nos = np.log(1. - sparse_mus)
        # Compute the logs of the dense mutation and no-mutation rates.
        log_mus = log_sparse_mus[self.sparse_pos]
        log_nos = log_sparse_nos[self.sparse_pos]
        # Update the denominator of each cluster using the sparse mus.
        self.log_denom[:] = calc_log_denom(log_sparse_mus,
                                           log_sparse_nos,
                                           self.bvec.min_mut_dist)[0]
        logger.debug(f"Computed log denominators of {self}:\n{self.log_denom}")

        # SUB-STEP E2: Compute for each cluster and observed bit vector
        #              the joint probability that another random bit
        #              vector would both come from the same cluster and
        #              have exactly the same set of bits.
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
            for j in range(self.bvec.n_pos_use):
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
                adjust_for_mut_j_given_k = log_mus[j, k] - log_nos[j, k]
                self.resps[k][self.bvec.muts[j]] += adjust_for_mut_j_given_k
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
        conc_params = 1. - np.random.default_rng().random(self.ncls)
        # Initialize cluster membership with a Dirichlet distribution.
        self.resps[:, :] = dirichlet.rvs(conc_params, self.bvec.n_bvec_use).T
        # Run EM until the log likelihood converges or the number of
        # iterations reaches max_iter, whichever happens first.
        self.converged = False
        for self.iter in range(1, self.max_iter + 1):
            logger.debug(f"Began iteration {self.iter} of {self}")
            with np.errstate(divide="ignore"):
                # Update the mutation rates and cluster proportions.
                self._max_step()
                # Update the cluster membership and log likelihood.
                self._exp_step()
            logger.info(f"Ended iteration {self.iter} of {self} "
                        f"(log likelihood = {self.log_like})")
            if not np.isfinite(self.log_like):
                raise ValueError(f"Log likelihood of {self} became "
                                 f"{self.log_like} on iteration {self.iter}")
            # Check for convergence.
            if self.delta_log_like < 0.:
                # The log likelihood should not decrease.
                logger.warning(f"Log likelihood of {self} decreased from "
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
        return pd.DataFrame(self.mus,
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
        return (f"EM Clustering of {self.bvec} into {self.ncls} clusters, "
                f"replicate run {self.rep}")

    def __str__(self):
        return self.description
