"""
DREEM Testing Core Module
Unit tests for all functions in package `core`.
------------------------------------------------------------------------
Auth: Matty
Date: 2023-06-05
========================================================================
"""

from itertools import chain, combinations, product
from string import printable
from typing import Generator, Sequence
import unittest as ut

import numpy as np

from .mu import calc_mu_adj, calc_mu_obs
from .rel import (IRREC, MATCH, DELET,
                  INS_5, INS_3, INS_8, MINS5, MINS3, ANY_8,
                  SUB_A, SUB_C, SUB_G, SUB_T, SUB_N,
                  ANY_N, INDEL, NOCOV,
                  MIN_QUAL, MAX_QUAL,
                  validate_relvec, iter_relvecs_q53, iter_relvecs_all,
                  relvec_to_read, as_sam)
from .seq import BASES, BASES_ARR, RBASE, DNA, RNA, expand_degenerate_seq
from .sim import rand_dna


# General functions ####################################################

def has_close_muts(bitvec: np.ndarray, min_gap: int):
    """ Return True if the bit vector has two mutations separated by
    fewer than `min_gap` non-mutated bits, otherwise False. """
    if bitvec.ndim != 1:
        raise ValueError(f"bitvec must have 1 dimension, but got {bitvec.ndim}")
    if min_gap < 0:
        raise ValueError(f"min_gap must be ≥ 0, but got {min_gap}")
    if min_gap == 0:
        # Two mutations cannot be separated by fewer than 0 positions.
        return False
    # Close mutations are separated from each other by less than min_gap
    # non-mutated bits. Equivalently, the difference in their positions
    # (their distance) is < min_gap + 1, or ≤ min_gap. These distances
    # are computed as the differences (using np.diff) in the positions
    # of consecutive mutations (using np.flatnonzero).
    dists = np.diff(np.flatnonzero(bitvec))
    return dists.size > 0 and np.min(dists) <= min_gap


def label_close_muts(bitvecs: np.ndarray, min_gap: int):
    """ Return a 1D vector that is True for every row in `bitvecs` that
    has two mutations that are too close, and otherwise False. """
    if bitvecs.ndim != 2:
        raise ValueError(
            f"bitvects must have 2 dimensions, but got {bitvecs.ndim}")
    return np.array([has_close_muts(bitvec, min_gap) for bitvec in bitvecs])


def drop_close_muts(bitvecs: np.ndarray, min_gap: int):
    """ Return a new array without every row in `bitvecs` that has two
    mutations that are too close. """
    return bitvecs[np.logical_not(label_close_muts(bitvecs, min_gap))]


# General function tests ###############################################

class TestTestHasCloseMuts(ut.TestCase):
    """ Test function `test.has_close_muts`. """

    def test_no_muts(self):
        """ Test that a bit vector with no mutations returns False. """
        bitvec = np.zeros(5, dtype=bool)
        for g in range(bitvec.size):
            self.assertFalse(has_close_muts(bitvec, g))

    def test_one_mut(self):
        """ Test that a bit vector with one mutation returns False. """
        for i in range(5):
            bitvec = np.zeros(5, dtype=bool)
            bitvec[i] = 1
            for g in range(bitvec.size):
                self.assertFalse(has_close_muts(bitvec, g))

    def test_more_muts(self):
        """ Test every bit vector with 2 - 5 mutations. """
        bitvecs_gaps = {
            # 2 mutations (n = 10)
            (0, 0, 0, 1, 1): 0,
            (0, 0, 1, 0, 1): 1,
            (0, 0, 1, 1, 0): 0,
            (0, 1, 0, 0, 1): 2,
            (0, 1, 0, 1, 0): 1,
            (0, 1, 1, 0, 0): 0,
            (1, 0, 0, 0, 1): 3,
            (1, 0, 0, 1, 0): 2,
            (1, 0, 1, 0, 0): 1,
            (1, 1, 0, 0, 0): 0,
            # 3 mutations (n = 10)
            (0, 0, 1, 1, 1): 0,
            (0, 1, 0, 1, 1): 0,
            (0, 1, 1, 0, 1): 0,
            (0, 1, 1, 1, 0): 0,
            (1, 0, 0, 1, 1): 0,
            (1, 0, 1, 0, 1): 1,
            (1, 0, 1, 1, 0): 0,
            (1, 1, 0, 0, 1): 0,
            (1, 1, 0, 1, 0): 0,
            (1, 1, 1, 0, 0): 0,
            # 4 mutations (n = 5)
            (0, 1, 1, 1, 1): 0,
            (1, 0, 1, 1, 1): 0,
            (1, 1, 0, 1, 1): 0,
            (1, 1, 1, 0, 1): 0,
            (1, 1, 1, 1, 0): 0,
            # 5 mutations (n = 1)
            (1, 1, 1, 1, 1): 0,
        }
        for bitvec, bit_gap in bitvecs_gaps.items():
            for min_gap in range(len(bitvec)):
                self.assertEqual(has_close_muts(np.array(bitvec), min_gap),
                                 bit_gap < min_gap)


# Module: mu ###########################################################

class TestMuCalcMuAdj(ut.TestCase):
    """ Test function `mu.calc_mu_adj`. """

    def test_mu_multiplex(self):
        """ Test that running 1 - 5 clusters simultaneously produces the
        same results as running each cluster separately. """
        n_pos = 16
        min_k, max_k = 1, 5
        min_g, max_g = 0, 4
        max_m = 0.1
        # Test each number of clusters (k).
        for k in range(min_k, max_k + 1):
            # Test each minimum gap between mutations (g).
            for g in range(min_g, max_g + 1):
                with self.subTest(k=k, g=g):
                    # Generate random observed mutation rates.
                    mus_obs = np.random.default_rng().random((n_pos, k)) * max_m
                    # Adjust all rates simultaneously.
                    mus_adj_sim = calc_mu_adj(mus_obs, g)
                    # Adjust the rates of each cluster (i) separately.
                    mus_adj_sep = np.empty_like(mus_obs)
                    for i in range(k):
                        mus_obs_i = mus_obs[:, i].reshape((n_pos, 1))
                        mus_adj_i = calc_mu_adj(mus_obs_i, g).reshape(n_pos)
                        mus_adj_sep[:, i] = mus_adj_i
                    # Compare the results.
                    self.assertTrue(np.allclose(mus_adj_sim, mus_adj_sep))

    def test_inv_calc_mu_obs(self):
        """ Test that this function inverts `mu.calc_mu_obs`. """
        n_pos = 16
        min_k, max_k = 1, 5
        min_g, max_g = 0, 4
        max_m = 0.2
        # Test each number of clusters (k).
        for k in range(min_k, max_k + 1):
            # Generate random real mutation rates.
            mus = np.random.default_rng().random((n_pos, k)) * max_m
            # Test each minimum gap between mutations (g).
            for g in range(min_g, max_g + 1):
                with self.subTest(k=k, g=g):
                    # Compute the observed mutation rates.
                    mus_obs = calc_mu_obs(mus, g)
                    # Adjust the observed mutation rates.
                    mus_adj = calc_mu_adj(mus_obs, g)
                    # Test if adjusted and initial mutation rates match.
                    self.assertTrue(np.allclose(mus_adj, mus))


class TestMuCalcMuObs(ut.TestCase):
    """ Test function `mu.calc_mu_obs`. """

    @ut.skip("Takes a long time to run: burdensome while debugging other tests")
    def test_obs_empirical(self):
        """ Test that this function accurately predicts the mutation
        rates that are actually observed when simulated bit vectors are
        filtered to remove mutations that are too close. """
        n_pos = 10
        min_g, max_g = 0, 4
        min_m, max_m = 0.01, 0.1
        # Choose the number of vectors to simulate as follows:
        # The number of mutations at each position in the simulated bit
        # vectors follows a binomial distribution, whose std. dev. is
        # sqrt(p * (1 - p) * n).
        # The proportion of mutations is this quantity divided by n:
        # sqrt(p * (1 - p) / n).
        # Choosing a tolerance of 3 std. dev. around the mean yields
        # 3 * sqrt(p * (1 - p) / n) ≤ tol
        # Solving for n (the number of vectors) gives
        # n ≥ p * (1 - p) / (tol / 3)^2 = p * (1 - p) * (2 / tol)^2
        nstdev = 3.  # number of standard deviations on each side
        tol = 5.e-4  # absolute tolerance for np.allclose
        n_vec = round(max_m * (1. - max_m) * (nstdev / tol) ** 2)
        # Generate random real mutation rates.
        mus = min_m + np.random.default_rng().random(n_pos) * (max_m - min_m)
        # Generate random bit vectors with the expected mutation rates.
        bvecs = None
        while bvecs is None or not np.allclose(np.mean(bvecs, axis=0), mus,
                                               atol=tol, rtol=0.):
            bvecs = np.less(np.random.default_rng().random((n_vec, n_pos)), mus)
        # Test each minimum gap between mutations (g).
        for g in range(min_g, max_g + 1):
            with self.subTest(g=g):
                # Drop bit vectors with mutations too close.
                bvecs_g = drop_close_muts(bvecs, g)
                # Compute the empirically observed mutation rates.
                mus_obs_emp = np.mean(bvecs_g, axis=0)
                # Predict the observed mutation rates with calc_mu_obs.
                mus_obs_prd = calc_mu_obs(mus.reshape((-1, 1)), g).reshape(-1)
                # Compare the empirical and predicted mutation rates.
                self.assertTrue(np.allclose(mus_obs_emp, mus_obs_prd,
                                            atol=tol, rtol=0.))

    def test_mu_multiplex(self):
        """ Test that running 1 - 5 clusters simultaneously produces the
        same results as running each cluster separately. """
        n_pos = 16
        min_k, max_k = 1, 5
        min_g, max_g = 0, 4
        max_m = 0.2
        # Test each number of clusters (k).
        for k in range(min_k, max_k + 1):
            # Test each minimum gap between mutations (g).
            for g in range(min_g, max_g + 1):
                with self.subTest(k=k, g=g):
                    # Generate random real mutation rates.
                    mus = np.random.default_rng().random((n_pos, k)) * max_m
                    # Adjust all rates simultaneously.
                    mus_obs_sim = calc_mu_obs(mus, g)
                    # Adjust the rates of each cluster (i) separately.
                    mus_obs_sep = np.empty_like(mus)
                    for i in range(k):
                        mus_i = mus[:, i].reshape((n_pos, 1))
                        mus_obs_i = calc_mu_obs(mus_i, g).reshape(n_pos)
                        mus_obs_sep[:, i] = mus_obs_i
                    # Compare the results.
                    self.assertTrue(np.allclose(mus_obs_sim, mus_obs_sep))

    def test_inv_calc_mu_adj(self):
        """ Test that this function inverts `mu.calc_mu_adj`. """
        n_pos = 16
        min_k, max_k = 1, 5
        min_g, max_g = 0, 4
        max_m = 0.1
        # Test each number of clusters (k).
        for k in range(min_k, max_k + 1):
            # Generate random observed mutation rates.
            mus_obs = np.random.default_rng().random((n_pos, k)) * max_m
            # Test each minimum gap between mutations (g).
            for g in range(min_g, max_g + 1):
                with self.subTest(k=k, g=g):
                    # Compute the adjusted mutation rates.
                    mus_adj = calc_mu_adj(mus_obs, g)
                    # Recompute the observed mutation rates.
                    mus_reobs = calc_mu_obs(mus_adj, g)
                    # Compare observed and reobserved mutation rates.
                    self.assertTrue(np.allclose(mus_obs, mus_reobs))


# Module: rel ##########################################################

class TestRelConstants(ut.TestCase):
    """ Test constants of `rel` module. """

    def test_primary_codes(self):
        """ Test the primary relation codes. """
        for exp, code in enumerate([MATCH, DELET, INS_5, INS_3,
                                    SUB_A, SUB_C, SUB_G, SUB_T]):
            self.assertEqual(code, 2 ** exp)

    def test_derived_codes(self):
        """ Test the derived relation codes. """
        self.assertEqual(IRREC, 0)
        self.assertEqual(MINS5, 5)
        self.assertEqual(MINS3, 9)
        self.assertEqual(INS_8, 12)
        self.assertEqual(ANY_8, 13)
        self.assertEqual(INDEL, 14)
        self.assertEqual(SUB_N, 240)
        self.assertEqual(ANY_N, 241)
        self.assertEqual(NOCOV, 255)


class TestRelValidateRelVec(ut.TestCase):
    """ Test function `rel.validate_relvecs`. """

    def test_matches(self):
        """ Test that a relation vector of all matches is valid. """
        self.assertIsNotNone(validate_relvec(np.ones(8, np.uint8)))

    def test_nocov_margin(self):
        """ Test that a relation vector with non-covered positions on
        its margins is valid. """
        relvec = np.full(8, NOCOV, np.uint8)
        relvec[4] = 1
        self.assertIsNotNone(validate_relvec(relvec))

    def test_nocov_middle(self):
        """ Test that a relation vector with non-covered positions in
        the middle is invalid. """
        relvec = np.full(8, MATCH, np.uint8)
        relvec[4] = NOCOV
        self.assertRaises(ValueError, validate_relvec, relvec)

    def test_blank(self):
        """ Test that an entirely blank relation vector is invalid. """
        self.assertRaises(ValueError, validate_relvec,
                          np.full(8, NOCOV, np.uint8))

    def test_not_array(self):
        """ Test that a non-array relation vector is invalid. """
        self.assertRaises(TypeError, validate_relvec,
                          np.ones(8, np.uint8).tolist())

    def test_not_uint8(self):
        """ Test that a non-uint8 relation vector is invalid. """
        self.assertRaises(TypeError, validate_relvec, np.ones(8, np.int8))


class TestRelIterRelvecsQ53(ut.TestCase):
    """ Test function `rel.iter_relvecs_q53`. """

    @staticmethod
    def iter_rels(seq: str, low_qual: Sequence[int] = (),
                  end5: int | None = None, end3: int | None = None):
        """ Convenience function to run `rel.iter_relvecs_q53` from a
        sequence of str and return a list of lists of ints. """
        return list(map(np.ndarray.tolist,
                        iter_relvecs_q53(DNA(seq.encode()),
                                         low_qual, end5, end3)))

    def test_type(self):
        """ Test that the result is a Generator of NumPy arrays. """
        self.assertTrue(isinstance(iter_relvecs_q53(DNA(b"A")), Generator))
        self.assertTrue(all(isinstance(relvec, np.ndarray)
                            for relvec in iter_relvecs_q53(DNA(b"A"))))
        self.assertIs(list(iter_relvecs_q53(DNA(b"A")))[0].dtype.type, np.uint8)

    def test_a(self):
        """ Test with sequence 'A'. """
        self.assertEqual(self.iter_rels("A"),
                         [[MATCH],
                          [MATCH + INS_3],
                          [MATCH + INS_5 + INS_3],
                          [MATCH + INS_5],
                          [SUB_C],
                          [SUB_C + INS_3],
                          [SUB_C + INS_5 + INS_3],
                          [SUB_C + INS_5],
                          [SUB_G],
                          [SUB_G + INS_3],
                          [SUB_G + INS_5 + INS_3],
                          [SUB_G + INS_5],
                          [SUB_T],
                          [SUB_T + INS_3],
                          [SUB_T + INS_5 + INS_3],
                          [SUB_T + INS_5]])

    def test_c(self):
        """ Test with sequence 'C'. """
        self.assertEqual(self.iter_rels("C"),
                         [[SUB_A],
                          [SUB_A + INS_3],
                          [SUB_A + INS_5 + INS_3],
                          [SUB_A + INS_5],
                          [MATCH],
                          [MATCH + INS_3],
                          [MATCH + INS_5 + INS_3],
                          [MATCH + INS_5],
                          [SUB_G],
                          [SUB_G + INS_3],
                          [SUB_G + INS_5 + INS_3],
                          [SUB_G + INS_5],
                          [SUB_T],
                          [SUB_T + INS_3],
                          [SUB_T + INS_5 + INS_3],
                          [SUB_T + INS_5]])

    def test_g(self):
        """ Test with sequence 'G'. """
        self.assertEqual(self.iter_rels("G"),
                         [[SUB_A],
                          [SUB_A + INS_3],
                          [SUB_A + INS_5 + INS_3],
                          [SUB_A + INS_5],
                          [SUB_C],
                          [SUB_C + INS_3],
                          [SUB_C + INS_5 + INS_3],
                          [SUB_C + INS_5],
                          [MATCH],
                          [MATCH + INS_3],
                          [MATCH + INS_5 + INS_3],
                          [MATCH + INS_5],
                          [SUB_T],
                          [SUB_T + INS_3],
                          [SUB_T + INS_5 + INS_3],
                          [SUB_T + INS_5]])

    def test_t(self):
        """ Test with sequence 'T'. """
        self.assertEqual(self.iter_rels("T"),
                         [[SUB_A],
                          [SUB_A + INS_3],
                          [SUB_A + INS_5 + INS_3],
                          [SUB_A + INS_5],
                          [SUB_C],
                          [SUB_C + INS_3],
                          [SUB_C + INS_5 + INS_3],
                          [SUB_C + INS_5],
                          [SUB_G],
                          [SUB_G + INS_3],
                          [SUB_G + INS_5 + INS_3],
                          [SUB_G + INS_5],
                          [MATCH],
                          [MATCH + INS_3],
                          [MATCH + INS_5 + INS_3],
                          [MATCH + INS_5]])

    def test_low_qual(self):
        """ Test with each low-quality base. """
        for base, low_qual in zip("ACGT", [225, 209, 177, 113]):
            self.assertEqual(self.iter_rels(base, [1]),
                             [[low_qual],
                              [low_qual + INS_3],
                              [low_qual + INS_5 + INS_3],
                              [low_qual + INS_5]])

    def test_low_qual_invalid(self):
        """ Test that invalid low-qual positions raise ValueError. """
        seq = "ACGT"
        for n in range(1, len(seq) + 1):
            self.assertTrue(isinstance(self.iter_rels(seq[: n], [1]), list))
            self.assertTrue(isinstance(self.iter_rels(seq[: n], [n]), list))
            self.assertRaises(ValueError, self.iter_rels, seq[: n], [0])
            self.assertRaises(ValueError, self.iter_rels, seq[: n], [n + 1])

    def test_xxax(self):
        """ Test that bases with no coverage are marked. """
        self.assertEqual(self.iter_rels("TCAG", end5=3, end3=3),
                         [[NOCOV, NOCOV, MATCH, NOCOV],
                          [NOCOV, NOCOV, MATCH + INS_3, NOCOV],
                          [NOCOV, NOCOV, MATCH + INS_5 + INS_3, NOCOV],
                          [NOCOV, NOCOV, MATCH + INS_5, NOCOV],
                          [NOCOV, NOCOV, SUB_C, NOCOV],
                          [NOCOV, NOCOV, SUB_C + INS_3, NOCOV],
                          [NOCOV, NOCOV, SUB_C + INS_5 + INS_3, NOCOV],
                          [NOCOV, NOCOV, SUB_C + INS_5, NOCOV],
                          [NOCOV, NOCOV, SUB_G, NOCOV],
                          [NOCOV, NOCOV, SUB_G + INS_3, NOCOV],
                          [NOCOV, NOCOV, SUB_G + INS_5 + INS_3, NOCOV],
                          [NOCOV, NOCOV, SUB_G + INS_5, NOCOV],
                          [NOCOV, NOCOV, SUB_T, NOCOV],
                          [NOCOV, NOCOV, SUB_T + INS_3, NOCOV],
                          [NOCOV, NOCOV, SUB_T + INS_5 + INS_3, NOCOV],
                          [NOCOV, NOCOV, SUB_T + INS_5, NOCOV]])

    def test_agg(self):
        """ Test with sequence 'AGG'. """
        rels = self.iter_rels("AGG")
        self.assertEqual(rels[: 47],
                         [[MATCH, SUB_A, SUB_A],
                          [MATCH + INS_3, SUB_A, SUB_A],
                          [MATCH + INS_5 + INS_3, SUB_A + INS_3, SUB_A],
                          [MATCH + INS_3, SUB_A + INS_5, SUB_A + INS_3],
                          [MATCH + INS_3, SUB_A, SUB_A + INS_5],
                          [MATCH + INS_5, SUB_A + INS_3, SUB_A],
                          [MATCH + INS_5, SUB_A + INS_5 + INS_3, SUB_A + INS_3],
                          [MATCH + INS_5, SUB_A + INS_3, SUB_A + INS_5],
                          [MATCH, SUB_A + INS_5, SUB_A + INS_3],
                          [MATCH, SUB_A + INS_5, SUB_A + INS_5 + INS_3],
                          [MATCH, SUB_A, SUB_A + INS_5],
                          [MATCH, SUB_A, SUB_C],
                          [MATCH + INS_3, SUB_A, SUB_C],
                          [MATCH + INS_5 + INS_3, SUB_A + INS_3, SUB_C],
                          [MATCH + INS_3, SUB_A + INS_5, SUB_C + INS_3],
                          [MATCH + INS_3, SUB_A, SUB_C + INS_5],
                          [MATCH + INS_5, SUB_A + INS_3, SUB_C],
                          [MATCH + INS_5, SUB_A + INS_5 + INS_3, SUB_C + INS_3],
                          [MATCH + INS_5, SUB_A + INS_3, SUB_C + INS_5],
                          [MATCH, SUB_A + INS_5, SUB_C + INS_3],
                          [MATCH, SUB_A + INS_5, SUB_C + INS_5 + INS_3],
                          [MATCH, SUB_A, SUB_C + INS_5],
                          [MATCH, SUB_A, MATCH],
                          [MATCH + INS_3, SUB_A, MATCH],
                          [MATCH + INS_5 + INS_3, SUB_A + INS_3, MATCH],
                          [MATCH + INS_3, SUB_A + INS_5, MATCH + INS_3],
                          [MATCH + INS_3, SUB_A, MATCH + INS_5],
                          [MATCH + INS_5, SUB_A + INS_3, MATCH],
                          [MATCH + INS_5, SUB_A + INS_5 + INS_3, MATCH + INS_3],
                          [MATCH + INS_5, SUB_A + INS_3, MATCH + INS_5],
                          [MATCH, SUB_A + INS_5, MATCH + INS_3],
                          [MATCH, SUB_A + INS_5, MATCH + INS_5 + INS_3],
                          [MATCH, SUB_A, MATCH + INS_5],
                          [MATCH, SUB_A, SUB_T],
                          [MATCH + INS_3, SUB_A, SUB_T],
                          [MATCH + INS_5 + INS_3, SUB_A + INS_3, SUB_T],
                          [MATCH + INS_3, SUB_A + INS_5, SUB_T + INS_3],
                          [MATCH + INS_3, SUB_A, SUB_T + INS_5],
                          [MATCH + INS_5, SUB_A + INS_3, SUB_T],
                          [MATCH + INS_5, SUB_A + INS_5 + INS_3, SUB_T + INS_3],
                          [MATCH + INS_5, SUB_A + INS_3, SUB_T + INS_5],
                          [MATCH, SUB_A + INS_5, SUB_T + INS_3],
                          [MATCH, SUB_A + INS_5, SUB_T + INS_5 + INS_3],
                          [MATCH, SUB_A, SUB_T + INS_5],
                          [MATCH, SUB_A, DELET],
                          [MATCH + INS_3, SUB_A, DELET],
                          [MATCH, SUB_C, SUB_A]])


class TestRelIterRelvecsAll(ut.TestCase):
    """ Test function `rel.iter_relvecs_all`. """

    def assert_equal(self, ref: DNA, expects: list):
        """ Check that the expected and actual results match. """
        for exp, res in zip(chain(*expects),
                            iter_relvecs_all(ref),
                            strict=True):
            with self.subTest(exp=exp, res=res):
                self.assertTrue(np.all(exp == res))

    def test_length_1(self):
        """ Test with all length-1 DNA sequences. """
        for ref in expand_degenerate_seq(b"N"):
            expects = [
                iter_relvecs_q53(ref, [], 1, 1),
                iter_relvecs_q53(ref, [1], 1, 1),
            ]
            self.assert_equal(ref, expects)

    def test_length_2(self):
        """ Test with all length-2 DNA sequences. """
        for ref in expand_degenerate_seq(b"NN"):
            expects = [
                iter_relvecs_q53(ref, [], 1, 1),
                iter_relvecs_q53(ref, [1], 1, 1),
                iter_relvecs_q53(ref, [], 1, 2),
                iter_relvecs_q53(ref, [1], 1, 2),
                iter_relvecs_q53(ref, [2], 1, 2),
                iter_relvecs_q53(ref, [1, 2], 1, 2),
                iter_relvecs_q53(ref, [], 2, 2),
                iter_relvecs_q53(ref, [2], 2, 2),
            ]
            self.assert_equal(ref, expects)

    def test_length_3(self):
        """ Test with all length-3 DNA sequences. """
        for ref in expand_degenerate_seq(b"NNN"):
            expects = [
                iter_relvecs_q53(ref, [], 1, 1),
                iter_relvecs_q53(ref, [1], 1, 1),
                iter_relvecs_q53(ref, [], 1, 2),
                iter_relvecs_q53(ref, [1], 1, 2),
                iter_relvecs_q53(ref, [2], 1, 2),
                iter_relvecs_q53(ref, [1, 2], 1, 2),
                iter_relvecs_q53(ref, [], 1, 3),
                iter_relvecs_q53(ref, [1], 1, 3),
                iter_relvecs_q53(ref, [2], 1, 3),
                iter_relvecs_q53(ref, [3], 1, 3),
                iter_relvecs_q53(ref, [1, 2], 1, 3),
                iter_relvecs_q53(ref, [1, 3], 1, 3),
                iter_relvecs_q53(ref, [2, 3], 1, 3),
                iter_relvecs_q53(ref, [1, 2, 3], 1, 3),
                iter_relvecs_q53(ref, [], 2, 2),
                iter_relvecs_q53(ref, [2], 2, 2),
                iter_relvecs_q53(ref, [], 2, 3),
                iter_relvecs_q53(ref, [2], 2, 3),
                iter_relvecs_q53(ref, [3], 2, 3),
                iter_relvecs_q53(ref, [2, 3], 2, 3),
                iter_relvecs_q53(ref, [], 3, 3),
                iter_relvecs_q53(ref, [3], 3, 3),
            ]
            self.assert_equal(ref, expects)


class TestRelRelvecToRead(ut.TestCase):
    """ Test function `rel.relvec_to_read`. """

    def assert_equal(self, ref: DNA,
                     relvecs: list[list[int]],
                     expects: list[tuple[bytes, bytes, bytes, int, int]]):
        """ Assert that the actual and expected outputs match. """
        for relvec, expect in zip(relvecs, expects, strict=True):
            with self.subTest(relvec=relvec, expect=expect):
                self.assertEqual(relvec_to_read(ref, np.array(relvec,
                                                              dtype=np.uint8),
                                                MAX_QUAL, MIN_QUAL),
                                 expect)

    def assert_raise(self, ref: DNA,
                     relvecs: list[list[int]],
                     error: type[Exception],
                     regex: str):
        """ Assert that the relation vectors raise an exception. """
        for relvec in relvecs:
            with self.subTest(relvec=relvec):
                self.assertRaisesRegex(error, regex, relvec_to_read,
                                       ref, np.array(relvec,
                                                     dtype=np.uint8),
                                       MAX_QUAL, MIN_QUAL)

    def test_all_match(self):
        """ Test when the read has four matching bases. """
        ref = DNA(b"ACGT")
        relvecs = [[MATCH, MATCH, MATCH, MATCH]]
        expects = [(b"ACGT", b"IIII", b"4=", 1, 4)]
        self.assert_equal(ref, relvecs, expects)

    def test_nocov_valid(self):
        """ Test when the read does not cover one or both ends of the
        reference. """
        ref = DNA(b"ACGT")
        relvecs = [
            [NOCOV, MATCH, MATCH, MATCH],
            [MATCH, MATCH, MATCH, NOCOV],
            [NOCOV, MATCH, MATCH, NOCOV],
            [NOCOV, NOCOV, MATCH, MATCH],
            [MATCH, MATCH, NOCOV, NOCOV],
        ]
        expects = [
            (b"CGT", b"III", b"3=", 2, 4),
            (b"ACG", b"III", b"3=", 1, 3),
            (b"CG", b"II", b"2=", 2, 3),
            (b"GT", b"II", b"2=", 3, 4),
            (b"AC", b"II", b"2=", 1, 2),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_nocov_middle_invalid(self):
        """ Test when the read does not cover a middle position. """
        ref = DNA(b"ACGT")
        relvecs = [
            [MATCH, NOCOV, MATCH, MATCH],
            [MATCH, MATCH, NOCOV, MATCH],
            [MATCH, NOCOV, NOCOV, MATCH],
        ]
        self.assert_raise(ref, relvecs, ValueError,
                          "Expected [0-9]+ base calls")

    def test_nocov_all_invalid(self):
        """ Test when the read does not cover any positions. """
        ref = DNA(b"ACGT")
        relvecs = [[NOCOV, NOCOV, NOCOV, NOCOV]]
        self.assert_raise(ref, relvecs, ValueError,
                          "Relation vector is blank")

    def test_low_qual_valid(self):
        """ Test when the read has a low-quality base. """
        ref = DNA(b"ACGT")
        relvecs = [
            [ANY_N - SUB_A, MATCH, MATCH, MATCH],
            [MATCH, ANY_N - SUB_C, MATCH, MATCH],
            [MATCH, MATCH, ANY_N - SUB_G, MATCH],
            [MATCH, MATCH, MATCH, ANY_N - SUB_T],
        ]
        expects = [
            (b"NCGT", b"!III", b"1M3=", 1, 4),
            (b"ANGT", b"I!II", b"1=1M2=", 1, 4),
            (b"ACNT", b"II!I", b"2=1M1=", 1, 4),
            (b"ACGN", b"III!", b"3=1M", 1, 4),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_low_qual_invalid(self):
        """ Test when the read has an invalid low-quality base. """
        ref = DNA(b"ACGT")
        relvecs = [
            [ANY_N, MATCH, MATCH, MATCH],
            [MATCH, ANY_N, MATCH, MATCH],
            [MATCH, MATCH, ANY_N, MATCH],
            [MATCH, MATCH, MATCH, ANY_N],
        ]
        self.assert_raise(ref, relvecs, ValueError,
                          f"Invalid relation {ANY_N}")

    def test_subst_valid(self):
        """ Test when the read has a substitution. """
        ref = DNA(b"ACGT")
        relvecs = [
            [SUB_C, MATCH, MATCH, MATCH],
            [SUB_G, MATCH, MATCH, MATCH],
            [SUB_T, MATCH, MATCH, MATCH],
            [MATCH, SUB_A, MATCH, MATCH],
            [MATCH, SUB_G, MATCH, MATCH],
            [MATCH, SUB_T, MATCH, MATCH],
            [MATCH, MATCH, SUB_A, MATCH],
            [MATCH, MATCH, SUB_C, MATCH],
            [MATCH, MATCH, SUB_T, MATCH],
            [MATCH, MATCH, MATCH, SUB_A],
            [MATCH, MATCH, MATCH, SUB_C],
            [MATCH, MATCH, MATCH, SUB_G],
        ]
        expects = [
            (b"CCGT", b"IIII", b"1X3=", 1, 4),
            (b"GCGT", b"IIII", b"1X3=", 1, 4),
            (b"TCGT", b"IIII", b"1X3=", 1, 4),
            (b"AAGT", b"IIII", b"1=1X2=", 1, 4),
            (b"AGGT", b"IIII", b"1=1X2=", 1, 4),
            (b"ATGT", b"IIII", b"1=1X2=", 1, 4),
            (b"ACAT", b"IIII", b"2=1X1=", 1, 4),
            (b"ACCT", b"IIII", b"2=1X1=", 1, 4),
            (b"ACTT", b"IIII", b"2=1X1=", 1, 4),
            (b"ACGA", b"IIII", b"3=1X", 1, 4),
            (b"ACGC", b"IIII", b"3=1X", 1, 4),
            (b"ACGG", b"IIII", b"3=1X", 1, 4),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_subst_invalid(self):
        """ Test when the read has an invalid substitution. """
        ref = DNA(b"ACGT")
        relvecs = [
            [SUB_A, MATCH, MATCH, MATCH],
            [MATCH, SUB_C, MATCH, MATCH],
            [MATCH, MATCH, SUB_G, MATCH],
            [MATCH, MATCH, MATCH, SUB_T],
        ]
        self.assert_raise(ref, relvecs, ValueError,
                          "Cannot substitute [ACGT] to itself")

    def test_delete_valid(self):
        """ Test when the read has deletions. """
        ref = DNA(b"ACGT")
        relvecs = [
            # 1 deletion (n=4)
            [DELET, MATCH, MATCH, MATCH],
            [MATCH, DELET, MATCH, MATCH],
            [MATCH, MATCH, DELET, MATCH],
            [MATCH, MATCH, MATCH, DELET],
            # 2 deletions (n=6)
            [DELET, DELET, MATCH, MATCH],
            [DELET, MATCH, DELET, MATCH],
            [DELET, MATCH, MATCH, DELET],
            [MATCH, DELET, DELET, MATCH],
            [MATCH, DELET, MATCH, DELET],
            [MATCH, MATCH, DELET, DELET],
        ]
        expects = [
            # 1 deletion (n=4)
            (b"CGT", b"III", b"1D3=", 1, 4),
            (b"AGT", b"III", b"1=1D2=", 1, 4),
            (b"ACT", b"III", b"2=1D1=", 1, 4),
            (b"ACG", b"III", b"3=1D", 1, 4),
            # 2 deletions (n=6)
            (b"GT", b"II", b"2D2=", 1, 4),
            (b"CT", b"II", b"1D1=1D1=", 1, 4),
            (b"CG", b"II", b"1D2=1D", 1, 4),
            (b"AT", b"II", b"1=2D1=", 1, 4),
            (b"AG", b"II", b"1=1D1=1D", 1, 4),
            (b"AC", b"II", b"2=2D", 1, 4),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_all_delete_invalid(self):
        """ Test when the read is entirely deleted. """
        ref = DNA(b"ACGT")
        relvecs = [[DELET, DELET, DELET, DELET]]
        self.assert_raise(ref, relvecs, ValueError,
                          "Read contained no bases")

    def test_insert_valid(self):
        """ Test when the read has insertions. """
        ref = DNA(b"ACGT")
        relvecs = [
            # 1 insertion (n=5)
            [MINS3, MATCH, MATCH, MATCH],
            [MINS5, MINS3, MATCH, MATCH],
            [MATCH, MINS5, MINS3, MATCH],
            [MATCH, MATCH, MINS5, MINS3],
            [MATCH, MATCH, MATCH, MINS5],
            # 2 insertions, 1 base apart (n=4)
            [ANY_8, MINS3, MATCH, MATCH],
            [MINS5, ANY_8, MINS3, MATCH],
            [MATCH, MINS5, ANY_8, MINS3],
            [MATCH, MATCH, MINS5, ANY_8],
            # 2 insertions, 2 bases apart (n=3)
            [MINS3, MINS5, MINS3, MATCH],
            [MINS5, MINS3, MINS5, MINS3],
            [MATCH, MINS5, MINS3, MINS5],
            # 2 insertions, 3 bases apart (n=2)
            [MINS3, MATCH, MINS5, MINS3],
            [MINS5, MINS3, MATCH, MINS5],
            # 2 insertions, 4 bases apart (n=1)
            [MINS3, MATCH, MATCH, MINS5],
            # 3 insertions, 1 base apart (n=3)
            [ANY_8, ANY_8, MINS3, MATCH],
            [MINS5, ANY_8, ANY_8, MINS3],
            [MATCH, MINS5, ANY_8, ANY_8],
            # 4 insertions, 1 base apart (n=2)
            [ANY_8, ANY_8, ANY_8, MINS3],
            [MINS5, ANY_8, ANY_8, ANY_8],
            # 5 insertions, 1 base apart (n=1)
            [ANY_8, ANY_8, ANY_8, ANY_8],
            # 1 insertion next to no coverage (n=3)
            [NOCOV, MINS3, MATCH, NOCOV],
            [NOCOV, MINS5, MINS3, NOCOV],
            [NOCOV, MATCH, MINS5, NOCOV],
            # 2 insertions, 1 base apart, next to no coverage (n=2)
            [NOCOV, ANY_8, MINS3, NOCOV],
            [NOCOV, MINS5, ANY_8, NOCOV],
            # 3 insertions, 1 base apart, next to no coverage (n=1)
            [NOCOV, ANY_8, ANY_8, NOCOV],
        ]
        expects = [
            # 1 insertion (n=5)
            (b"NACGT", b"IIIII", b"1I4=", 1, 4),
            (b"ANCGT", b"IIIII", b"1=1I3=", 1, 4),
            (b"ACNGT", b"IIIII", b"2=1I2=", 1, 4),
            (b"ACGNT", b"IIIII", b"3=1I1=", 1, 4),
            (b"ACGTN", b"IIIII", b"4=1I", 1, 4),
            # 2 insertions, 1 base apart (n=4)
            (b"NANCGT", b"IIIIII", b"1I1=1I3=", 1, 4),
            (b"ANCNGT", b"IIIIII", b"1=1I1=1I2=", 1, 4),
            (b"ACNGNT", b"IIIIII", b"2=1I1=1I1=", 1, 4),
            (b"ACGNTN", b"IIIIII", b"3=1I1=1I", 1, 4),
            # 2 insertions, 2 bases apart (n=3)
            (b"NACNGT", b"IIIIII", b"1I2=1I2=", 1, 4),
            (b"ANCGNT", b"IIIIII", b"1=1I2=1I1=", 1, 4),
            (b"ACNGTN", b"IIIIII", b"2=1I2=1I", 1, 4),
            # 2 insertions, 3 bases apart (n=2)
            (b"NACGNT", b"IIIIII", b"1I3=1I1=", 1, 4),
            (b"ANCGTN", b"IIIIII", b"1=1I3=1I", 1, 4),
            # 2 insertions, 4 bases apart (n=1)
            (b"NACGTN", b"IIIIII", b"1I4=1I", 1, 4),
            # 3 insertions, 1 base apart (n=3)
            (b"NANCNGT", b"IIIIIII", b"1I1=1I1=1I2=", 1, 4),
            (b"ANCNGNT", b"IIIIIII", b"1=1I1=1I1=1I1=", 1, 4),
            (b"ACNGNTN", b"IIIIIII", b"2=1I1=1I1=1I", 1, 4),
            # 4 insertions, 1 base apart (n=2)
            (b"NANCNGNT", b"IIIIIIII", b"1I1=1I1=1I1=1I1=", 1, 4),
            (b"ANCNGNTN", b"IIIIIIII", b"1=1I1=1I1=1I1=1I", 1, 4),
            # 5 insertions, 1 base apart (n=1)
            (b"NANCNGNTN", b"IIIIIIIII", b"1I1=1I1=1I1=1I1=1I", 1, 4),
            # 1 insertion next to no coverage (n=3)
            (b"NCG", b"III", b"1I2=", 2, 3),
            (b"CNG", b"III", b"1=1I1=", 2, 3),
            (b"CGN", b"III", b"2=1I", 2, 3),
            # 2 insertions, 1 base apart, next to no coverage (n=2)
            (b"NCNG", b"IIII", b"1I1=1I1=", 2, 3),
            (b"CNGN", b"IIII", b"1=1I1=1I", 2, 3),
            # 3 insertions, 1 base apart, next to no coverage (n=1)
            (b"NCNGN", b"IIIII", b"1I1=1I1=1I", 2, 3),
        ]
        self.assert_equal(ref, relvecs, expects)

    def test_insert_dangling_invalid(self):
        """ Test when the read has a dangling insertion (a 5' or 3'
        insertion without a corresponding 3' or 5' insertion). """
        ref = DNA(b"ACGT")
        relvecs = [
            [MINS5, MATCH, MATCH, MATCH],
            [MATCH, MINS5, MATCH, MATCH],
            [MATCH, MATCH, MINS5, MATCH],
            [MATCH, MINS3, MATCH, MATCH],
            [MATCH, MATCH, MINS3, MATCH],
            [MATCH, MATCH, MATCH, MINS3],
            [MINS5, MATCH, MINS3, MATCH],
            [MATCH, MINS5, MATCH, MINS3],
            [NOCOV, MINS5, MATCH, NOCOV],
            [NOCOV, MATCH, MINS3, NOCOV],
        ]
        self.assert_raise(ref, relvecs, ValueError,
                          "3' ins")

    def test_insert_bare_invalid(self):
        """ Test when the read has bare insertions (with no underlying
        relationship). """
        ref = DNA(b"ACGT")
        relvecs = [
            # 1 bare insertion
            [INS_3, MATCH, MATCH, MATCH],
            [MINS5, INS_3, MATCH, MATCH],
            [MATCH, MINS5, INS_3, MATCH],
            [MATCH, MATCH, MINS5, INS_3],
            [INS_5, MINS3, MATCH, MATCH],
            [MATCH, INS_5, MINS3, MATCH],
            [MATCH, MATCH, INS_5, MINS3],
            [MATCH, MATCH, MATCH, INS_5],
            [NOCOV, INS_3, MATCH, NOCOV],
            [NOCOV, MINS5, INS_3, NOCOV],
            [NOCOV, INS_5, MINS3, NOCOV],
            [NOCOV, MATCH, INS_5, NOCOV],
            # 2 bare insertions
            [INS_8, MINS3, MATCH, MATCH],
            [MINS5, INS_8, MINS3, MATCH],
            [MATCH, MINS5, INS_8, MINS3],
            [MATCH, MATCH, MINS5, INS_8],
            [NOCOV, INS_8, MINS3, NOCOV],
            [NOCOV, MINS5, INS_8, NOCOV],
        ]
        self.assert_raise(ref, relvecs, ValueError,
                          "Invalid relation 0")

    def test_insert_deletion_invalid(self):
        """ Test when the read has an insertion next to a deletion. """
        ref = DNA(b"ACGT")
        relvecs = [
            [INS_3 + DELET, MATCH, MATCH, MATCH],
            [MINS5, INS_3 + DELET, MATCH, MATCH],
            [MATCH, MINS5, INS_3 + DELET, MATCH],
            [MATCH, MATCH, MINS5, INS_3 + DELET],
            [DELET + INS_5, MINS3, MATCH, MATCH],
            [MATCH, DELET + INS_5, MINS3, MATCH],
            [MATCH, MATCH, DELET + INS_5, MINS3],
            [MATCH, MATCH, MATCH, DELET + INS_5],
        ]
        self.assert_raise(ref, relvecs, ValueError,
                          "Relation .+ is del and ins")

    def test_insert_non_match_valid(self):
        """ Test when the read has insertions next to substitutions or
        low-quality base calls. """
        ref = DNA(b"ACGT")
        relvecs = [
            # 1 insertion next to 1 substitution
            [INS_3 + SUB_C, MATCH, MATCH, MATCH],
            [SUB_C + INS_5, MINS3, MATCH, MATCH],
            [SUB_C, MINS5, MINS3, MATCH],
            [MINS3, SUB_T, MATCH, MATCH],
            [MINS5, INS_3 + SUB_T, MATCH, MATCH],
            [MATCH, SUB_T + INS_5, MINS3, MATCH],
            [MATCH, SUB_T, MINS5, MINS3],
            # 1 insertion next to 1 low-quality base call
            [INS_3 + ANY_N - SUB_A, MATCH, MATCH, MATCH],
            [ANY_N - SUB_A + INS_5, MINS3, MATCH, MATCH],
            [ANY_N - SUB_A, MINS5, MINS3, MATCH],
            [MINS3, ANY_N - SUB_C, MATCH, MATCH],
            [MINS5, INS_3 + ANY_N - SUB_C, MATCH, MATCH],
            [MATCH, ANY_N - SUB_C + INS_5, MINS3, MATCH],
            [MATCH, ANY_N - SUB_C, MINS5, MINS3],
        ]
        expects = [
            # 1 insertion next to 1 substitution
            (b"NCCGT", b"IIIII", b"1I1X3=", 1, 4),
            (b"CNCGT", b"IIIII", b"1X1I3=", 1, 4),
            (b"CCNGT", b"IIIII", b"1X1=1I2=", 1, 4),
            (b"NATGT", b"IIIII", b"1I1=1X2=", 1, 4),
            (b"ANTGT", b"IIIII", b"1=1I1X2=", 1, 4),
            (b"ATNGT", b"IIIII", b"1=1X1I2=", 1, 4),
            (b"ATGNT", b"IIIII", b"1=1X1=1I1=", 1, 4),
            # 1 insertion next to 1 low-quality base call
            (b"NNCGT", b"I!III", b"1I1M3=", 1, 4),
            (b"NNCGT", b"!IIII", b"1M1I3=", 1, 4),
            (b"NCNGT", b"!IIII", b"1M1=1I2=", 1, 4),
            (b"NANGT", b"II!II", b"1I1=1M2=", 1, 4),
            (b"ANNGT", b"II!II", b"1=1I1M2=", 1, 4),
            (b"ANNGT", b"I!III", b"1=1M1I2=", 1, 4),
            (b"ANGNT", b"I!III", b"1=1M1=1I1=", 1, 4),
        ]
        self.assert_equal(ref, relvecs, expects)


class TestRelAsSam(ut.TestCase):
    """ Test function `rel.as_sam`. """

    def test_line_in_sam_format(self):
        line = as_sam(b"FS10000136:77:BPG61616-2310:1:1101:1000:1300", 99,
                      "SARS2_FSE", 1, 42, b"151=", "=", 133, 283,
                      DNA(b"CCCTGTGGGTTTTACACTTAAAAACACAGTCTGTACCGTCTGCGGTATGTG"
                          b"GAAAGGTTATGGCTGTAGTTGTGATCAACTCCGCGAACCCATGCTTCAGTC"
                          b"AGCTGATGCACAATCGTTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCC"),
                      b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                      b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                      b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:")
        expect = (b"FS10000136:77:BPG61616-2310:1:1101:1000:1300\t99\tSARS2_FSE"
                  b"\t1\t42\t151=\t=\t133\t283\t"
                  b"CCCTGTGGGTTTTACACTTAAAAACACAGTCTGTACCGTCTGCGGTATGTGGAAAGGTT"
                  b"ATGGCTGTAGTTGTGATCAACTCCGCGAACCCATGCTTCAGTCAGCTGATGCACAATCG"
                  b"TTTTTAAACGGGTTTGCGGTGTAAGTGCAGCCC\t"
                  b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                  b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF"
                  b"FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:\n")
        self.assertEqual(line, expect)


# Module: seq ##########################################################

class TestSeqConstants(ut.TestCase):
    """ Test constants of `seq` module. """

    def test_bases(self):
        """ Test that `BASES` contains the four DNA letters. """
        self.assertEqual(BASES, b"ACGT")

    def test_rbase(self):
        """ Test that `RBASE` contains the four RNA letters. """
        self.assertEqual(RBASE, b"ACGU")

    def test_bases_arr(self):
        """ Test that `BASES` contains the four DNA ASCII codes. """
        self.assertTrue(isinstance(BASES_ARR, np.ndarray))
        self.assertIs(BASES_ARR.dtype.type, np.uint8)
        self.assertEqual(BASES_ARR.tolist(), [65, 67, 71, 84])


class TestSeqDna(ut.TestCase):
    """ Test class `DNA`. """

    def test_valid(self):
        """ Test whether valid DNA sequences can be created. """
        for length in range(1, 5):
            for bases in product(*(["ACGT"] * length)):
                dna = DNA("".join(bases).encode())
                self.assertEqual(len(dna), length)
                self.assertEqual(dna.decode(), "".join(bases))

    def test_slice(self):
        """ Test slicing DNA sequences. """
        dnaseq = "GATTACA"
        dna = DNA(dnaseq.encode())
        for i, j in combinations(range(len(dna) + 1), r=2):
            subseq = dna[i: j]
            self.assertTrue(isinstance(subseq, DNA))
            self.assertEqual(subseq.decode(), dnaseq[i: j])

    def test_reverse_complement(self):
        """ Test reverse complementing DNA sequences. """
        seqs = ["ACGT", "GTCAGCTGCATGCATG", "TAAAGTGGGGGGACATCATCATACT"]
        recs = ["ACGT", "CATGCATGCAGCTGAC", "AGTATGATGATGTCCCCCCACTTTA"]
        for seq, rec in zip(seqs, recs, strict=True):
            with self.subTest(seq=seq, rec=rec):
                fwd = DNA(seq.encode())
                rev = DNA(rec.encode())
                self.assertTrue(isinstance(fwd.rc, DNA))
                self.assertTrue(isinstance(rev.rc, DNA))
                self.assertEqual(fwd.rc, rev)
                self.assertEqual(rev.rc, fwd)
                self.assertEqual(fwd.rc.rc, fwd)
                self.assertEqual(rev.rc.rc, rev)

    def test_transcribe(self):
        """ Test transcribing DNA sequences. """
        dseqs = ["ACGT", "GTCAGCTGCATGCATG", "TAAAGTGGGGGGACATCATCATACT"]
        rseqs = ["ACGU", "GUCAGCUGCAUGCAUG", "UAAAGUGGGGGGACAUCAUCAUACU"]
        for dna, rna in zip(dseqs, rseqs, strict=True):
            with self.subTest(dna=dna, rna=rna):
                tr = DNA(dna.encode()).tr()
                self.assertTrue(isinstance(tr, RNA))
                self.assertEqual(tr, RNA(rna.encode()))

    def test_invalid_bases(self):
        """ Test whether invalid characters raise ValueError. """
        for char in printable:
            if char not in "ACGT":
                self.assertRaises(ValueError, DNA, char.encode())

    def test_zero(self):
        """ Test that zero-length DNA sequences raise ValueError. """
        self.assertRaises(ValueError, DNA, b"")


class TestSeqRna(ut.TestCase):
    """ Test class `RNA`. """

    def test_valid(self):
        """ Test whether valid RNA sequences can be created. """
        for length in range(1, 5):
            for bases in product(*(["ACGU"] * length)):
                rna = RNA("".join(bases).encode())
                self.assertEqual(len(rna), length)
                self.assertEqual(rna.decode(), "".join(bases))

    def test_slice(self):
        """ Test slicing RNA sequences. """
        rnaseq = "GAUUACA"
        rna = RNA(rnaseq.encode())
        for i, j in combinations(range(len(rna) + 1), r=2):
            subseq = rna[i: j]
            self.assertTrue(isinstance(subseq, RNA))
            self.assertEqual(subseq.decode(), rnaseq[i: j])

    def test_reverse_complement(self):
        """ Test reverse complementing RNA sequences. """
        seqs = ["ACGU", "GUCAGCUGCAUGCAUG", "UAAAGUGGGGGGACAUCAUCAUACU"]
        recs = ["ACGU", "CAUGCAUGCAGCUGAC", "AGUAUGAUGAUGUCCCCCCACUUUA"]
        for seq, rec in zip(seqs, recs, strict=True):
            with self.subTest(seq=seq, rec=rec):
                fwd = RNA(seq.encode())
                rev = RNA(rec.encode())
                self.assertTrue(isinstance(fwd.rc, RNA))
                self.assertTrue(isinstance(rev.rc, RNA))
                self.assertEqual(fwd.rc, rev)
                self.assertEqual(rev.rc, fwd)
                self.assertEqual(fwd.rc.rc, fwd)
                self.assertEqual(rev.rc.rc, rev)

    def test_reverse_transcribe(self):
        """ Test reverse transcribing RNA sequences. """
        rseqs = ["ACGU", "GUCAGCUGCAUGCAUG", "UAAAGUGGGGGGACAUCAUCAUACU"]
        dseqs = ["ACGT", "GTCAGCTGCATGCATG", "TAAAGTGGGGGGACATCATCATACT"]
        for rna, dna in zip(rseqs, dseqs, strict=True):
            with self.subTest(rna=rna, dna=dna):
                rt = RNA(rna.encode()).rt()
                self.assertTrue(isinstance(rt, DNA))
                self.assertEqual(rt, DNA(dna.encode()))

    def test_invalid_bases(self):
        """ Test whether invalid characters raise ValueError. """
        for char in printable:
            if char not in "ACGU":
                self.assertRaises(ValueError, RNA, char.encode())

    def test_zero(self):
        """ Test that zero-length RNA sequences raise ValueError. """
        self.assertRaises(ValueError, RNA, b"")


class TestSeqExpandDegenerateSeq(ut.TestCase):
    """ Test function `seq.expand_degenerate_seq`. """

    def test_zero_degenerate(self):
        """ Test that the original sequence is returned. """
        self.assertEqual(list(expand_degenerate_seq(b"ACGT")),
                         [b"ACGT"])

    def test_one_degenerate(self):
        """ Test that one sequence is returned for each DNA base. """
        self.assertEqual(list(expand_degenerate_seq(b"ACNT")),
                         [b"ACAT", b"ACCT", b"ACGT", b"ACTT"])

    def test_two_degenerate(self):
        """ Test that one sequence is returned for every combination of
        two DNA bases. """
        self.assertEqual(list(expand_degenerate_seq(b"NCGN")),
                         [b"ACGA", b"ACGC", b"ACGG", b"ACGT",
                          b"CCGA", b"CCGC", b"CCGG", b"CCGT",
                          b"GCGA", b"GCGC", b"GCGG", b"GCGT",
                          b"TCGA", b"TCGC", b"TCGG", b"TCGT"])


# Module: sim ##########################################################

class TestSimRandDna(ut.TestCase):
    """ Test function `sim.rand_dna`. """

    def test_type(self):
        """ Test that the type of the return value is DNA. """
        self.assertIs(type(rand_dna(1)), DNA)

    def test_length(self):
        """ Test that the length of the DNA sequence is as expected. """
        for length in range(1, 10):
            with self.subTest(length=length):
                self.assertEqual(len(rand_dna(length)), length)

    def test_invalid_length(self):
        """ Test that lengths ≤ 0 raise ValueError. """
        for length in range(0, -10, -1):
            with self.subTest(length=length):
                self.assertRaises(ValueError, rand_dna, length)


if __name__ == "__main__":
    ut.main()
