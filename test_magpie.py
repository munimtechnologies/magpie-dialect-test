# Copyright 2026 MV2, Munim, Inc. Licensed under the Apache License, Version 2.0 (Apache-2.0).
# Made with AI agents working in MV2's name.
"""Unit tests for magpie_test.py (numpy; GLMM fitting tests need scipy and skip without it).
Run: python -m unittest -v test_magpie"""
import math
import unittest

import numpy as np

import magpie_test as mt

try:
    import scipy  # noqa: F401
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False


def row(f, g, t, a, b, c):
    return {"ID": f, "group": g, "call": t, "fl01": str(a), "group01": str(b), "exgroup01": str(c)}


class TestLoader(unittest.TestCase):
    def test_build(self):
        rows = [row("B", "X", "DS LH", 1, 0, 1), row("B", "X", "NL", 0, 1, 1),
                row("A", "Y", "DS  LH", 0, 1, 0), row("A", "Y", "NL", 1, 1, 0)]
        d = mt.build(rows)
        self.assertEqual(d["ids"], ["A", "B"])
        self.assertEqual(d["types"], ["DS LH", "NL"])
        self.assertEqual(d["groups"], ["X", "Y"])
        self.assertEqual(d["gidx"].tolist(), [1, 0])
        self.assertEqual(d["F"].tolist(), [[0, 1], [1, 0]])
        self.assertEqual(d["G"].tolist(), [[1, 1], [0, 1]])
        self.assertEqual(d["E"].tolist(), [[0, 0], [1, 1]])
        self.assertTrue(d["seen"].all())

    def test_build_rejects(self):
        with self.assertRaises(ValueError):
            mt.build([row("A", "X", "NL", 1, 0, 0), row("A", "X", "NL", 1, 0, 0)])
        with self.assertRaises(ValueError):
            mt.build([row("A", "X", "NL", "NA", 0, 0)])
        with self.assertRaises(ValueError):
            mt.build([row("A", "X", "NL", 1, 0, 0), row("A", "Y", "DS", 1, 0, 0)])

    def test_group_records(self):
        G = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        R = mt.group_records(G, np.array([0, 0, 1]), 2)
        self.assertEqual(R.tolist(), [[1, 1, 0], [0, 0, 1]])


class TestOverlap(unittest.TestCase):
    def test_overlap_and_d(self):
        F = np.array([[1, 1, 0, 0], [0, 0, 1, 1], [1, 0, 1, 0]])
        R = np.array([[1, 1, 0, 0], [0, 0, 1, 1], [1, 0, 0, 0]])
        M = mt.overlap_matrix(F, R)
        np.testing.assert_allclose(M, [[1.0, 0.0, 0.5], [0.0, 1.0, 0.0], [0.5, 0.5, 0.5]])
        self.assertAlmostEqual(mt.d_stat(M, np.array([0, 1, 2])), 1.75 / 3, places=12)

    def test_empty_repertoire_raises(self):
        with self.assertRaises(ValueError):
            mt.overlap_matrix(np.array([[0, 0]]), np.array([[1, 0]]))

    def test_perm_p(self):
        self.assertAlmostEqual(mt.perm_p(0.5, [0.1, 0.6, 0.5, 0.2]), 0.6)
        self.assertIsNone(mt.perm_p(None, [0.1]))

    def test_t1_detects_perfect_group_sets(self):
        R = np.zeros((4, 12), dtype=np.int64)
        for g in range(4):
            R[g, 3 * g:3 * g + 3] = 1
        assign = np.array([0, 0, 1, 1, 2, 2, 3, 3])
        F = R[assign]
        r = mt.t1_test(F, R, assign, 999, np.random.default_rng(1))
        self.assertAlmostEqual(r["D"], 1.0)
        self.assertLessEqual(r["p"], 0.05)
        self.assertEqual(r["reading"], "OWN_GROUP_EXCESS")


class TestOrder(unittest.TestCase):
    def test_bag_info(self):
        bag_id, multi = mt.bag_info(["A B", "B A", "A C"])
        self.assertEqual(bag_id.tolist(), [0, 0, 1])
        self.assertEqual(multi.tolist(), [True, True, False])

    def test_d_ord(self):
        types = ["A B", "B A", "A C"]
        bag_id, multi = mt.bag_info(types)
        F = np.array([[1, 0, 1], [0, 1, 0]])
        R = np.array([[1, 0, 0], [0, 1, 1]])
        A, B = mt.order_matrices(F, R, bag_id, multi)
        self.assertEqual(A.tolist(), [[1, 0], [0, 1]])
        self.assertEqual(B.tolist(), [[1, 1], [1, 1]])
        self.assertEqual(mt.d_ord(A, B, np.array([0, 1])), (1.0, 2))
        self.assertEqual(mt.d_ord(A, B, np.array([1, 0])), (-1.0, 2))

    def test_t2_uninformative_on_tiny_data(self):
        types = ["A B", "B A", "A C"]
        F = np.array([[1, 0, 1], [0, 1, 0]])
        R = np.array([[1, 0, 0], [0, 1, 1]])
        r = mt.t2_test(F, R, np.array([0, 1]), types, 19, np.random.default_rng(2))
        self.assertEqual(r["multi_order_types_produced"], 2)
        self.assertEqual(r["reading"], "UNINFORMATIVE")


class TestPlantAndDonor(unittest.TestCase):
    def test_plant_t1_keeps_size(self):
        F = np.array([[1, 0, 0, 0]])
        R = np.array([[0, 1, 1, 0]])
        Fp = mt.plant_t1(F, R, np.array([0]), 1.0, np.random.default_rng(3))
        self.assertEqual(int(Fp.sum()), 1)
        self.assertEqual(int(Fp[0, 0]), 0)
        self.assertEqual(int(Fp[0, 1] + Fp[0, 2]), 1)

    def test_plant_t2_keeps_size(self):
        types = ["A B", "B A", "A C"]
        bag_id, multi = mt.bag_info(types)
        F = np.array([[1, 0, 1]])
        R = np.array([[0, 1, 0]])
        Fp = mt.plant_t2(F, R, np.array([0]), bag_id, multi, np.random.default_rng(4))
        self.assertEqual(int(Fp.sum()), 2)
        self.assertEqual(int(Fp[0, 2]), 1)

    def test_donor_map(self):
        gidx = np.array([0, 0, 1, 2])
        p = mt.donor_map(gidx, np.random.default_rng(5))
        self.assertEqual(sorted(p.tolist()), [0, 1, 2, 3])
        self.assertTrue(np.all(gidx[p] != gidx))


class TestRanks(unittest.TestCase):
    def test_rankdata_ties(self):
        self.assertEqual(mt.rankdata([10, 20, 20, 5]).tolist(), [2.0, 3.5, 3.5, 1.0])

    def test_spearman(self):
        self.assertAlmostEqual(mt.spearman([1, 2, 3], [3, 2, 1]), -1.0)
        rho, p = mt.spearman_perm_p([1, 2, 3, 4, 5], [2, 4, 6, 8, 10], 9999,
                                    np.random.default_rng(6))
        self.assertAlmostEqual(rho, 1.0)
        self.assertLess(p, 0.05)


class TestGLMM(unittest.TestCase):
    def setUp(self):
        self.y = np.array([1, 0, 1, 1, 0, 0], dtype=float)
        self.x = np.array([1, 0, 1, 0, 0, 1], dtype=float)
        self.cl = np.array([0, 0, 1, 1, 2, 2])
        self.oh = np.zeros((3, 6))
        self.oh[self.cl, np.arange(6)] = 1.0
        z, w = np.polynomial.hermite.hermgauss(mt.GH_NODES)
        self.z, self.lw = z, np.log(w)

    def test_tiny_sigma_is_logistic(self):
        th = [0.2, 0.5, -6.0]
        eta = 0.2 + 0.5 * self.x
        logistic = float(np.sum(np.logaddexp(0, eta) - self.y * eta))
        got = mt.glmm_nll(th, self.y, self.x, self.oh, self.cl, self.z, self.lw)
        self.assertAlmostEqual(got, logistic, places=3)

    def test_large_sigma_matches_brute_force(self):
        s = 2.0
        th = [0.2, 0.5, math.log(s)]
        eta = 0.2 + 0.5 * self.x
        u = np.arange(-25.0, 25.0, 0.001)
        dens = np.exp(-u ** 2 / (2 * s * s)) / math.sqrt(2 * math.pi * s * s)
        total = 0.0
        for c in range(3):
            idx = self.cl == c
            lin = eta[idx][:, None] + u[None, :]
            ll = (self.y[idx][:, None] * lin - np.logaddexp(0, lin)).sum(axis=0)
            total += math.log(float(np.sum(np.exp(ll) * dens) * 0.001))
        got = mt.glmm_nll(th, self.y, self.x, self.oh, self.cl, self.z, self.lw)
        self.assertAlmostEqual(got, -total, places=3)

    @unittest.skipUnless(HAVE_SCIPY, "scipy not installed")
    def test_fit_recovers_effect(self):
        rng = np.random.default_rng(7)
        n_c, per = 10, 117
        cl = np.repeat(np.arange(n_c), per)
        x = (rng.random(n_c * per) < 0.3).astype(float)
        u = rng.normal(0, 0.5, n_c)
        y = (rng.random(n_c * per) < 1 / (1 + np.exp(-(-1.5 + 2.0 * x + u[cl])))).astype(float)
        r = mt.glmm_fit(y, x, cl, n_c)
        self.assertGreater(r["b1"], 1.0)
        self.assertLess(r["b1"], 3.0)
        self.assertLess(r["p"], 0.001)
        base = mt.glmm_fit(y, None, cl, n_c)
        self.assertLess(r["aicc"], base["aicc"])


if __name__ == "__main__":
    unittest.main()
