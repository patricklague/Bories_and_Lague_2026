import unittest

import numpy as np

from get_pmf import block_mean_se, block_statistics


class BlockStatisticsTests(unittest.TestCase):
    def setUp(self):
        self.z = np.array([0.5, 10.5, 30.5, 31.5])
        self.pmfs = np.array([[2., 4., 8., 12.],
                              [-3., 0., 2., 4.],
                              [12., 10., 14., 16.]])

    def statistics(self, pmfs):
        return block_statistics(self.z, pmfs, 30., 10000, 2025)

    def test_independent_references_mean_and_exact_bootstrap_se(self):
        result = self.statistics(self.pmfs)
        references = self.pmfs[:, 2:].mean(axis=1)
        aligned = self.pmfs - references[:, None]
        mean, se = block_mean_se(aligned)
        np.testing.assert_allclose(result["block_references"], references)
        np.testing.assert_allclose(result["mean"], mean)
        np.testing.assert_allclose(result["se"], se, atol=1e-14)
        self.assertAlmostEqual(result["mean"][2:].mean(), 0.)
        self.assertEqual(result["method"], "exact (all 27 resamples)")
        raw_mean, raw_se = block_mean_se(self.pmfs)
        np.testing.assert_allclose(result["mean_raw"], raw_mean)
        np.testing.assert_allclose(result["se_raw"], raw_se)

    def test_independent_constant_offsets_do_not_change_final_mean_or_se(self):
        original = self.statistics(self.pmfs)
        shifted = self.statistics(self.pmfs + np.array([100., -30., 8.])[:, None])
        for key in ("mean", "se"):
            np.testing.assert_allclose(original[key], shifted[key], atol=1e-13)

    def test_common_water_support_excludes_incomplete_bin(self):
        pmfs = self.pmfs.copy()
        pmfs[0, 2] = np.nan
        result = self.statistics(pmfs)
        self.assertEqual(result["n_win"], 1)
        np.testing.assert_allclose(result["block_references"], pmfs[:, 3])
        self.assertTrue(np.isnan(result["mean"][2]))
        self.assertAlmostEqual(result["mean"][3], 0.)
        self.assertAlmostEqual(result["se"][3], 0.)

    def test_two_blocks_exact_bootstrap(self):
        result = self.statistics(self.pmfs[:2])
        aligned = self.pmfs[:2] - self.pmfs[:2, 2:].mean(axis=1)[:, None]
        np.testing.assert_allclose(result["se"], block_mean_se(aligned)[1])
        self.assertEqual(result["method"], "exact (all 4 resamples)")

    def test_missing_water_and_single_block_raise(self):
        with self.assertRaises(SystemExit):
            self.statistics(self.pmfs[:1])
        pmfs = self.pmfs.copy()
        pmfs[:, 2:] = np.nan
        with self.assertRaises(SystemExit):
            self.statistics(pmfs)


if __name__ == "__main__":
    unittest.main()
