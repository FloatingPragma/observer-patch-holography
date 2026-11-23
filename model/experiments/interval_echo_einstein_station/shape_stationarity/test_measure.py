"""Cross-check native resident response, source entropy, and separate geometry surveys."""

import json
from math import comb, log, log2
from pathlib import Path
import unittest

from .measure import resident_tag_strength


class ShapeStationarityTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.run_data = json.loads(Path(__file__).with_name("results.json").read_text())

    def test_stationary_population_and_fresh_word_counts(self):
        r = self.run_data
        self.assertEqual(resident_tag_strength(), 1.)
        self.assertEqual(len(r["regions"]), 12)
        for row in r["regions"] + r["same_source_regions"]:
            size = 12*len(row["source_cells"])
            m = len(row["active_source_cut_ports"])
            total = comb(size, 2)
            expected = sum(comb(m, k)*comb(size-m, 2-k)*log2(comb(m, k))/total for k in range(3))
            self.assertAlmostEqual(row["population"]["cut_bits"], expected)
            self.assertAlmostEqual(row["population"]["bulk_nats"], log(total)-log(2)*expected)
            self.assertEqual(total, row["native_pair_states"])
            self.assertGreater(row["native_pair_hopping_links"], total)
            self.assertEqual(row["resident_tag_launch_volume"], row["accessible_volume"]-2)
            self.assertAlmostEqual(row["resident_tag_aperture"], row["aperture"]-2*m/size)
            self.assertLess(row["extrapolation_relative_error"], .02)
            for sample in row["acquisitions"].values():
                self.assertEqual(sum(sample["exterior_counts"]), r["shots_per_independent_meter"])
                self.assertEqual(sum(sample["bulk_counts"]), r["shots_per_independent_meter"])

    def test_frozen_coefficients_across_shapes_and_matter(self):
        r = self.run_data
        self.assertAlmostEqual(r["fitted_nats_per_exterior_bit"], log(2), delta=.01)
        self.assertAlmostEqual(r["shape_changes"][0]["population"]["bulk_nats"], -.18093014669, places=10)
        self.assertAlmostEqual(r["shape_changes"][1]["population"]["bulk_nats"], -.15955332269, places=10)
        for key in ("shape_changes", "same_source_shape_queries", "within_shape_changes"):
            for x in r[key]:
                self.assertLess(abs(x["held"]["cut_residual_nats"]), .01)
                if key == "shape_changes":
                    self.assertEqual(x["held"]["aperture"], 2)
                else:
                    for channel in ("aperture", "resident_tag", "blank_flight"):
                        self.assertGreater(x["held"][f"{channel}_residual_nats"], .05)
                        if key == "within_shape_changes":
                            self.assertEqual(x["held"][channel], 0)
