"""Check the native sector, full observer read, equilibrium law and negative controls."""

import json
from math import comb, log
from pathlib import Path
import unittest

import numpy as np

from model.executable_model import federation
from .measure import apparatus, capture_labels, generator, repeated_endpoints, sampling_kernel, seam


class RegionalEquilibriumTest(unittest.TestCase):
    def test_actual_native_sector_and_fresh_word_meters(self):
        np.testing.assert_allclose(seam(), [[0, 0, 0, 0], [0, 32.5, -.5, 0],
                                            [0, -.5, 32.5, 0], [0, 0, 0, 0]])
        car = federation(0)
        ports, contacts, cut, words = apparatus(car, [0])
        labels, k, group_k, cut_ids, cut_words, gates = capture_labels(car, [0], ports, cut, words)
        self.assertEqual(gates["echo_read_undo"], 12*len(words))
        self.assertEqual(len(set(cut_ids)), len(cut_words))
        for j in range(3):
            self.assertTrue(np.all(np.bincount(labels)[group_k == j] == comb(len(cut), j)))
        p, _ = sampling_kernel(generator(car, words, contacts))
        for tilted in (False, True):
            endpoints, diagnostics = repeated_endpoints(p, k, tilted)
            weight = np.exp2(-k.astype(float)) if tilted else np.ones(len(k))
            self.assertGreater(diagnostics["spectral_gap"], 0)
            for row in endpoints:
                np.testing.assert_allclose(row, weight/weight.sum(), atol=1e-12)

    def test_archived_held_transfers_and_competing_estimators(self):
        r = json.loads(Path(__file__).with_name("results.json").read_text())
        self.assertEqual([x["board_cells"] for x in r["observations"]], [20]*3+[80]*3+[320]*3)
        self.assertAlmostEqual(r["fitted_nats_per_cut_only_bit"], log(2), delta=.003)
        for x in r["observations"]:
            m = len(x["cut_ports"])
            self.assertEqual(m, 2+len(x["cells"]))
            self.assertEqual(x["enrolled_words"], comb(12*len(x["cells"]), 2))
            direction = np.array(x["mean_energy_matched_transfer"])
            self.assertAlmostEqual(sum(direction), 0)
            self.assertAlmostEqual(direction @ x["class_mean_native_energy"], 0)
            self.assertAlmostEqual(x["unmonitored"]["purity"], 1, places=12)
            for mode, experiment in x["repair"].items():
                self.assertLess(experiment["worst_start_log10_tv_upper"], -10)
                data = experiment["acquisitions"]["holdout"]
                for name in ["interior_start_bulk_counts", "cut_start_edge_counts", "cut_only_counts"]:
                    self.assertEqual(sum(data[name]), r["preparations_per_meter"])
                for s in data["slopes"].values():
                    prediction = 0 if mode == "unbiased" else log(2)*s["occupied_cut_endpoints"]
                    self.assertAlmostEqual(s["cut_only_residual_nats"], prediction, delta=.05)
            fair = x["repair"]["unbiased"]["acquisitions"]["holdout"]["slopes"]
            self.assertGreater(abs(fair["matched_mean_energy"]["occupation_residual_nats"]), .4)
            constrained = fair["matched_mean_energy_and_cut_load"]
            self.assertEqual(constrained["occupied_cut_endpoints"], 0)
            self.assertGreater(constrained["bulk_nats"], .4)
            self.assertLess(abs(constrained["cut_only_residual_nats"]), .004)
            biased = x["repair"]["outward_rejection"]["acquisitions"]["holdout"]["slopes"]
            self.assertGreater(biased["0_to_2"]["cut_only_residual_nats"], 1.35)
