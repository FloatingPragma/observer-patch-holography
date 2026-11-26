"""Test the full seam columns, detector gates and observed passive-energy transfer."""

import json
import unittest

import numpy as np
from scipy.linalg import eigh

from . import calibration, measure
from .material import columns, hardware
from .readers import circuit


class PassiveDecay(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record = measure.run()
        cls.cal = calibration.run()

    def test_all_registered_receipts(self):
        for filename, current in (
            ("records.json", self.record),
            ("calibration.json", self.cal),
        ):
            saved = json.loads((measure.HERE / filename).read_text())
            measure.compare_records(saved, current)

    def test_complete_geometry_native_heads_and_cut(self):
        board, dimers, heads = hardware()
        s, r, uncoupled, join = measure.hamiltonians(len(dimers))
        for attached in (True, False):
            stock = dict(heads)
            if not attached:
                a, c = sorted([dimers[0][0], dimers[2][0]])
                stock[a, c] = 32
            native = measure.physical(board, dimers, stock)
            tested = np.diag(uncoupled + attached * join)
            np.testing.assert_allclose(
                native - native[0], tested - tested[0], atol=1e-12
            )
        self.assertLess(columns(), 1e-12)
        self.assertLess(columns(comparison=True), 1e-12)
        for component, cost in ((s, 32), (r, 33)):
            e, v = eigh(component)
            compiled = circuit(v, cost)
            np.testing.assert_allclose(
                compiled.T @ component @ compiled, np.diag(e), atol=1e-12
            )

    def test_light_formation_and_energy(self):
        r = self.record
        self.assertLess(r["reservoir_levels"][1], r["light_cut"])
        self.assertLess(r["joined_levels"][1], r["light_cut"])
        self.assertLess(r["cold_band"]["all_time_cold_band_amplitude_bound"], 1e-8)
        for arm in r["arms"]:
            for t in arm["settings"]:
                self.assertEqual(np.sum(t["joint"]), measure.SHOTS)
                self.assertAlmostEqual(
                    t["total_energy"], arm["initial_energy"]["mean"], places=9
                )
        for meter in self.cal["signed_receivers"]:
            self.assertEqual(meter["coherent_amplitude_error"], 0)
            for i, target in ((0, 7), (1, 0), (3, 0), (4, 1)):
                self.assertEqual(meter["word_00_to_22"][i], target)

    def test_simultaneous_registered_contrasts(self):
        arm = {
            (a["source"], a["joined"]): a["settings"][5] for a in self.record["arms"]
        }
        self.assertTrue(all(v["wait"] == 16 for v in arm.values()))
        radius = np.sqrt(np.log(2 * 2048 / 0.01) / (2 * measure.SHOTS))
        raw = arm["parent", True]
        for control in (arm["ground", True], arm["parent", False]):
            self.assertGreater(
                (raw["joint"][2][1] - control["joint"][2][1]) / measure.SHOTS
                - 2 * radius,
                0.01,
            )
            self.assertGreater(
                (raw["far_minus"] - control["far_minus"]) / measure.SHOTS - 2 * radius,
                0.01,
            )
        changed = self.cal["changed_law"]
        for i in range(len(measure.WAITS)):
            self.assertLess(
                abs(
                    changed[0]["settings"][i]["far_minus"]
                    - changed[1]["settings"][i]["far_minus"]
                )
                / measure.SHOTS,
                2 * radius,
            )
