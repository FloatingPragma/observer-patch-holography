"""Check the physical route, complete registered readout and false-loss discrimination."""

import json
import unittest
import numpy as np
from scipy.linalg import expm
from scipy.stats import beta
from ..passive_decay.measure import compare_records, SHOTS
from ..passive_decay.readers import charge
from . import joint, measure, precision
from .material import route


def interval(count, family):
    side = 0.01 / (2 * family)
    return (
        0.0 if count == 0 else beta.ppf(side, count, SHOTS - count + 1),
        1.0 if count == SHOTS else beta.ppf(1 - side, count + 1, SHOTS - count),
    )


def bank(name):
    return json.loads((measure.HERE / (name + ".json")).read_text())


class ExtendedDecay(unittest.TestCase):
    def test_routes_and_local_controls(self):
        board, d, heads, _, geometry = route((3,), True)
        self.assertEqual(geometry["raw_seams_to_each_end"], [4, 4])
        self.assertFalse(
            {p // 12 for pair in d[:2] for p in pair} & {p // 12 for p in d[-1]}
        )
        self.assertNotEqual(d[2][0] // 12, d[2][1] // 12)
        h = measure.rebuild(board, d, heads, False, True)
        cut = measure.rebuild(board, d, heads, True, True)
        words = np.arange(32)
        np.testing.assert_array_equal(
            np.diag(cut - h), ((words >> 1) & 1) * (words & 1)
        )
        self.assertEqual(
            charge(board, (*d[:2], d[-1]))[0]["word_00_to_22"][:5], [7, 0, 3, 0, 1]
        )
        # A direct material exponential factorizes at the interior cut for all
        # source and third-party reference inputs. The far dimer stays on X.
        last = np.array([[0.0, -0.5], [-0.5, 0.0]])
        np.testing.assert_allclose(
            expm(-0.63j * cut),
            np.kron(expm(-0.63j * cut[::2, ::2]), expm(-0.63j * last)),
            atol=1e-12,
        )

    def test_replay(self):
        compare_records(bank("records"), measure.run())
        for name, params in precision.CASES.items():
            compare_records(bank(name), measure.run(*params))
        compare_records(bank("joint"), joint.run())

    def test_finite_aperture_and_folds(self):
        folded = bank("records")["panels"]
        self.assertLess(
            folded[0]["energy"][1]["loss_gain_other"][0] / SHOTS,
            folded[0]["false_loss_bound"],
        )
        for p, q in zip(folded[0]["flight"], folded[1]["flight"], strict=True):
            self.assertLess(
                abs(
                    sum(r[1] for r in p["core_by_far_x"])
                    - sum(r[1] for r in q["core_by_far_x"])
                )
                / SHOTS,
                0.01,
            )
        for name in ("four", "three", "crossing"):
            source, idle, cut = bank(name)["panels"]
            for idx in range(2):

                def loss(panel):
                    return panel["energy"][idx]["loss_gain_other"][0]

                self.assertGreater(
                    interval(loss(source), 12)[0],
                    max(
                        source["false_loss_bound"],
                        interval(loss(idle), 12)[1],
                        interval(loss(cut), 12)[1],
                    ),
                )
        positive, idle, cut = bank("crossing")["panels"]

        def far(panel):
            return sum(row[1] for row in panel["flight"][5]["core_by_far_x"])

        radius = np.sqrt(np.log(2 * 100 / 0.01) / (2 * SHOTS))
        self.assertGreater((far(positive) - far(idle)) / SHOTS - 2 * radius, 0.01)
        self.assertEqual(far(cut), 0)

    def test_same_launch_and_changed_law(self):
        out = bank("joint")
        src, idle, cut = out["joint"]
        for arm in (src, idle, cut):
            self.assertEqual(
                sum(arm["loss_core"]) + sum(arm["gain_core"]) + arm["other"], SHOTS
            )
        hits = sum(src["loss_core"][1:3])
        self.assertGreater(
            interval(hits, 6)[0],
            max(
                bank("crossing")["panels"][0]["false_loss_bound"],
                interval(sum(idle["loss_core"][1:3]), 6)[1],
                interval(sum(cut["loss_core"][1:3]), 6)[1],
            ),
        )
        a, b = out["changed"]
        for key in ("loss", "far_minus"):
            lo, hi = interval(a[key], 6)
            other_lo, other_hi = interval(b[key], 6)
            self.assertLess(max(hi - other_lo, other_hi - lo), 0.01)
