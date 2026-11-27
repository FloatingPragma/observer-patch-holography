"""Check preparation, exact range equivalence, entropy slopes and held counts."""

from fractions import Fraction
from itertools import product
import json
import math

import numpy as np

from . import measure


def test_writer_is_one_positive_three_tag_lottery():
    _, _, _, source = measure.sectors.volume.radar_source_and_hand(
        measure.clock.federation(0)
    )
    lower, middle, upper = np.asarray(source["source_units"]).reshape(3, 9)
    condition = {-1: lower, 1: upper}
    q = np.outer(middle, middle)
    delta = np.outer(upper - middle, upper - middle)
    for s in (Fraction(-1, 2), Fraction(1, 4), Fraction(1, 2)):
        weights = {
            (a, b, c): (1 + s * a * (b - c)) / 8
            for a, b, c in product((-1, 1), repeat=3)
        }
        assert min(weights.values()) >= 0 and sum(weights.values()) == 1
        for i, j, rho in [(0, 1, s), (0, 2, -s), (1, 2, 0), (0, 0, 1)]:
            tab = sum(
                w * np.outer(condition[z[i]], condition[z[j]])
                for z, w in weights.items()
            )
            np.testing.assert_array_equal(tab, q + rho * delta)


def test_exact_full_board_ranges_and_addressed_entropy():
    assert measure.atlas_reflection() == [0, 2, 1, 4, 3, 5, 7, 6, 8, 9, 11, 10]
    height, bins, rows = measure.native_regions()
    assert (height, bins) == (23, 84)
    groups = measure.word_entropy.probabilities()
    base, conditional = measure.exact_derivative(groups)
    assert base == Fraction(947661, 1183693)
    # Compute the log derivative from ungrouped physical writer spectra.
    _, _, _, native = measure.sectors.volume.radar_source_and_hand(
        measure.clock.federation(0)
    )
    laws, _ = measure.sectors.shared_variance.source_words(native)
    reverse, q, _ = laws
    d = q - reverse
    a, b = q + d / 4, q - d / 4
    derivative = np.sum(d * (np.log(b) - np.log(a)))
    assert math.isclose(derivative, math.log(float(base)) / 144, rel_tol=3e-14)
    assert measure.word_entropy.endpoints(conditional)[1] < 0
    assert [pop for _, pop, _ in rows] == [
        [4, 1, 1],
        [6, 2, 2],
        [9, 4, 4],
        [12, 6, 6],
        [16, 9, 9],
    ]
    for _, (a, b, _), p in rows:
        assert p == Fraction(2 * a * b, 320 * 319)


def test_held_counts_and_prediction_at_all_five_radii():
    out = json.loads((measure.HERE / "results.json").read_text())
    assert out == measure.run()
    for bank in out["terminal_outcomes"].values():
        assert (
            sum(sum(row) for row in bank["first_admitted_shell_AB_AC"]) + bank["other"]
            == measure.TRIALS
        )
    gain = math.log(947661 / 1183693) / 8
    for row in out["fixed_board_nested_readings"]:
        s = row["entropy_slope_per_fourth"]
        e = float(Fraction(row["angular_energy_slope_per_fourth"]))
        si, ei = row["addressed_entropy_99"], row["angular_energy_score_99"]
        assert si[0] <= s[0] <= s[1] <= si[1] < 0 < ei[0] <= e <= ei[1]
        assert math.isclose(s[0], gain * e, rel_tol=1e-14)
        assert max(si[0], gain * ei[1]) < min(si[1], gain * ei[0])
    assert len(out["midpoint_readings"]) == 4
    assert out["angular_energy_slope_per_fourth_limit"] == "9/400"
    lim = out["entropy_slope_per_fourth_limit"]
    assert lim[0] < -0.0006254933847632132 < lim[1]
    base = Fraction(out["conditional_entropy_derivative_log_base_per_144"])
    z = (1 - base) / (1 + base)
    assert Fraction(9, 1600) * (z + z**3 / 3) > Fraction(1, 1600)
