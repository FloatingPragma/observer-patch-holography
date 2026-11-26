import hashlib
import json

import numpy as np

from . import measure
from .readout import examine


def ratio(numerator, denominator):
    assert denominator[0] * denominator[1] > 0
    ends = [n / d for n in numerator for d in denominator]
    return min(ends), max(ends)


def test_source_and_all_native_frequencies():
    path = measure.HERE / "results.json"
    assert (
        hashlib.sha256(path.read_bytes()).hexdigest()
        == "3c2f2d23c6dcfefaaa98ddd0ad4cfab7ec063248d4e7e63c4276ee7c492f2d66"
    )
    data = json.loads(path.read_text())
    measure.verify(data)
    for q in (-1, 1):
        for phase in measure.PHASES:
            z = measure.prepare(q, phase)
            np.testing.assert_allclose(
                z[1] / z[0], np.exp(1j * np.pi * phase / 4), atol=3e-14
            )


def test_held_phases_stock_interpolation_controls_and_material_gain():
    data = json.loads((measure.HERE / "results.json").read_text())
    readings = examine(data)
    assert readings["baselines"] == [[0, 1], [0, 38], [0, 39]]
    assert readings["four_tick_interval"][0] < 4 < readings["four_tick_interval"][1]
    shifts = readings["phase_changes_radians"]
    for stock in ("formed", "half"):
        ratios = [
            ratio(a, b)
            for a, b in zip(shifts["1"][stock], shifts["-1"][stock], strict=True)
        ]
        assert 0 < ratios[0][0] < ratios[0][1] < ratios[1][0] < ratios[1][1]
    for k in (0, 1):
        equivalent = [
            ratio(shifts[q]["formed"][k], shifts[q]["ordinary"][k]) for q in ("-1", "1")
        ]
        assert (
            equivalent[0][1] < equivalent[1][0] or equivalent[1][1] < equivalent[0][0]
        )
    for q in ("1", "-1"):
        for s in ("formed", "half", "ordinary"):
            np.testing.assert_allclose(
                shifts[q][s],
                data["readout"]["phase_changes_radians"][q][s],
                atol=1e-13,
                rtol=0,
            )


def test_certified_local_phase_response():
    data = json.loads((measure.HERE / "certificate.json").read_text())
    slopes = data["phase_slopes_radians_per_pairing"]
    assert data["minimum_probability_margin"] > 3e-11
    for s in ("product", "formed"):
        r1 = ratio(slopes["1"][s][0], slopes["-1"][s][0])
        r2 = ratio(slopes["1"][s][1], slopes["-1"][s][1])
        assert 0.11 < r1[0] <= r1[1] < 0.12 < 0.28 < r2[0] <= r2[1] < 0.32


def test_shutter_uncertainty_preserves_material_discrepancy():
    from . import timing

    clocks = json.loads((measure.HERE / "timing.json").read_text())
    eta = timing.audit(clocks)["per_probability_error"]
    assert eta < 2.55e-13
    raw = json.loads((measure.HERE / "results.json").read_text())
    changes = examine(raw, eta)["phase_changes_radians"]
    for k, lower, upper in ((0, 0.104, 0.127), (1, 0.274, 0.319)):
        gain = ratio(changes["1"]["formed"][k], changes["-1"]["formed"][k])
        assert lower < gain[0] < gain[1] < upper


def test_variable_wavenumber_on_the_same_two_source_directions():
    from mpmath import iv

    from . import dispersion

    receipt = json.loads((measure.HERE / "results.json").read_text())
    empirical = dispersion.counts_certificate(receipt, 2.55e-13)["determinant"]
    exact = json.loads((measure.HERE / "certificate.json").read_text())[
        "variable_wavenumber"
    ]
    assert 3.94e-5 < empirical[0] < exact["native"]["determinant"][0]
    assert exact["native"]["determinant"][1] < empirical[1] < 6.28e-5
    # A genuinely variable common-length model must pass the necessary test.
    # F=1 at the reference, so Im(F_destination/F_reference) is its phase slope.
    vectors = {}
    for q, k, dk in ((-1, 2, [1, 3]), (1, 7, [4, -1])):
        base = [0, 1, -1, 5]
        tangents = ([0, -1, 3, 8], [0, 2, -5, 1])
        vectors[q] = {"formed": [iv.mpc(1) for _ in base]}
        for name, d, slope in zip(("product", "ordinary"), tangents, dk, strict=True):
            vectors[q][name] = [
                iv.mpc(1, slope * v + k * u) for v, u in zip(base, d, strict=True)
            ]
    low, high = dispersion.rule_out(vectors)["determinant"]
    assert low <= 0 <= high
