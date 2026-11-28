from fractions import Fraction as Q
import json

import numpy as np

from . import kernel_budget as budget
from . import kernel_certificate as exact
from . import newton_transport as newton


def test_density_limit_at_corners_and_continuum_chord():
    edges = np.array([0.0, 1.0, 3.0])
    weights = np.array([0.25, 0.75])
    chart = edges, weights
    table = exact.AnnularMap(chart)
    offsets = np.array([0.0, 1e-11, 0.2, 0.9, 2.0, 4.0, 7.0])
    for r in (0.7, 1.0, 1.5, 3.0):
        ref = [float(table.density(Q(str(r)), Q(str(h)))) for h in offsets]
        np.testing.assert_allclose(
            4 * np.pi * budget.density(r, offsets, chart), ref, atol=3e-14
        )
    rings = np.array([0.83, 1.67, 2.77])
    knots = budget.corners(rings, edges)
    a, b = knots[1:-1], knots[2:]
    mid = 2 * a * b / (a + b)

    def matrix(h):
        return budget.bank(dict(radius_kpc=rings), chart, h, 1.0)

    np.testing.assert_allclose(matrix(mid), (matrix(a) + matrix(b)) / 2, atol=2e-12)


def test_exact_local_witness_and_separation():
    inputs, charts, _, masses = budget.load_inputs()
    saved = json.loads((budget.HERE / "kernel_budget.json").read_bytes())
    n = "NGC2841"
    for name, fit in saved["individual"][n]["fits"].items():
        witness = exact.witness(inputs[n], charts[n], masses[n], fit)
        separation = exact.dual(
            inputs[n], charts[n], masses[n], fit["dual_ray"], name == "signed"
        )
        assert witness["verified_rings"] == len(inputs[n]["observed_kms"])
        assert 0 < witness["upper"] - separation["lower"] < 5
    fit = budget.local_inverse(
        np.array([[0.0, -1.0], [1.0, 0.0]]),
        np.array([1.0, -2.0]),
        np.array([2.0, 2.0]),
        False,
    )
    assert not fit["feasible"] and fit["dual_ray"] == ["1", "0"]


def test_newtonian_enclosure_including_origin_and_cubic_corners():
    edges = [Q(0), Q(1), Q(4)]
    seed = [Q(1, 4), Q(3, 4)]
    radius = Q(2)
    events = newton.cubic_events(radius, edges, seed)
    coeff = [Q(0)] * 4
    knots = sorted(events)
    chart = np.array(edges, float), np.array(seed, float)

    def fraction(x):
        return sum(
            w * newton.enclosed(radius, x, a, b)
            for w, a, b in zip(seed, edges[:-1], edges[1:], strict=True)
        )

    for lo, hi in zip(knots, knots[1:]):
        coeff = [a + b for a, b in zip(coeff, events[lo], strict=True)]
        bound = newton.max_ratio(coeff, lo, hi)
        for k in range(12):
            x = lo + (hi - lo) * k / 12
            assert fraction(x) <= bound
            assert x * fraction(x) == sum(v * x**i for i, v in enumerate(coeff))
            given = newton.cdf_matrix(
                [float(radius)], np.array([0.0, float(x if x else hi)]), chart
            )
            np.testing.assert_allclose(
                given[0, 0 if x == 0 else 1], float(fraction(x)), atol=2e-15
            )


def test_held_rotation_speeds_do_not_change_learned_profile():
    inputs = {
        k: dict(
            radius_kpc=[1.0, 2.0],
            observed_kms=[4.0, 5.0],
            error_kms=[0.1, 0.3],
            baryons_kms=[1.0, 1.0],
        )
        for k in ("a", "b", "c")
    }
    banks = {k: np.array([[1.0, 2.0], [2.0, 1.0]]) for k in inputs}
    for signed in (False, True):
        before = budget.train(banks, inputs, "c", signed)
        inputs["c"]["observed_kms"] = [1000.0, 5000.0]
        after = budget.train(banks, inputs, "c", signed)
        assert before == after
        inputs["c"]["observed_kms"] = [4.0, 5.0]
