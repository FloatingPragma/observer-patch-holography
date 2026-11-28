"""Independent geometry, transfer, and observational checks on the published experiments."""

import json
from itertools import combinations

import numpy as np

from openwave.interval_echo.row_19.selected_collar import (
    finite_dose_curvature as survey,
)
from openwave.interval_echo.row_19.selected_collar import midpoint_ranging as midpoint
from openwave.interval_echo.row_19.selected_collar import (
    link_source_contraction as source_decoder,
)


def record(name):
    return json.loads((survey.lab.HERE / (name + ".json")).read_bytes())


def euclidean_curvature(ranges, tri):
    # Reconstruct coordinates and project the two opposite vertices into the
    # plane perpendicular to each edge. This uses neither the angle formula
    # nor the volume derivative used by the interval calculation.
    squares = dict(
        zip(
            survey.EDGES, np.exp(2 * (survey.clock_map() @ np.log(ranges))), strict=True
        )
    )
    lengths = {e: np.sqrt(v) for e, v in squares.items()}
    total = 2 * np.pi * sum(lengths.values())
    volume = 0
    for tet in survey.TETRA:
        if len(set(tet) & set(tri)) != 2:
            continue
        a, *others = tet
        gram = np.array(
            [
                [
                    (
                        squares[a, v]
                        + squares[a, w]
                        - (0 if v == w else squares[tuple(sorted((v, w)))])
                    )
                    / 2
                    for w in others
                ]
                for v in others
            ]
        )
        coords = dict(
            zip(tet, np.vstack((np.zeros(3), np.linalg.cholesky(gram))), strict=True)
        )
        volume += np.sqrt(np.linalg.det(gram)) / 6
        for i, j in combinations(tet, 2):
            u = (coords[j] - coords[i]) / lengths[i, j]
            offsets = [coords[k] - coords[i] for k in tet if k not in (i, j)]
            v, w = [c - c.dot(u) * u for c in offsets]
            angle = np.arctan2(np.linalg.norm(np.cross(v, w)), v.dot(w))
            total -= lengths[i, j] * angle
    return total / np.cbrt(volume)


def test_scale_invariance_and_schlafli_derivative(monkeypatch):
    iv = survey.iv
    monkeypatch.setattr(iv, "dps", 32)
    ones = [iv.mpf(1)] * 30
    tetra_volume = np.sqrt(2) / 12
    expected = (30 * np.pi - 54 * np.arccos(1 / 3)) / np.cbrt(9 * tetra_volume)
    for r in survey.along(ones, ones, at_reference=True).values():
        assert r["d_normalized"].a <= 0 <= r["d_normalized"].b
        np.testing.assert_allclose(float(r["normalized"].mid), expected)
    rng = np.random.default_rng(1345)
    means = rng.uniform(0.97, 1.03, 30)
    slope = rng.uniform(-0.5, 0.5, 30)
    step = 1e-5
    by_filling = survey.along(
        list(map(iv.mpf, means)), list(map(iv.mpf, slope)), at_reference=True
    )
    for tri, result in by_filling.items():
        contrast = (
            euclidean_curvature(means + step * slope, tri)
            - euclidean_curvature(means - step * slope, tri)
        ) / (2 * step)
        np.testing.assert_allclose(
            float(result["d_normalized"].mid), contrast, rtol=2e-7, atol=1e-8
        )


def test_force_and_equal_entropy_controls_leave_no_common_calibration():
    rows = record("finite_dose_curvature")["fillings"]

    def sign(v):
        return 1 if v[0] > 0 else -1 if v[1] < 0 else 0

    assert len(rows) == 6
    assert all(sign(r["equal_entropy_first_variation"]["residual"]) for r in rows)
    assert all(sign(r["normalized"]["held_minus_equal_entropy_ordinary"]) for r in rows)
    assert (
        sum(
            sign(r["normalized"]["held_from_components"])
            * sign(r["normalized"]["ordinary_from_components"])
            == -1
            for r in rows
        )
        == 5
    )
    forces = [
        r["resolved_port_instrument"][str(2**76)]
        for r in record("finite_dose_material")["readings"]
        if r["channel"][1] == "marked_force"
    ]
    low, high = [f["ordinary_energy_equivalent"] for f in forces]
    assert 0 < low[0] < low[1] < high[0] < high[1]
    for f in forces:
        assert all(sum(h) == 2**76 for h in f["plus_minus_miss"].values())


def test_transfer_decoder_including_degenerate_simplex_vertex():
    axes = [
        {"independent": [1, 1], "compensator": [0, 0], "held": [0.25, 0.25]},
        {"independent": [0, 0], "compensator": [1, 1], "held": [0.5, 0.5]},
    ]
    known = source_decoder.infer(axes, [2, 3])
    np.testing.assert_allclose(known["mixture_intervals"], [[0.25, 0.25], [0.5, 0.5]])
    np.testing.assert_allclose(known["entropy_lower"], 2)
    for row in axes:
        row["held"] = [0, 0]
    zero = source_decoder.infer(axes, [2, 3])
    assert abs(zero["entropy_lower"]) < 1e-15
    measurement = record("link_source_contraction")
    prior = record("finite_pairing_dose")
    for value, (lo, hi) in zip(
        [prior["doses"][-1]["dose"], prior["selected_fraction"]],
        measurement["blind_source_bounds"]["mixture_intervals"],
        strict=True,
    ):
        assert lo <= value <= hi
    assert all(
        sum(h) == measurement["starts_per_transfer"]
        for r in measurement["links"]
        for h in r["plus_minus_miss"].values()
    )
    assert (
        measurement["blind_source_bounds"]["entropy_lower"]
        <= prior["held_entropy_gain"][1]
    )


def test_null_certificate_and_observable_links():
    d = record("curvature_source_null")
    chosen = d["strongest_prepared_test"]
    for ref, key in [
        ("prepared", chosen),
        ("exact_equilibrium_box", chosen + "_stationarity"),
    ]:
        cert = d["certificates"][ref][key]
        assert (
            cert["constraint_count"] == 15
            and cert["entropy_per_unit_true_tangent"][0] > 2.4e-7
        )
    assert d["source_contrast_from_counts"][0] > 0
    raw = d["interior_link_measurements"]
    contrasts = []
    for entry, propagated in zip(
        raw["links"], d["independently_propagated_material_contrasts"][8:], strict=True
    ):
        lo, hi = entry["contrast"]
        assert lo <= propagated <= hi
        assert all(
            sum(arm) == raw["exposure"] for arm in entry["plus_minus_miss_by_arm"]
        )
        contrasts.append((lo, hi))
    assert len(contrasts) == 10 and sum(lo > 0 or hi < 0 for lo, hi in contrasts) == 5


def test_noncommuting_projection_bound_and_midpoint_scope():
    rng = np.random.default_rng(2701)
    for _ in range(32):
        psi = rng.normal(size=5) + 1j * rng.normal(size=5)
        psi /= np.linalg.norm(psi)
        survivor = psi.copy()
        free_probs = []
        for _ in range(10):
            m = rng.normal(size=5) + 1j * rng.normal(size=5)
            m /= np.linalg.norm(m)
            free_probs.append(abs(np.vdot(m, psi)) ** 2)
            survivor -= m * np.vdot(m, survivor)
            assert 1 - np.vdot(survivor, survivor).real <= 4 * sum(free_probs) + 1e-12
        assert np.isclose(
            midpoint.lower_range(np.array(free_probs)[:, None], 0)[0],
            3 - np.minimum(1, 4 * np.cumsum(free_probs)).sum() / 4,
        )
    for arm in record("midpoint_ranging")["results"].values():
        assert len(arm["edges"]) == 15
        assert min(e["mixture_ratio_lower"] for e in arm["edges"]) > 1.8815
        assert sum(e["affine_excess_lower"] > 0 for e in arm["edges"]) == 3
