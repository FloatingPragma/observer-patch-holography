"""Independent raw-bin, complex quadrature, and held-alternative checks."""

import json

import numpy as np
import pytest

from . import apparatus as a
from . import measure as m
from .sampling import multinomial


def test_eight_equal_bins_at_the_single_rounding_boundary():
    rng = np.random.default_rng(a.SEED)
    # Two gain calibrations and both zero-doublet shells precede this bank.
    for weights in ([1, 1, 0], [3, 1, 0], [1, 1, 0], [1, 1, 0]):
        multinomial(rng, a.SHOTS, np.array(weights))
    # The first two bins leave 196799 trials. In the third bin (n+1)/6
    # has BTPE mode 32800; two separately rounded operations give 32799.
    weights = np.array([1] * 8 + [0])
    np.testing.assert_array_equal(
        multinomial(rng, a.SHOTS, weights),
        [32641, 32704, 32611, 33008, 32546, 32755, 32901, 32978, 0],
    )


@pytest.fixture(scope="module")
def bank():
    stored = json.loads((m.HERE / "results.json").read_text())
    m.compare(stored, json.loads(json.dumps(m.run())))
    return stored


def test_all_stations_controls_and_quartic(bank):
    geo = bank["geometry"]
    assert geo["owners"] == {"source": 0, "receiver": 13}
    assert geo["address_clicks"] == [
        [word, word, "010101"] for word in ("000101", "010001", "010100")
    ]
    assert len(bank["quartic"]) == 14
    for p in bank["quartic"]:
        raw = p["raw"]
        assert sum(n for i, j, k, n in raw) == a.SHOTS
        assert sum(n for i, j, k, n in p["blocked"] if j == 3) == a.SHOTS
        num = sum(
            n * (1, 2, 0)[k] * (-1) ** ((i == 2) + (j == 2)) for i, j, k, n in raw
        )
        den = sum(n for i, j, k, n in raw if i == j == 2)
        assert np.isclose(p["measured"], num / den)
        norm = sum(abs(complex(*v)) ** 2 for v in p["h"])
        prediction = 1.5 * ((1 + norm) if p["repeated_address"] else (norm - 1) ** 2)
        assert p["interval"][0] <= prediction <= p["interval"][1]
        if np.isclose(norm, 1):
            if p["repeated_address"]:
                assert p["interval"][0] > 2.8
            else:
                assert p["interval"][0] < 0 < p["interval"][1] < 0.1
    assert bank["vacuum_phase"]["X"] == [[1, a.SHOTS]]


def test_all_complex_component_returns(bank):
    panels = {p["name"]: p for p in bank["lines"]}
    for panel in panels.values():
        denominator = sum(
            n for i, j, k, n in panel["reference_raw"] if i == j == k == 2
        )
        assert (
            sum(row[-1] for row in panel["blocked_reference"] if row[1] == 3) == a.SHOTS
        )
        signed = np.zeros((2, 2), complex)
        for fair, i, j, axes, raw, blocked in panel["shots_by_setting"]:
            assert sum(r[-1] for r in raw) == a.SHOTS
            assert sum(r[-1] for r in blocked if r[1] == 3) == a.SHOTS
            signed_total = 0
            for p, q, r, g, n in raw:
                if p not in (i, 2) or q not in (j, 2) or r not in (0, 2) or g == 2:
                    continue
                negative = (p == 2) + (q == 2) + (r == 2) + (g == 1)
                signed_total += (-1) ** negative * n
            signed[fair, i] += 1j ** axes.count("Y") * signed_total / (8 * denominator)
        for row in panel["estimates"]:
            first, second = signed[int(row["fair_gain"])]
            ratio = (
                first - second if row["orientation"] == "epsilon" else first + second
            )
            assert abs(ratio - complex(*row["ratio"])) < 1e-12
            if row["orientation"] == "epsilon":
                expected = [0, 0] if row["fair_gain"] else panel["prediction"]
                assert all(
                    low <= z <= high
                    for (low, high), z in zip(row["intervals"], expected, strict=True)
                )
    assert panels["exchanged"]["estimates"][0]["intervals"][0][1] < -0.4
    assert panels["exchanged"]["estimates"][1]["intervals"][0][0] > 0.4
    assert panels["parallel"]["estimates"][1]["intervals"][1][0] > 0.4
    assert panels["held"]["estimates"][0]["intervals"][1][0] > 0.2
    assert panels["held"]["estimates"][1]["intervals"][1][1] < 0.2


def test_referenced_chart_and_continuous_oriented_read():
    def doublet(z):
        ref = a.read(z, site="receiver")[2]
        result = []
        for j in (0, 1):
            xy = []
            for axis in "XY":
                p = a.read(z, j, axis, "receiver")
                xy.append((p[j] - p[2]) / (2 * ref))
            result.append(complex(*xy))
        return np.array(result)

    psi = doublet((0.5 + 2j, -1.7j))
    h = doublet((-2 + 1j, 1 + 1j))
    np.testing.assert_allclose(psi, [0.5 + 2j, -1.7j], atol=2e-12)
    np.testing.assert_allclose(h, [-2 + 1j, 1 + 1j], atol=2e-12)
    u = a.scalar.execute(a.scalar.field_word((0.2, -0.3, 0.5, 0)))[:2, :2]
    np.testing.assert_allclose(
        np.linalg.det(np.stack([u @ psi, u @ h])),
        np.linalg.det(np.stack([psi, h])),
        atol=2e-12,
    )
    assert a.ratio_interval(0, 0, 0.01) is None
