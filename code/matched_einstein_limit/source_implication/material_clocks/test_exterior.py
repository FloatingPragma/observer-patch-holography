import json
from pathlib import Path

import numpy as np

from . import exterior


def test_exact_exterior_clocks_and_one_receiver():
    report = exterior.report()
    recorded = json.loads(
        Path(exterior.__file__).with_name("exterior.json").read_text()
    )
    assert report["boards"] == recorded["boards"]
    assert [row["clock_count"] for row in report["boards"]] == [0, 24, 104]
    assert report["ordinary_parked_energy_slope"][0] > 0.29756
    assert report["unpair_parked_energy_slope"] == 0
    joint = report["one_receiver"]
    assert joint["simultaneous_99_percent_frequency_half_width"] < 0.03
    np.testing.assert_allclose(
        [r["differential_angular_frequency"] for r in joint["experiments"]],
        [0, -np.sqrt(5) / 4],
        atol=1e-13,
    )
    for row, ref in zip(
        joint["experiments"], recorded["one_receiver"]["experiments"], strict=True
    ):
        assert row["counts"] == ref["counts"]


def test_integrated_pressure_with_compact_conserved_stress():
    # A nonspherical mixed-sign spatial stress, compact in the zero collar.
    n = 17
    xyz = np.indices((n, n, n)) - n // 2
    x, y, z = xyz
    mark = ((abs(x) < 3) & (abs(y) < 3) & (abs(z) < 3)) * (x * y + z * z)

    def deriv(a, i):
        return np.roll(a, -1, i) - np.roll(a, 1, i)

    hess = [[deriv(deriv(mark, i), j) for j in range(3)] for i in range(3)]
    lap = sum(hess[i][i] for i in range(3))
    stress = [[hess[i][j] - (lap if i == j else 0) for j in range(3)] for i in range(3)]
    assert max(np.max(abs(a)) for row in stress for a in row) > 1
    for j in range(3):
        residual = sum(deriv(stress[i][j], i) for i in range(3))
        assert not residual.any()
        for k in range(3):
            assert int((xyz[k] * residual).sum()) == -2 * int(stress[k][j].sum()) == 0
    # A cut through the source has a nonvanishing surface term.
    assert np.any(np.cumsum(stress[0][0].sum(axis=(0, 2))) != 0)
