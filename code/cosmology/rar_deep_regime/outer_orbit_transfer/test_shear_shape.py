"""Independent covariance and active-face tests for the exterior-orbit shear check."""

import itertools
import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from cosmology.rar_deep_regime.outer_orbit_transfer import shear_shape


def min_over_box_faces(y, covariance, lower, upper):
    """Enumerate the 3^6 faces and solve the quadratics on their affine hulls."""
    precision = np.linalg.inv(covariance)
    best = np.inf
    for face in itertools.product((-1, 0, 1), repeat=len(y)):
        fixed = np.asarray(face) != 0
        coords = np.where(np.asarray(face) == -1, lower, upper)
        free = ~fixed
        if free.any():
            lhs = precision[np.ix_(free, free)]
            rhs = precision[np.ix_(free, fixed)] @ (coords[fixed] - y[fixed])
            coords[free] = y[free] - np.linalg.solve(lhs, rhs)
        if np.all(coords >= lower - 1e-10) and np.all(coords <= upper + 1e-10):
            best = min(best, (coords - y) @ precision @ (coords - y))
    return best


@pytest.mark.parametrize("name", ["KiDS-bright", "GAMA"])
def test_forward_fit_and_face_enumeration(name):
    raw = json.loads(shear_shape.lensing_profile.RECEIPT.read_text())["profiles"][name]
    low, centre, high, signal, covariance = shear_shape.six_bins(raw)
    receipt = json.loads((shear_shape.HERE / "shear_shape.json").read_text())[
        "surveys"
    ][name]
    replay = shear_shape.run()["surveys"][name]
    central = receipt["central_bin_best_fit"]
    offset = shear_shape.shear(centre, 0)
    slope = np.sqrt(centre) / shear_shape.lensing_profile.FOUR_G
    amplitude = (
        slope
        @ np.linalg.solve(covariance, signal - offset)
        / (slope @ np.linalg.solve(covariance, slope))
    )
    fit = signal - offset - amplitude * slope
    assert central["chi_square"] == pytest.approx(
        fit @ np.linalg.solve(covariance, fit)
    )
    assert central["a0_m_s2"] == pytest.approx(amplitude**2, abs=1e-23)
    tested = set()
    for group, witnesses in receipt["witnesses"].items():
        for label, row in witnesses.items():
            actual = replay["witnesses"][group][label]
            residual = signal - np.array(row["centres"]["predicted_msun_psqpc"])
            expected = residual @ np.linalg.solve(covariance, residual)
            assert row["centres"]["chi_square"] == pytest.approx(expected)
            assert actual["centres"]["chi_square"] == pytest.approx(expected)
            assert expected >= central["chi_square"] - 1e-12
            if row["a0_m_s2"] in tested:
                continue
            tested.add(row["a0_m_s2"])
            a = np.sqrt(row["a0_m_s2"])
            best = min_over_box_faces(
                signal,
                covariance,
                shear_shape.shear(low, a),
                shear_shape.shear(high, a),
            )
            for source in (row, actual):
                assert source["arbitrary_within_bin_weights"][
                    "min_squared_covariance_distance"
                ] == pytest.approx(best)
            assert best <= expected
