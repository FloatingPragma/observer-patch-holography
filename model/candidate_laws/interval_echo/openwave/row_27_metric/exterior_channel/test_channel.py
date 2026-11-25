"""Recompute the gate occupation by sparse actions on complete independent sharp starts."""

import json
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply

from . import measure


def test_integer_certificates_and_independent_gate_evolution():
    certified = measure.certified_gate()
    for occupied, (lo, hi) in enumerate(certified):
        pairs = np.arange(10)
        belt = np.zeros((10, 10))
        belt[pairs, (pairs + 1) % 10] = belt[(pairs + 1) % 10, pairs] = -0.5
        mediator = belt.copy()
        mediator[5, 5] = -32 * occupied
        start = expm_multiply(-12j * csr_matrix(mediator), np.eye(10)[0])
        joint = np.kron(mediator, np.eye(10)) + np.kron(np.eye(10), belt)
        joint[2, 2] -= 32
        out = expm_multiply(
            -4j * csr_matrix(joint), np.kron(start[:, None], np.eye(10))
        )
        answer = np.abs(out[2::10]) ** 2
        assert lo < answer.sum() / 10 < hi
        assert hi - lo < 3e-7
    assert certified[1][0] - certified[0][1] > 0.01543


def test_complete_uniform_census_controls_and_held_prediction():
    out = json.loads((measure.HERE / "results.json").read_text())
    fresh = measure.run()
    assert len(out["rows"]) == len(fresh["rows"]) == 80
    np.testing.assert_allclose(
        out["certified_gate_probability_empty_then_occupied"],
        fresh["certified_gate_probability_empty_then_occupied"],
        rtol=0,
        atol=0,
    )
    for experiment in (out, fresh):
        n = experiment["shots"]
        rows = experiment["rows"]
        for row in rows:
            bank = np.array([*row["primary_site_counts"], row["held_site_counts"]])
            assert np.all(bank.sum(axis=1) == n)
            estimate = (1 - bank[:2].sum(axis=0) / n) / 8
            residual = max(abs(bank[2] / n - estimate))
            assert residual == row["held_max_error"]
            assert residual < experiment["held_forecast_radius"]
            if row["k"] == 0 or row["tag"] == 2 or row["parked"]:
                assert (
                    np.max(abs(bank / n - 0.1))
                    < experiment["simultaneous_probability_radius"]
                )
        primary = {
            (row["time"], row["pole"]): row
            for row in rows
            if row["tag"] == 1 and row["k"] == 8 and not row["parked"]
        }
        eps = experiment["simultaneous_probability_radius"]
        for pole in (0, 1):
            gate = primary[4, pole]["primary_site_counts"][0][2] / n
            lo, hi = experiment["certified_gate_probability_empty_then_occupied"][pole]
            assert lo - eps < gate < hi + eps
        for arm in (0, 1):
            gap = (
                primary[2, 1]["primary_site_counts"][arm][2]
                - primary[2, 0]["primary_site_counts"][arm][2]
            ) / n
            assert gap < -2 * eps
        held_gap = (
            primary[2, 1]["held_site_counts"][2] - primary[2, 0]["held_site_counts"][2]
        ) / n
        assert held_gap > 2 * eps
