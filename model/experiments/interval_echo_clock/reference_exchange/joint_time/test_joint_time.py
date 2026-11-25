import json

import numpy as np
from scipy.linalg import fractional_matrix_power

from model.executable_model import federation
from model.experiments.interval_echo_clock.modular_reference import reference_reader
from model.experiments.interval_echo_clock.reference_exchange.joint_time import measure


def records():
    return json.loads((measure.HERE / "records.json").read_bytes())


def test_joint_evolution_from_density_matrix():
    launches, settings, _ = reference_reader.instruments(federation(0))
    sharp = launches[0][1][:2]
    source = np.kron(np.outer(sharp, sharp.conj()), np.eye(2) / 2)
    bars = [(u.conj().T @ e @ u)[:2, :2] for e, u in settings]
    basis = np.kron(measure.X_BASIS, measure.X_BASIS)
    report = records()
    for name, weights in [
        ("shared_exchange", [5 / 12, 1 / 4, 1 / 4, 1 / 12]),
        ("independent_exchange", [4 / 9, 2 / 9, 2 / 9, 1 / 9]),
    ]:
        rho = basis @ np.diag(weights) @ basis.T
        for row in report["inference"][name]["held"]:
            theta = row["accepted"] * measure.ANGLE / np.log(2)
            # Spectral power on the full four-dimensional tensor product.
            lam, v = np.linalg.eigh(rho)
            u = (v * (lam ** (1j * theta))) @ v.conj().T
            pred = measure.outcomes(
                u @ source @ u.conj().T, bars[row["analyzer"]], bars[2]
            ).reshape(2, 2)
            np.testing.assert_array_less(
                abs(np.array(row["forecast"]) - pred), row["joint_radius"]
            )
        np.testing.assert_allclose(
            fractional_matrix_power(rho, 0.5) @ fractional_matrix_power(rho, 0.5),
            rho,
            atol=1e-12,
        )


def test_trial_accounting_and_counted_clock_control():
    raw = records()["records"]
    for data in raw["calibration"].values():
        assert sum(data) == raw["shots"]
    for mode, bins in raw["clock_counts"].items():
        data = np.array(bins)
        assert data.sum() == raw["shots"] and data.shape == (17, 4)
        foil = mode.startswith("contrast_only")
        quarter = mode.endswith("_1")
        for n, row in enumerate(data):
            total = row.sum()
            expected = (
                1 + np.cos((0 if foil else n * measure.ANGLE) + quarter * np.pi / 2)
            ) / 2
            assert abs(row[2:].sum() / total - expected) <= 2 * measure.margin(total)
            if foil and not quarter:
                assert row[:2].sum() == 0


def test_record_only_decisions_and_reference_choice():
    data = records()
    from model.candidate_laws.interval_echo.openwave.electron_field.replay import (
        compare,
    )

    decoded = measure.summarize(data["records"])
    compare(decoded, data["inference"], "clock decode")
    assert sum(r["disjoint_entries"] for r in decoded["shared_exchange"]["held"]) == 22
    assert (
        sum(r["disjoint_entries"] for r in decoded["independent_exchange"]["held"]) == 0
    )
    for k, v in decoded.items():
        row = next(r for r in v["held"] if r["accepted"] == 8 and r["analyzer"] == 0)
        separated = (
            abs(row["forecast_covariance"] - row["measured_covariance"])
            > row["forecast_covariance_error"] + row["measured_covariance_error"]
        )
        assert separated == (k == "shared_exchange")
