import json
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse.linalg import expm_multiply

from . import hidden, history, measure, readout

HERE = Path(__file__).parent


@pytest.mark.parametrize("level", [0, 1])
def test_uncoupled_material_record_at_full_band_exposures(level):
    panel = readout.panel(level, epsilon=0)
    gaps = np.array(panel["slow_gaps"])
    weights = np.array(panel["spectral_probabilities"])
    overlap = panel["rates"]["position_score_overlap"]
    for row in panel["appointments"]:
        value = weights @ np.sin(gaps * row["time"])
        derivative = weights @ (gaps * np.cos(gaps * row["time"]))
        for name, factor in (("band_clock", 1), ("band_position", overlap)):
            np.testing.assert_allclose(row[name], [factor * value, 0, 0, 0], atol=1e-11)
            np.testing.assert_allclose(
                row[name + "_derivative"], [factor * derivative, 0, 0, 0], atol=1e-12
            )


def test_release_against_direct_port_probabilities_and_local_clock():
    h, b, g, u, f, gap, _, _, edge = readout.modes(0)
    rates = readout.initial_reads(b, g, u, f, gap, edge)
    gen = -1j * readout.reciprocal(h, b, readout.EPSILON, readout.OMEGA)
    dt = 1e-5
    for index, multiplier in enumerate([0.8, 1.25]):
        qubit = np.eye(2)[:, index]
        launch = np.column_stack(
            [np.kron(qubit, (g + phase * u) / np.sqrt(2)) for phase in (1j, -1j, 1)]
        )
        records = []
        for sign in (-1, 1):
            t = sign * dt
            evolved = expm_multiply(t * gen, launch, traceA=t * gen.diagonal().sum())
            joint = evolved.reshape(2, len(g), 3)
            terminal = (abs(joint) ** 2).sum(axis=0)
            np.testing.assert_allclose(terminal.sum(axis=0), 1, atol=1e-13)
            position = f @ (terminal[:, 0] - terminal[:, 1]) / 2
            # Measure the two resolved analyzer populations; all complement
            # outcomes score zero, including amplitude scattered out of the band.
            click = [
                sum(abs(joint[:, :, 2] @ (g + k * 1j * u) / np.sqrt(2)) ** 2)
                for k in (1, -1)
            ]
            records.append([position, click[0] - click[1]])
        slope = (np.array(records[1]) - records[0]) / (2 * dt)
        expected = [
            rates["position_rate"]
            * (1 + (multiplier - 1) * rates["fraction_of_position_rate_on_bridge"]),
            gap * (1 + (multiplier - 1) * rates["fractional_clock_response_to_bridge"]),
        ]
        np.testing.assert_allclose(slope, expected, atol=2e-10)
        for lo, hi in ((0, 8), (8, 9)):
            local = readout.local_effect(lo, hi, dt, readout.EPSILON, readout.OMEGA)
            np.testing.assert_allclose(
                qubit @ local @ qubit / dt, (hi - lo) * multiplier, atol=3e-8
            )


def test_signed_source_transfer_and_receipt():
    stored = json.loads((HERE / "results.json").read_text())
    fresh = measure.run()
    assert len(fresh["source_comparisons"]) == len(stored["source_comparisons"]) == 32
    for old, new in zip(stored["readout_panels"], fresh["readout_panels"], strict=True):
        for a, b in zip(old["appointments"], new["appointments"], strict=True):
            assert set(a) == set(b)
            for key in a:
                np.testing.assert_allclose(a[key], b[key], atol=1e-10)
    for a, b in zip(
        stored["source_comparisons"], fresh["source_comparisons"], strict=True
    ):
        for key in (
            "cells",
            "forming_charge",
            "pump_charge",
            "pump_time",
            "collar_to_ordinary_modular_ratio",
        ):
            assert a[key] == pytest.approx(b[key], abs=1e-11)
        np.testing.assert_allclose(
            list(a["initial_bridge_multiplier"].values()),
            list(b["initial_bridge_multiplier"].values()),
            atol=1e-12,
        )
        for x, y in zip(a["band_crossings"], b["band_crossings"], strict=True):
            for name in ("band_clock", "band_position"):
                for key in x[name]:
                    assert x[name][key] == pytest.approx(
                        y[name][key], rel=1e-7, abs=2e-10
                    )
    selected = next(
        r
        for r in fresh["source_comparisons"]
        if (r["cells"], r["forming_charge"], r["pump_charge"], r["pump_time"])
        == (80, -1, -1, 4)
    )
    for row in selected["band_crossings"]:
        for key in ("band_clock", "band_position"):
            assert row[key]["collar"] < 0 < row[key]["prediction"]


def test_rank_certificate_on_the_shared_formed_reference():
    stored = json.loads((HERE / "hidden.json").read_text())
    checked = hidden.certificate()
    assert checked["rank_with_modular_read"] == 4
    assert checked["rank_qubit"] == 3
    assert checked["inverse_error_upper"] < 1e-3
    assert checked["selected_columns"] == stored["selected_columns"]
    np.testing.assert_allclose(
        checked["minor_midpoint"], stored["minor_midpoint"], atol=1e-10
    )


def test_history_recovers_the_fixed_modular_functional_from_counts():
    computed = history.run()
    stored = json.loads((HERE / "history.json").read_text())
    assert computed["rank"] == 9 and computed["inverse_error_upper"] < 0.01
    assert computed["settings"] == stored["settings"]
    assert computed["uniform_decoder_error_upper"] < 0.0001
    weights = np.array(computed["modular_decoder_weights"])
    np.testing.assert_allclose(weights, stored["modular_decoder_weights"], atol=2e-6)
    for old, fresh in zip(stored["reads"], computed["reads"], strict=True):
        n = old["all_trials_per_arm"]
        counts = np.array(old["plus_counts"])
        assert counts.shape == (24, 2) and np.all((0 <= counts) & (counts <= n))
        decoded = float(weights @ (2 * (counts[:, 1] - counts[:, 0]) / n))
        assert decoded == pytest.approx(old["decoded_modular_nats"], abs=1e-7)
        assert old["exact_population_prediction"] == pytest.approx(
            fresh["exact_population_prediction"], abs=1e-8
        )
        assert abs(decoded - fresh["direct_modular_nats"]) <= (
            computed["uniform_decoder_error_upper"]
            + computed["simultaneous_count_error_99pct"]
        )


def test_ordinary_calibration_with_independently_counted_collar_charge():
    tomography = json.loads((HERE / "history.json").read_text())
    charge = tomography["reads"][0]["decoded_modular_nats"]
    width = (
        tomography["uniform_decoder_error_upper"]
        + tomography["simultaneous_count_error_99pct"]
    )
    run = next(
        r
        for r in json.loads((HERE / "results.json").read_text())["source_comparisons"]
        if (r["cells"], r["forming_charge"], r["pump_charge"], r["pump_time"])
        == (80, -1, -1, 4)
    )
    for name in ("band_position", "band_clock"):
        row = run["band_crossings"][0][name]
        lower, upper = (
            row["ordinary"] * (charge + np.array([-width, width])) * 3 / np.log(2)
        )
        assert row["collar"] < 0 < lower < row["prediction"] < upper


def test_held_signs_by_independent_dense_population_flight():
    from scipy.linalg import eigh

    h, bond, g, u, score, gap, _, _, _ = readout.modes(1)
    values, vectors = eigh(
        readout.reciprocal(h, bond, readout.EPSILON, readout.OMEGA).toarray(),
        driver="evd",
    )
    run = next(
        r
        for r in json.loads((measure.sources.HERE / "results.json").read_text())["runs"]
        if (r["level"], r["forming_charge"], r["pump_charge"], r["time"])
        == (1, -1, -1, 4)
    )
    launch = np.column_stack([(g + phase * u) / np.sqrt(2) for phase in (1j, -1j, 1)])
    coords = vectors.T @ np.kron(np.eye(2), launch)
    t = np.pi / (6 * gap)
    dt = 1e-4
    records = {}
    for offset in (-dt, 0, dt):
        wave = (
            vectors @ (np.exp(-1j * (t + offset) * values)[:, None] * coords)
        ).reshape(2, len(g), 2, 3)
        for name, arm in run["arms"].items():
            rho = np.array(arm["density_real"]) + 1j * np.array(arm["density_imag"])
            # Retain both launch signs and the full terminal probability distribution.
            distribution = np.einsum(
                "spia,ij,spja->pa",
                wave[:, :, :, :2],
                rho,
                wave[:, :, :, :2].conj(),
                optimize=True,
            ).real
            np.testing.assert_allclose(distribution.sum(axis=0), 1, atol=2e-12)
            position = score @ (distribution[:, 0] - distribution[:, 1]) / 2
            counts = []
            for sign in (1, -1):
                amplitude = np.einsum(
                    "p,spi->si", (g + sign * 1j * u) / np.sqrt(2), wave[:, :, :, 2]
                )
                counts.append(
                    np.einsum("si,ij,sj->", amplitude, rho, amplitude.conj()).real
                )
            records[offset, name] = np.array([position, counts[0] - counts[1]])
    slope = (records[dt, "reference"] - records[-dt, "reference"]) / (2 * dt)
    computed = (records[0, "decorrelated"] - records[0, "reference"]) / (t * slope)
    forecast = (
        run["transfer"]["modular_ratio"]
        * (records[0, "ordinary"] - records[0, "reference"])
        / (t * slope)
    )
    np.testing.assert_allclose(computed, [-0.000353929980, -0.000089370005], rtol=2e-6)
    np.testing.assert_allclose(forecast, [0.006361454411, 0.001597977687], rtol=2e-6)
