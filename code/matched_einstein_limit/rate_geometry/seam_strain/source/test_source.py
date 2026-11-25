import json

import numpy as np
import pytest
from scipy import sparse

from . import certificate as cert
from . import measure as m


def archived(name):
    return json.loads((m.HERE / name).read_text())


def test_native_rail_energy_and_short_time_strain():
    for level in (0, 1):
        board = m.formation.federation(level)
        route = m.formation.region(board)
        start = m.launch(board, route)
        eye_rails = np.diag(np.repeat((16, 18), 4))
        for sign in (-1, 1):
            a, h2, b2 = m.system(board, route, sign)
            difference = h2 / 2 - m.formation.hamiltonian(board, route, sign, 8)
            scalar = difference.diagonal()[0]
            assert np.max(abs(difference - scalar * sparse.eye(h2.shape[0]))) < 1e-12
            n = h2.shape[0]
            # Real isometries are stored with factor two.  B acts solely on rail.
            np.testing.assert_array_equal(start[:n].T @ (b2 @ start[:n]), 2 * eye_rails)
            diff = np.diag(np.repeat((-1, 1), 4)) / 12
            for op in (h2, b2):
                assert np.trace(start[:n].T @ (op @ start[:n]) @ diff) == pytest.approx(
                    4 / 3
                )
            # i[H,Y] from the Schrödinger derivative, without a time-difference fit.
            v = start / 2
            dv = -1j * (a @ v) / 40
            cross = dv[n:].conj().T @ v[:n] + v[n:].conj().T @ dv[:n]
            yprime = 1j * (cross - cross.conj().T)
            np.testing.assert_allclose(
                yprime, -2 * (9 / 40) * (eye_rails / 2), atol=4e-15, rtol=0
            )
            assert np.trace(yprime @ diff).real == pytest.approx(-3 / 20)


def test_replay_source_charge_controls_and_held_predictions():
    regenerated = m.run()
    saved = archived("results.json")
    assert len(regenerated["runs"]) == len(saved["runs"]) == 32
    for data, ref in zip(regenerated["formations"], saved["formations"]):
        for state in ("formed", "decorrelated"):
            p = np.array(data["head_pauli_probabilities"][state])
            assert p.min() >= 0 and p.shape == (3, 3, 2, 2)
            np.testing.assert_allclose(p.sum(axis=(-2, -1)), 1, atol=1e-14, rtol=0)
            np.testing.assert_allclose(
                p, ref["head_pauli_probabilities"][state], atol=2e-10, rtol=0
            )
        p, q = (
            np.array(data["head_pauli_probabilities"][t])
            for t in ("formed", "decorrelated")
        )
        for axis in (-1, -2):
            np.testing.assert_allclose(
                p.sum(axis=axis), q.sum(axis=axis), atol=1e-14, rtol=0
            )
        assert np.max(abs(p - q)) > 1e-6
    assert len(regenerated["zero_epsilon_controls"]) == 16
    assert (
        max(
            r["largest_minus_x_projection_norm"]
            for r in regenerated["zero_epsilon_controls"]
        )
        < 2e-11
    )
    for data, ref in zip(regenerated["runs"], saved["runs"]):
        assert all(
            data[k] == ref[k]
            for k in (
                "level",
                "forming_charge",
                "pump_charge",
                "time",
                "cells",
                "route",
                "epsilon",
                "omega",
            )
        )
        arms = data["arms"]
        for name, arm in arms.items():
            stored = ref["arms"][name]
            for k in (
                "bare_energy",
                "whole_bridge_energy",
                "coupled_energy",
                "bloch",
                "density_real",
                "density_imag",
                "modular_nats",
            ):
                np.testing.assert_allclose(arm[k], stored[k], atol=5e-9, rtol=0)
            r = np.array(arm["density_real"]) + 1j * np.array(arm["density_imag"])
            assert np.min(np.linalg.eigvalsh(r)) >= -1e-13
            reconstructed = (
                np.eye(2) + sum(v * p for v, p in zip(arm["bloch"], m.formation.PAULI))
            ) / 2
            np.testing.assert_allclose(r, reconstructed, atol=5e-12, rtol=0)

        def delta(name, k):
            return arms[name][k] - arms["reference"][k]

        for energy in ("bare_energy", "whole_bridge_energy", "coupled_energy"):
            assert abs(delta("decorrelated", energy)) < 2e-12
        assert delta("ordinary", "bare_energy") == pytest.approx(1 / 3, abs=2e-12)
        assert delta("ordinary", "coupled_energy") == pytest.approx(41 / 120, abs=2e-12)
        assert delta("ordinary", "modular_nats") == pytest.approx(
            np.log(2) / 3, abs=2e-12
        )
        assert delta("decorrelated", "modular_nats") > 0.0019
        for k, t in data["transfer"].items():
            np.testing.assert_allclose(t, ref["transfer"][k], atol=3e-8, rtol=0)


def test_rational_certificate_recovers_archived_moments_and_noncollinearity():
    data = cert.certificate(archived("results.json"))
    saved = archived("certificate.json")
    assert len(data["checks"]) == len(saved["checks"]) == 32
    for row, ref in zip(data["checks"], saved["checks"]):
        assert all(
            row[k] == ref[k] for k in ("level", "forming_charge", "pump_charge", "time")
        )
        # The ordinary-line bound uses the exact Gram center. The modular
        # certificate uses the current platform's checked eigenspace proposal.
        line = "collar_distance_from_ordinary_line_lower"
        assert row[line] == ref[line] and row[line] > 1.99e-6
        for check in (row, ref):
            assert check["modular_residual_max_component_lower"] > 1e-4
        assert all(e < 1e-9 for e in row["moment_absolute_errors"].values())
        assert all(e < 1e-9 for e in ref["moment_absolute_errors"].values())
        assert (
            row["nominal_column_error"] < 1e-9
            and row["fixed_point_column_error"] < 2.2e-10
        )
    for state, ref in zip(data["source"], saved["source"]):
        assert (state["cells"], state["forming_charge"]) == (
            ref["cells"],
            ref["forming_charge"],
        )
        for record in (state, ref):
            a, b = record["collar_modular_nats"]
            assert 0 < a <= b and b - a < 4e-5
            assert record["nominal_formed_trace_norm_error"] < 2e-12
        a, b = state["collar_modular_nats"]
        x, y = ref["collar_modular_nats"]
        assert max(a, x) <= min(b, y)
