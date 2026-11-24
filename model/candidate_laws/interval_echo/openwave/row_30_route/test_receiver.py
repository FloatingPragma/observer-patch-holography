"""Check the complete temporal filter and a remote outcome against its leakage."""

import json

import numpy as np
from scipy.linalg import expm

from model.candidate_laws.interval_echo.openwave.row_30_route.receiver import HERE, PULSES


def test_coherent_positive_holds_compose_to_advertised_filter():
    h = np.array([[-4.1, 0.3], [0.3, -3.65]])
    detuning = h + 4 * np.eye(2)
    held = expm(-1j * detuning / 8)
    k_plus = (np.eye(2) + held) / 2
    k_minus = (np.eye(2) - held) / 2
    np.testing.assert_allclose(
        k_plus.conj().T @ k_plus + k_minus.conj().T @ k_minus, np.eye(2), atol=1e-14
    )
    energies, vecs = np.linalg.eigh(detuning)
    f = np.cos(energies / 16) ** PULSES * np.exp(-1j * PULSES * energies / 16)
    np.testing.assert_allclose(
        np.linalg.matrix_power(k_plus, PULSES), (vecs * f) @ vecs.T, atol=1e-12
    )


def test_all_launches_and_gap_transmission_witness():
    d = json.loads((HERE / "receiver.json").read_text())
    for r in d["records"]:
        assert sum(r["accepted_by_cell"]) + r["rejected"] == d["shots"]
    selected = {
        r["route"]: r
        for r in d["records"]
        if r["degree"]
        and r["source_index"] == 97
        and r["frequency"] == -4
        and r["extra_hold"] == 0
    }
    assert selected["shortest"]["accepted_by_cell"][101] == 2220
    assert sum(selected["four"]["accepted_by_cell"]) == 0
    assert np.cos(1 / 16) ** (2 * PULSES) < 1.25e-14
