"""Independent small-domain instrument checks; no acquisition of the cell menu."""

import numpy as np
import pytest
from scipy.linalg import expm

from ...contrast_memory_spherical_v1 import law as native
from .coherent_radar import first_hit, labels, pair_records
from .law import collar_states


def _source(contact=.7):
    h0 = np.array([[.2, .9, .1], [.9, -.3, .6], [.1, .6, .5]])
    h1 = h0 + np.diag([0, contact, 0])
    return h0, h1


def _direct(h0, h1, sender, receiver, ticks, dt, kappa):
    """Explicit two-sector state, projected branches and a separate capture mode."""
    d = len(h0)
    # The final mode is an isolated scalar capture, with no contact potential.
    hs = []
    for r, h in enumerate((h0, h1)):
        expanded = np.zeros((d + 1, d + 1))
        expanded[:d, :d] = h
        hs.append(expanded + r * kappa * np.eye(d + 1))
    h = np.block([[hs[0], np.zeros_like(hs[0])], [np.zeros_like(hs[0]), hs[1]]])
    u = expm(-1j * dt * h)
    state = np.zeros(2 * (d + 1), dtype=complex)
    state[sender] = state[d + 1 + sender] = 1 / np.sqrt(2)
    records = {}

    def project(value, station):
        mask = np.zeros(len(value), dtype=bool)
        mask[[station, d + 1 + station]] = True
        return value * mask, value * ~mask

    for j in range(1, ticks + 1):
        outgoing, state = project(u @ state, receiver)
        return_state = outgoing.copy()
        for k in range(1, ticks + 1):
            returning, return_state = project(u @ return_state, sender)
            captured = np.zeros_like(returning)
            captured[d] = returning[sender]
            captured[2 * d + 1] = returning[d + 1 + sender]
            for _ in range(2 * ticks - j - k):
                captured = u @ captured
            records[j, k] = captured
        for _ in range(ticks - j):
            return_state = u @ return_state
        records[j, 0] = return_state
    for _ in range(ticks):
        state = u @ state
    records[0, 0] = state

    result = np.zeros((2, len(records), 2))
    for index, label in enumerate(labels(ticks)):
        arms = records[tuple(label)].reshape(2, d + 1)
        for axis, phase in enumerate((1, -1j)):
            result[axis, index] = [np.linalg.norm(arms[0] + sign * phase * arms[1]) ** 2 / 2
                                   for sign in (1, -1)]
    return result


@pytest.mark.parametrize("ticks", [1, 2, 4])
def test_complete_history_and_deadline_against_direct_projection(ticks):
    h0, h1 = _source()
    dt, kappa = .37, .83
    u = [expm(-1j * dt * h) for h in (h0, h1)]
    actual = pair_records(u, 0, 2, ticks, dt, kappa)
    expected = _direct(h0, h1, 0, 2, ticks, dt, kappa)
    np.testing.assert_allclose(actual["coherent"], expected, atol=2e-14, rtol=0)
    for name in ("coherent", "dephased", "independent_reprepared"):
        assert actual[name].min() >= -1e-14
        np.testing.assert_allclose(actual[name].sum(axis=(1, 2)), 1, atol=2e-14, rtol=0)


def test_first_hit_amplitudes_against_quantum_renewal():
    u = expm(-.41j * _source()[0])
    ticks = 7
    for a, b in ((0, 1), (1, 0), (0, 2)):
        direct, _ = first_hit(u, a, b, ticks)
        powers = [np.linalg.matrix_power(u, n) for n in range(ticks + 1)]
        renewal = []
        for n in range(1, ticks + 1):
            renewal.append(powers[n][b, a] - sum(powers[n - m][b, b] * renewal[m - 1]
                                                for m in range(1, n)))
        np.testing.assert_allclose(direct, renewal, atol=2e-14, rtol=0)


def test_imaginary_quadrature_sign_and_contact_free_factorization():
    h0, _ = _source(0)
    dt, ticks, kappa = .3, 3, .4
    u = expm(-1j * dt * h0)
    actual = pair_records([u, u], 0, 2, ticks, dt, kappa)
    history = actual["coherent"][0].sum(axis=-1)
    signed = actual["coherent"][..., 0] - actual["coherent"][..., 1]
    angle = 2 * ticks * dt * kappa
    np.testing.assert_allclose(signed[0], history * np.cos(angle), atol=2e-14, rtol=0)
    np.testing.assert_allclose(signed[1], -history * np.sin(angle), atol=2e-14, rtol=0)
    np.testing.assert_allclose(history, actual["independent_reprepared"][0].sum(axis=-1),
                               atol=2e-14, rtol=0)


def test_retained_memory_differs_from_independent_repreparation():
    h0, h1 = _source(2)
    actual = pair_records([expm(-.8j * h) for h in (h0, h1)], 0, 2, 4, .8, .7)
    history = actual["coherent"].sum(axis=-1)
    np.testing.assert_allclose(history, actual["dephased"].sum(axis=-1), atol=2e-14, rtol=0)
    assert np.max(abs(history - actual["independent_reprepared"].sum(axis=-1))) > 1e-4
    assert np.max(abs(actual["coherent"] - actual["dephased"])) > 1e-3
    assert np.max(abs(actual["coherent"] - actual["dephased"]).sum(axis=(1, 2)) / 2) <= .5 + 1e-14


def test_fixed_deadline_contact_perturbation_bound():
    dt, ticks, kappa = .31, 3, .9
    h0, h1 = _source(.4)
    baseline = pair_records([expm(-1j * dt * h) for h in (h0, h1)], 0, 2, ticks, dt, kappa)["coherent"]
    for change in (-.02, .03, .1):
        perturbed = pair_records([expm(-1j * dt * h0), expm(-1j * dt * (h1 + np.diag([0, change, 0])))],
                                 0, 2, ticks, dt, kappa)["coherent"]
        tv = abs(perturbed - baseline).sum(axis=(1, 2)) / 2
        assert np.max(tv) <= min(1, 2 * ticks * dt * abs(change)) + 2e-14


def test_native_two_arm_source_and_written_clock_reduction():
    exchange, square = native.local_term_parts()
    arms = collar_states(1)
    np.testing.assert_allclose(square @ arms, arms @ np.diag([0, 1]), atol=1e-13, rtol=0)
    np.testing.assert_allclose((exchange + 8 * square) @ arms, arms @ np.diag([0, 8]),
                               atol=1e-13, rtol=0)
    clock = np.column_stack([collar_states(sign)[:, 1] for sign in (1, -1)])
    np.testing.assert_allclose(square @ clock, clock, atol=1e-13, rtol=0)
    np.testing.assert_allclose((exchange + 8 * square) @ clock, clock @ np.diag([8, 9]),
                               atol=1e-13, rtol=0)
