import json

import numpy as np
import pytest

from model.benchmarks.eight_slot_torus.electron_wall import phase_clock as pc
from model.benchmarks.eight_slot_torus.electron_wall import walk as w
from model.benchmarks.eight_slot_torus.electron_wall.cyclotron import magnetic_block, oscillator_sources


def test_fourier_counter_against_coherent_native_services():
    n, length, ticks, grade, q = 5, 2, 32, -1, 1
    conn = w.links((n, 1, 1))
    conn[1, :, 0, 0] = np.exp(2j * np.pi * np.arange(n) / n)
    source = next(oscillator_sources(n, length, grade))
    a = source.reshape(n, length, 2) * [1, 1j] / np.sqrt(2)
    b = grade * a[:, ::-1] * [1, -1]
    state = w.prepare(np.stack([a, b], axis=1)[:, None, None])
    history = []
    for _ in range(ticks):
        history.append(state.ravel())
        state = w.step(state, conn, q)
    counter = np.fft.fft(history, axis=0) / ticks
    read = (abs(counter) ** 2).sum(axis=1)
    roots, vectors = np.linalg.eigh(magnetic_block(n, length, q, grade).toarray())
    forecast = np.zeros(ticks)
    for root, weight in zip(roots, abs(vectors.T @ source) ** 2):
        phase = np.arccos(root / 12.)
        for sign in (-1, 1):
            forecast += weight * pc.fourier_probabilities(sign * phase, ticks) / 2
    np.testing.assert_allclose(read, forecast, atol=1e-13)


def test_complete_counts_and_finite_resource_gates():
    receipt = json.loads((pc.HERE / "clock_counts.json").read_text())
    lines = (json.loads((pc.HERE / "cyclotron_lines.json").read_text())
             + json.loads((pc.HERE / "weak_lines.json").read_text()))
    result = pc.verify(receipt, pc.SPEC, lines)
    assert result["log10_family_error_bound"] < -8
    assert result["series"][0]["sampled"]
    assert not result["series"][-1]["sampled"]
    for bad in ([[0, 1], [0, 2]], [[0, 1.5], [1, 1.5]], [[-3, 3]]):
        with pytest.raises(ValueError):
            pc.median_rate(bad, 8, 3)
    with pytest.raises(ValueError):
        pc.verify({}, pc.SPEC)
    for missing in ([], lines[:-1], lines[2:], lines[:2] + lines):
        with pytest.raises(ValueError, match="incomplete magnetic clock series"):
            pc.verify(receipt, pc.SPEC, missing)


def test_bin_gate_rejects_counts_in_an_impossible_frequency_port():
    # One eigenphase at pi/2 deposits every trial in folded bin zero.
    phases, weights = np.array([np.pi / 2]), np.array([1.])
    pc.check_bins([[0, 2049]], 32, 2049, phases, weights)
    with pytest.raises(ArithmeticError, match="retained bin"):
        pc.check_bins([[1, 2049]], 32, 2049, phases, weights)
