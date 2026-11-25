"""Independent full-word evolution, timed-coherence identities and record replay."""

import itertools
import json

import numpy as np

from model.candidate_laws.interval_echo.openwave.row_23_coupling.timed_scattering import (
    measure as m,
)


def full_step(words, size, state, strength):
    """Apply each disjoint matching as its full matrix on unordered port words."""
    for p, theta in enumerate(m.pair.THETA):
        gate = m.native_matrix(theta, strength)
        pairs = [
            tuple(sorted(((x + p) % size, (x + p + 1) % size)))
            for x in range(0, size, 2)
        ]
        ops = np.ones((len(words), len(words)), complex)
        for a, b in pairs:
            pattern = [(int(a in word), int(b in word)) for word in words]
            local = np.eye(4, dtype=complex)
            local[np.ix_([1, 2], [1, 2])] = gate
            index = [2 * u + v for u, v in pattern]
            ops *= local[np.ix_(index, index)]
        state = ops @ state
    return state


def test_full_port_words_and_complex_timing():
    n = 5
    size = 2 * n
    rng = np.random.default_rng(2304)
    words = list(itertools.combinations(range(size), 2))
    r, b = map(np.array, zip(*(m.pair.canonical(n, *word) for word in words)))
    v = rng.normal(size=(size - 1, 2)) + 1j * rng.normal(size=(size - 1, 2))
    w = 2

    def lift(state, k):
        return np.exp(2j * np.pi * k * b / n)[:, None] * state[r] / np.sqrt(n)

    v /= np.linalg.norm(v)
    amp = m.native_matrix(m.static.DWELL, 8)[0, 0]
    for strength in (8, 0):
        spec = m.pair_spectrum(n, w, v, m.pair.transfers(n, strength), amp)
        U = full_step(words, size, np.eye(len(words), dtype=complex), strength)
        np.testing.assert_allclose(U.conj().T @ U, np.eye(len(words)), atol=2e-13)
        for j in range(size):
            click = (amp - 1) * np.exp(-2j * np.pi * j * np.array(words) / size).sum(
                axis=1
            ) / size + (j == 0)
            for t in m.DELAYS:
                u = np.linalg.matrix_power(U, t)
                early = u @ (click[:, None] * lift(v, w))
                late = click[:, None] * (u @ lift(v, w))
                np.testing.assert_allclose(
                    spec[t][:, j], m.norms(early, late), atol=1e-13
                )


def test_one_period_and_control_invariance():
    """All contrast in one full period reduces to a final diagonal word phase."""
    size = 10
    words = list(itertools.combinations(range(size), 2))
    eye = np.eye(len(words), dtype=complex)
    native, switched = (full_step(words, size, eye, k) for k in (8, 0))
    ratio = native @ switched.conj().T
    np.testing.assert_allclose(ratio, np.diag(np.diag(ratio)), atol=5e-13)
    np.testing.assert_allclose(abs(np.diag(ratio)), 1, atol=6e-13)
    collision_words = {tuple(sorted((x, (x + 1) % size))) for x in range(1, size, 2)}
    separate_phase = np.exp(-8j * 8 * m.pair.THETA[1])
    expected = np.array(
        [1.0 if word in collision_words else separate_phase for word in words]
    )
    np.testing.assert_allclose(np.diag(ratio), expected, rtol=0, atol=5e-13)
    amp = m.native_matrix(m.static.DWELL, 8)[0, 0]
    native, switched = (m.single_spectrum(size, k, amp) for k in (8, 0))
    for t in m.DELAYS:
        np.testing.assert_allclose(native[t], switched[t], atol=5e-13)


def test_archived_receipts_and_controls():
    results = json.loads((m.HERE / "results.json").read_text())
    size = results["ports"]
    with np.load(m.HERE / "counts.npz") as raw:
        assert len(raw.files) == 3 * 2 * 4 * 4
        for key, row in results["rows"].items():
            for t in m.DELAYS:
                rec = {}
                for k in (8, 0):
                    banks = [
                        raw[f"{key}_{k}_{t}_{arm}"]
                        for arm in ("early", "late", "real", "imag")
                    ]
                    for bank, num in zip(
                        banks, [size, size, 2 * size, 2 * size], strict=True
                    ):
                        assert (
                            bank.size == num + 2
                            and bank.min() >= 0
                            and bank.sum() == 2**40
                        )
                    rec[k] = m.read_banks(banks, size)
                reading = m.compare(rec[8], rec[0], row["radius_interval"])
                saved = row["delays"][str(t)]
                for metric, value in reading.items():
                    if value is None:
                        assert saved[metric] is None
                    else:
                        np.testing.assert_allclose(
                            value, saved[metric], rtol=5e-10, atol=5e-12
                        )
                if key == "single" or t < 2:
                    assert reading["detected_bins"] == 0 and reading["max_upper"] < 0.03
                else:
                    assert reading["first_above"] is not None
