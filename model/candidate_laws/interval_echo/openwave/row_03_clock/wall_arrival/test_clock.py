"""Independent physical orientation calculation, and checks on held records."""

import json

import numpy as np
import pytest
from scipy.linalg import eigh

from .apparatus import Strand, word_law
from .measure import DELAYS, HERE, LAWS, SHOTS, read


def physical(n, j, h):
    states = np.arange(1 << n)
    spins = 1 - 2 * ((states[:, None] >> np.arange(n)) & 1)
    z = np.column_stack([np.ones(len(states)), spins, -np.ones(len(states))])
    signs = z[:, :-1] * z[:, 1:]
    mat = np.diag(j * signs.sum(axis=1))
    for i in range(n):
        mat[states, states ^ (1 << i)] -= h
    return signs, spins, eigh(mat)


def evolve(eig, state, t):
    es, vs = eig
    return vs @ (np.exp(-1j * t * es) * (vs.T @ state))


@pytest.mark.parametrize("j,h", [(8, 0.5), (4, 0.5), (0, 7.5)])
def test_independent_original_orientation_words(j, h):
    n = 7
    prep_j = 4 if j == 4 else 8
    signs, spins, prep = physical(n, prep_j, 0.5)
    _, _, flight = physical(n, j, h)
    rng = np.random.default_rng(332)
    f = rng.normal(size=n + 1) + 1j * rng.normal(size=n + 1)
    f /= np.linalg.norm(f)
    state = np.zeros(1 << n, dtype=complex)
    for address, amplitude in enumerate(f):
        idx = np.flatnonzero(((signs == 1).sum(axis=1) == 1) & (signs[:, address] == 1))
        assert len(idx) == 1
        state[idx[0]] = amplitude
    initial = evolve(prep, state, 2)
    times = [0, 0.031, 1.37, 4.5]
    expected = []
    for t in times:
        panel = []
        for outcome in (-1, 1):
            released = evolve(flight, initial * (spins[:, 0] == outcome), t)
            for a in (-1, 1):
                for b in (-1, 1):
                    panel.append(
                        np.sum(
                            abs(released[(spins[:, 0] == a) & (spins[:, 6] == b)]) ** 2
                        )
                    )
        expected.append(panel)
    g = Strand(n, prep_j).enroll(f)
    actual = word_law(g, Strand(n, j, h), times)
    np.testing.assert_allclose(actual, expected, atol=2e-12, rtol=0)


def test_archive():
    with np.load(HERE / "counts.npz") as raw:
        for n, name, *_ in LAWS:
            bank = raw[f"n{n}_{name}_words"]
            assert bank.shape == (len(DELAYS), 8)
            assert np.all(bank >= 0) and np.all(bank.sum(axis=1) == SHOTS)
            at_trigger = bank[0].reshape(2, 2, 2)
            assert at_trigger[0, 1].sum() == at_trigger[1, 0].sum() == 0
        summary = read(raw)
    saved = json.loads((HERE / "results.json").read_text())
    for observed, archived in zip(summary["recoil"], saved["recoil"], strict=True):
        np.testing.assert_allclose(
            list(observed.values()), list(archived.values()), atol=1e-11
        )
    for observed, archived in zip(
        summary["arrival_and_clock"], saved["arrival_and_clock"], strict=True
    ):
        for key in (
            "first_amplitude",
            "seventh_amplitude",
            "excess",
            "minimum_count_per_trigger",
        ):
            np.testing.assert_allclose(observed[key], archived[key], atol=1e-11)
        if observed["law"] == "native":
            assert (
                max(observed["first_amplitude"] + observed["seventh_amplitude"]) < 0.01
            )
            assert observed["excess"] > 5 * observed["excess_sampler_se"]
        if observed["law"] == "planted":
            assert min(observed["first_amplitude"]) > 0.97


def test_many_body_boundary_and_rational_lines():
    from fractions import Fraction

    from .eighth_port import gap
    from .rational_lines import PRIMES, certificate

    lo, hi = gap()
    assert (lo, hi) == (Fraction(1523386, 100000), Fraction(1523402, 100000))
    for m, p in PRIMES.items():
        prime, checks, full = certificate(m)
        assert prime == p and full and all(degree == 0 for _, degree in checks)


def test_enrolled_source_and_sampled_word_laws():
    from .apparatus import source

    with np.load(HERE / "counts.npz") as raw:
        for n, name, j, after_j, after_h in LAWS:
            prefix = f"n{n}_{name}_"
            near = source(n)
            assert near[0] == near[-1] == 0
            initial = Strand(n, j)
            enrolled = initial.enroll(near)
            # Fixed samples across the recorded duration, including the zero-delay bin.
            idx = np.array([0, 37, 114, 1011, 2047])
            probabilities = word_law(enrolled, Strand(n, after_j, after_h), DELAYS[idx])
            observed = raw[prefix + "words"][idx] / SHOTS
            allowed = (
                8 * np.sqrt(probabilities * (1 - probabilities) / SHOTS) + 2 / SHOTS
            )
            assert np.all(abs(observed - probabilities) < allowed)
            far = initial.enroll(source(n, mirror=True))
            p = (1 + np.array([enrolled[0, n + 1], far[0, n + 1]])) / 2
            assert np.all(
                abs(raw[prefix + "arrivals"] / SHOTS - p)
                < 8 * np.sqrt(p * (1 - p) / SHOTS)
            )
