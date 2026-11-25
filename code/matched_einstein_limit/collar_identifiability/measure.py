"""Positive aggregate of conditional signed-word readings on the native three-level collar."""

import json
from fractions import Fraction
from math import comb, log, sqrt

import numpy as np
from scipy.fft import irfft, next_fast_len, rfft
from scipy.linalg import expm
from scipy.stats import binom

from model.experiments.interval_echo_einstein_station import backreaction
from model.experiments.material_inertia.dark_collar import label

# A directed cyclic word measured by the two registered local word readers.
FORWARD = np.roll(np.eye(3, dtype=int), 1, axis=1)
SCORE_SCALE = log(3)


def word_channel(exchange_rate):
    """Keep the local preparation fixed. Apply the echo, degree phase and terminal exchange."""
    law = backreaction.law
    exchange, contrast = law.local_term_parts()
    term = exchange_rate * exchange + backreaction.KAPPA * contrast
    phase = np.tile([np.exp(-1j * np.pi * a**2 / 8) for a, _ in law.SRC], 1 + law.N_LAB)
    reader = expm(-1j * np.pi / 2 * term) @ (phase[:, None] * law.echo())
    word_indices = [law.SRC.index(pair) for pair in ((1, 1), (0, 2), (2, 0))]
    loads = np.eye(law.N_SRC)[:, word_indices[0]]
    records = [label(m, parity) for m, parity in ((0, 1), (2, 1), (2, -1))]
    initial = np.column_stack([np.kron(word, loads) for word in records])
    probabilities = abs((reader @ initial)[word_indices])**2
    # Any other endpoint outcome would have to be included in the instrument.
    np.testing.assert_allclose(probabilities.sum(axis=0), 1., rtol=0, atol=5e-14)
    return probabilities


def source_inputs():
    source = backreaction.correlated_contacts()["ordered_three_levels"]
    arms = np.array([arm["source_units_of_1_over_36"] for arm in source["arms"]])
    energies = np.rint(source["source_energies"]).astype(int)
    np.testing.assert_allclose(energies, source["source_energies"], rtol=0, atol=1e-12)
    assert np.array_equal(arms[0] + arms[2], 2 * arms[1])
    np.testing.assert_array_equal(arms, [4*np.eye(3, dtype=int) + (6-t)*FORWARD + (2+t)*FORWARD.T
                                         for t in (-1, 0, 1)])
    return arms, energies


def signed_word_gaps(arms, energies):
    """All total-energy bins, including degeneracies, and both one-register marginals."""
    pair_energy = (energies[:, None] + energies).ravel()
    histograms = np.array([np.bincount(pair_energy, weights=arm.ravel()) for arm in arms])
    np.testing.assert_array_equal(histograms, np.tile(histograms[1], (3, 1)))
    np.testing.assert_array_equal(arms.sum(axis=1), np.full((3, 3), 12))
    np.testing.assert_array_equal(arms.sum(axis=2), np.full((3, 3), 12))
    score = np.log(4 / arms[1])
    delta = (arms[2] - arms[0]) / 72
    np.testing.assert_allclose(abs(np.sum(delta * score)), SCORE_SCALE / 12, rtol=0, atol=1e-14)
    return {int(e): int(histograms[1, e]) for e in np.flatnonzero(histograms[1])}


def persistence_transfer(arms):
    """All two-write word laws share marginals on every proper subset of the four records."""
    # Integer numerator at denominator 36 for a one-word probability direction.
    direction = (arms[2] - arms[0]) // 2
    np.testing.assert_array_equal(arms[2] - arms[0], 2 * direction)
    two = np.multiply.outer(direction, direction)
    for leg in range(4):
        np.testing.assert_array_equal(two.sum(axis=leg), np.zeros((3, 3, 3), dtype=int))
    denominator = 36**2
    trace_distance = Fraction(int(abs(two).sum()), 2 * denominator)
    assert trace_distance == Fraction(1, 72)
    reference = np.multiply.outer(arms[1], arms[1])
    coincidence = FORWARD[:, :, None, None] == FORWARD[None, None, :, :]
    # K=2 and N=1: a difference in trials has zero score; equal trials have score (log 3)/2.
    folded_slope = SCORE_SCALE * np.sum(two * coincidence) / (2 * denominator)
    target = np.log(3) / 144
    np.testing.assert_allclose(folded_slope, target, rtol=0, atol=1e-14)
    effect = np.multiply.outer(FORWARD, FORWARD)
    forward_covariance = Fraction(int(np.sum(two * effect)), denominator)
    assert forward_covariance == Fraction(1, 144)
    slopes = {}
    for s in (-1., -.5, 0., .5, 1.):
        rho = (reference + s*two) / denominator
        assert np.all(rho > 0)
        measured = -np.sum((two / denominator) * np.log(rho))
        exact = -log((36 + s) * (4 + s) / (12 - s)**2) / 144
        np.testing.assert_allclose(measured, exact, rtol=0, atol=5e-16)
        slopes[str(s)] = float(measured)
    return dict(two_write_trace_distance=str(trace_distance), forward_covariance=str(forward_covariance),
                signed_modular_slopes=slopes, folded_amplitude_slope=float(folded_slope))



def click_laws(arms, channel):
    return np.array([np.sum(FORWARD * (channel @ (arm / 36) @ channel.T)) for arm in arms])


def mean_absolute_signs(n):
    """E abs(T_1 + ... + T_n) as a rational for independent, fair signed tags."""
    return Fraction(n * comb(n - 1, (n - 1) // 2), 2**(n - 1))


def finite_aggregate_mean(sources, repeats, p_minus, p_plus):
    """Average |sum_i(mean_j Y_ij - p0)| over all readouts, using finite convolutions."""
    k = np.arange(repeats + 1)
    one_source = (binom.pmf(k, repeats, p_minus) + binom.pmf(k, repeats, p_plus)) / 2
    trials = sources * repeats
    fft_size = next_fast_len(trials + 1)
    chances = irfft(rfft(one_source, fft_size)**sources, fft_size)[:trials + 1]
    if chances.min() < -1e-12 or abs(chances.sum() - 1) > 2e-11:
        raise ArithmeticError("The complete folded count distribution lost probability.")
    centered_counts = np.arange(trials + 1) - trials * (p_minus + p_plus) / 2
    return float(SCORE_SCALE * (chances @ abs(centered_counts)) / repeats)


def run():
    arms, energies = source_inputs()
    histogram = signed_word_gaps(arms, energies)
    native = click_laws(arms, word_channel(1))
    no_exchange = click_laws(arms, word_channel(0))
    np.testing.assert_allclose(native, [7/12, 1/2, 5/12], rtol=0, atol=5e-14)
    np.testing.assert_allclose(no_exchange, np.full(3, no_exchange[1]), rtol=0, atol=5e-14)
    rows = []
    for n in [1, 4, 16, 64]:
        signs = float(mean_absolute_signs(n))
        for k in [1, 64, 4096]:
            paired = finite_aggregate_mean(n, k, native[0], native[2])
            fresh = finite_aggregate_mean(n, k, native[1], native[1])
            muted = finite_aggregate_mean(n, k, no_exchange[0], no_exchange[2])
            floor = finite_aggregate_mean(n, k, no_exchange[1], no_exchange[1])
            np.testing.assert_allclose(muted, floor, atol=3e-11)
            target = SCORE_SCALE * signs / 12
            bound = SCORE_SCALE * sqrt(n * 35 / (144*k))
            assert -3e-11 <= paired - target <= bound + 3e-11
            rows.append(dict(sources=n, reads_per_source=k, paired=paired, fresh=fresh,
                             no_reversal_exchange=muted, no_reversal_fresh=floor, infinite_reads=target,
                             tag_relative_scatter=sqrt(n / signs**2 - 1)))
    return dict(energy_multiplicities_per_36=histogram, native_forward=native.tolist(),
                no_reversal_forward=no_exchange.tolist(), two_writes=persistence_transfer(arms), comparisons=rows)


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
