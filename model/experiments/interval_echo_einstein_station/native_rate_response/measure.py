"""Test a fractional local rate with two load preparations and one terminal port read."""

from itertools import product

import numpy as np
from scipy.linalg import block_diag
from scipy.sparse.linalg import expm_multiply

from model.experiments.interval_echo_einstein_station import backreaction as seam


def integer_native_parts():
    """Return eight times each native term in the complete memory-load basis."""
    law = seam.law
    write = law._write_pattern()
    loads = np.eye(law.N_SRC, dtype=int)
    twice_projection = 2 * np.eye(write.shape[0], dtype=int) - write @ write.T
    reflection, contrast = law._label_ops()
    parts = []
    for twice_record in (np.eye(law.N_LAB, dtype=int) - reflection, 2 * contrast):
        generator = np.kron(twice_record, loads)
        mixing = write.T @ generator @ twice_projection
        np.testing.assert_array_equal(mixing, np.zeros_like(mixing))
        blank = 2 * write.T @ generator @ write
        written = twice_projection @ generator @ twice_projection
        parts.append(block_diag(blank, written))
    return parts


def double_tag_blocks():
    """Certify the invariant two-position block and its actions on all 11 memories."""
    law = seam.law
    exchange, contrast = integer_native_parts()
    native = exchange + int(seam.KAPPA) * contrast
    positions = [law.SRC.index(pair) for pair in ((2, 0), (0, 2))]
    columns = np.arange(1 + law.N_LAB)[:, None] * law.N_SRC + positions
    outside = np.setdiff1d(np.arange(law.DIM), columns.ravel())
    np.testing.assert_array_equal(native[np.ix_(outside, columns.ravel())], 0)
    left, right = columns.T
    blank = np.zeros((1 + law.N_LAB, 1 + law.N_LAB), dtype=int)
    blank[0, 0] = 1
    np.testing.assert_array_equal(native[np.ix_(right, left)], -4 * blank)
    np.testing.assert_array_equal(
        native[np.ix_(left, left)], native[np.ix_(right, right)]
    )
    return native, columns


def port_probability(generator, memory, load, duration, columns):
    state = np.zeros(generator.shape[0], dtype=complex)
    state[columns.ravel()] = np.kron(memory, load)
    final = expm_multiply(-1j * duration * generator, state)
    return float(np.sum(np.abs(final[columns[:, 1]]) ** 2))


def rate_experiment():
    native, columns = double_tag_blocks()
    native = native / 8
    sharp = np.array([1, 0])
    quadrature = np.array([1, 1j]) / np.sqrt(2)
    bare = np.eye(columns.shape[0])[0]
    # A coherent, nonstationary written memory spanning all ten labels.
    written = np.r_[
        0, np.exp(0.37j * np.arange(seam.law.N_LAB)) / np.sqrt(seam.law.N_LAB)
    ]
    rows = []
    for p, duration in product((0, 0.25, 0.5, 0.75, 1), (0.25, 0.75, np.pi / 2)):
        memory = np.sqrt(p) * bare + np.sqrt(1 - p) * written
        native_probs = [
            port_probability(native, memory, load, duration, columns)
            for load in (sharp, quadrature)
        ]
        # Same local Hilbert space, endpoint symmetry and route; both seam terms slowed.
        slowed_probs = [
            port_probability(p * native, bare, load, duration, columns)
            for load in (sharp, quadrature)
        ]
        np.testing.assert_allclose(
            native_probs,
            [p * np.sin(duration / 2) ** 2, (1 + p * np.sin(duration)) / 2],
            atol=3e-12,
            rtol=0,
        )
        np.testing.assert_allclose(
            slowed_probs,
            [np.sin(p * duration / 2) ** 2, (1 + np.sin(p * duration)) / 2],
            atol=3e-12,
            rtol=0,
        )
        rows.append((p, duration, native_probs, slowed_probs))
    return rows


if __name__ == "__main__":
    rows = rate_experiment()
    for p, duration, native, slowed in rows:
        if p == 0.5 and duration == np.pi / 2:
            print("Native coherent preparation", native)
            print("Slowed whole seam", slowed)
