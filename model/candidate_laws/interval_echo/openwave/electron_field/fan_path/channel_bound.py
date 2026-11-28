"""Certify two-channel admission, terminal decoding and separated-ring bounds."""

from fractions import Fraction
from math import isqrt

import numpy as np


def sqrt_bounds(value, scale=1 << 48):
    """Outward dyadic bounds for sqrt of a nonnegative rational."""
    value = Fraction(value)
    root = isqrt(value.numerator * scale**2 // value.denominator)
    return Fraction(root, scale), Fraction(root + 1, scale)


def dyadic_vectors(h, energies, columns, lower_cut, *, bits=40):
    """Enclose two distinct eigenvalues and both normalized real eigenvectors."""
    multiply = (2 * h).astype(np.int64)
    assert (multiply != 2 * h).nnz == 0
    scale = 2**bits
    trial = np.rint(columns * scale).astype(np.int64)
    assert np.max(abs(multiply).sum(axis=1)) * np.max(abs(trial)) < 2**60
    action = multiply @ trial
    norms = []
    centers, radii = [], []
    for i, e in enumerate(energies):
        nu = Fraction(round(float(e) * scale), scale)
        n2 = sum(int(z) ** 2 for z in trial[:, i])
        numerator = sum(
            (int(y) * nu.denominator - 2 * nu.numerator * int(z)) ** 2
            for y, z in zip(action[:, i], trial[:, i], strict=True)
        )
        res = Fraction(numerator, 4 * nu.denominator**2 * n2)
        norms.append(n2)
        centers.append(nu)
        radii.append(sqrt_bounds(res)[1])
    assert (
        centers[0] + radii[0]
        < centers[1] - radii[1]
        < centers[1] + radii[1]
        < lower_cut
    )
    error = Fraction()
    for i in range(2):
        gap = min(
            abs(centers[i] - centers[1 - i]) - radii[1 - i], lower_cut - centers[i]
        )
        assert gap > 0
        error += 2 * radii[i] / gap
    return trial, norms, centers, radii, error


def capture_floor(path, cell, trial, norms, frame_error):
    instrument = path.source(cell)
    sharp = np.sign(instrument).astype(np.int64)
    block_size = path.blocks[cell, 0].rows.stop - path.blocks[cell, 0].rows.start
    np.testing.assert_allclose(
        instrument * np.sqrt(block_size), sharp, atol=2e-15, rtol=0
    )
    np.testing.assert_array_equal(
        sharp.T @ sharp, block_size * np.eye(2, dtype=np.int64)
    )
    assert len(trial) * np.max(abs(trial)) < 2**60
    a, b, c, d = map(int, (trial.T @ sharp).ravel())
    denom = block_size * (isqrt(norms[0] * norms[1]) + 1)
    det_floor = Fraction(abs(a * d - b * c), denom)
    # Each normalized trial row is projected onto an isometric two-column instrument.
    _, frobenius_bound = sqrt_bounds(
        sum(
            Fraction(x * x + y * y, block_size * n2)
            for x, y, n2 in [(a, b, norms[0]), (c, d, norms[1])]
        )
    )
    singular = max(Fraction(0), det_floor / frobenius_bound - frame_error)
    # Perturb each of two columns by at most frame_error in the determinant.
    determinant = max(
        Fraction(0), det_floor - 2 * frobenius_bound * frame_error - frame_error**2
    )
    return singular, determinant


def separated_gram(path, trial, distance):
    """Exact numerator of W^T P W for the complement of joint spectator vacua."""
    gram = np.full((2, 2), Fraction(0), dtype=object)
    for block in path.blocks.values():
        axes = tuple(i for i in range(path.size) if abs(i - block.cell) >= distance)
        if not axes:
            continue
        words = trial[block.rows].astype(object)
        vac = words.reshape(*block.shape, 2).sum(axis=axes).reshape(-1, 2)
        weight = Fraction(1, int(np.prod([block.shape[i] for i in axes])))
        gram += words.T @ words - weight * (vac.T @ vac)
    return gram


def receiver_witness(path, trial, norms, error):
    """Rational witness on the receiver after tracing every spectator head."""
    scale = 2**40
    end = path.blocks[path.size - 1, 0]
    code_twice = (2 * path.code).astype(np.int64)
    np.testing.assert_array_equal(code_twice, 2 * path.code)
    # The hub is the penultimate tensor axis, immediately before the input.
    tensor = trial[end.rows].reshape(-1, len(code_twice), 2)
    meter = np.einsum("eo,bei->boi", code_twice, tensor).reshape(-1, 4).astype(object)
    choi = (meter.T @ meter).reshape(2, 2, 2, 2).transpose(0, 3, 2, 1).reshape(4, 4)
    denom = 8 * scale**2
    _, directions = np.linalg.eigh(np.asarray(choi / denom, dtype=float))
    probe = np.rint(directions[:, 0] * 2**24).astype(np.int64).astype(object)
    quotient = Fraction(int(probe @ choi @ probe), denom * int(probe @ probe))
    total_error = error
    for norm in norms:
        lo, hi = sqrt_bounds(Fraction(norm, scale**2))
        total_error += max(abs(1 - lo), abs(1 - hi))
    assert total_error < 1
    # Trace norm contracts under the partial trace, and Hilbert–Schmidt
    # norm is preserved by partial transpose. The rank-one outer products
    # before the trace differ in trace norm by at most (2+delta)*delta.
    negativity = max(Fraction(0), -quotient - 3 * total_error)
    return probe.tolist(), quotient, negativity


def certify_channel(path, h, energies, columns, lower_cut, energy_ground_floor):
    trial, norms, centers, radii, error = dyadic_vectors(
        h, energies, columns, lower_cut
    )
    source, source_det = capture_floor(path, 0, trial, norms, error)
    receiver, receiver_det = capture_floor(path, path.size - 1, trial, norms, error)
    probe, quotient, local_witness = receiver_witness(path, trial, norms, error)
    distant = []
    for d in range(1, path.size):
        g = separated_gram(path, trial, d)
        a, b = g[0, 0] / norms[0], g[1, 1] / norms[1]
        disc = (a - b) ** 2 + 4 * g[0, 1] ** 2 / (norms[0] * norms[1])
        upper = (a + b + sqrt_bounds(disc)[1]) / 2
        amp = sqrt_bounds(upper)[1]
        bound = min(
            Fraction(1),
            (amp + error) ** 2,
            centers[-1] + radii[-1] - energy_ground_floor,
        )
        distant.append(float(np.nextafter(float(bound), np.inf)))
    return {
        "dyadic_energy_centers": list(map(str, centers)),
        "eigenvalue_error_ceilings": list(map(str, radii)),
        "eigenvector_frame_error_upper": float(np.nextafter(float(error), np.inf)),
        "source_minimum_singular_value_lower": max(
            0.0, float(np.nextafter(float(source), -np.inf))
        ),
        "vacuum_and_receiver_minimum_singular_value_lower": max(
            0.0, float(np.nextafter(float(receiver), -np.inf))
        ),
        "all_waits_joint_probability_lower_for_every_input": max(
            0.0, float(np.nextafter(float((source * receiver) ** 2), -np.inf))
        ),
        "all_waits_unnormalized_reference_negativity_lower": max(
            0.0, float(np.nextafter(float(source_det * receiver_det / 2), -np.inf))
        ),
        "receiver_only_integer_witness_for_band_input": probe,
        "receiver_only_trial_witness_quotient": str(quotient),
        "receiver_only_band_input_negative_eigenvalue_magnitude_lower": max(
            0.0, float(np.nextafter(float(local_witness), -np.inf))
        ),
        "all_waits_receiver_only_negative_eigenvalue_magnitude_per_launch_lower": max(
            0.0, float(np.nextafter(float(local_witness * source**2), -np.inf))
        ),
        "all_states_separated_probability_upper_by_distance": distant,
    }
