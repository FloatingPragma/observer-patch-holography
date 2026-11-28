"""Resolved low-energy prefixes without resolving the admitted modes during flight."""

import argparse
import json
from fractions import Fraction
from itertools import product
from math import ceil, factorial, floor, isqrt, pi
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.linalg import eigh

from model.candidate_laws.interval_echo.openwave.electron_field.fan_path.spectral import (
    verified_inertia,
)
from model.candidate_laws.interval_echo.openwave.electron_field.fan_rotation.coherent_filter import (
    pi_interval,
)

from ..measure import (
    ACCEPTANCE,
    WAITS,
    geometry,
    local_certificate,
    native_generator,
    source_slots,
)

TARGETS = (1, 2, 4, 8, 16, 32)
PARITIES = tuple(product((1, -1), repeat=2))
PREPARED = (1, 1)
VECTOR_GRID = 1 << 53
ENERGY_GRID = 1 << 48
PHASE_GRID = 1 << 60
CERTIFICATE_WAIT = Fraction(32)


def sqrt_upper(value):
    value = Fraction(value)
    grid = 1 << 80
    return Fraction(isqrt(value.numerator * grid**2 // value.denominator) + 1, grid)


def upper_float(x):
    return float(np.nextafter(float(x), np.inf))


def lower_float(x):
    return float(np.nextafter(float(x), -np.inf))


def band_cuts(matrices, levels):
    spectra = sorted((ev, parity) for parity, bank in levels.items() for ev in bank)
    result = []
    for target in TARGETS:
        last = target
        while spectra[last][0] - spectra[last - 1][0] < 1e-8:
            last += 1
        a, b = (spectra[last - 1][0], spectra[last][0])
        low, high = (
            Fraction(round((a + (b - a) * f) * 2**25), 2**25)
            for f in (Fraction(1, 3), Fraction(2, 3))
        )
        assert a < low < high < b
        counts = {
            parity: [
                verified_inertia(matrices[parity] / 2, cut)["count_below"]
                for cut in (low, high)
            ]
            for parity in PARITIES
        }
        assert all(a == b for a, b in counts.values())
        assert sum(v[0] for v in counts.values()) == last
        result.append((target, low, high, {p: n[0] for p, n in counts.items()}))
    return result


def certify_frame(h2, values, vectors, high):
    """Enclose the spectral projector and its phase evolution with integer residuals."""
    vectors_int = np.rint(vectors * VECTOR_GRID).astype(np.int64)
    energies = [Fraction(round(x * ENERGY_GRID), ENERGY_GRID) for x in values]
    assert np.max(abs(vectors_int)) <= VECTOR_GRID
    assert np.max(np.asarray(abs(h2).sum(axis=1))) * VECTOR_GRID < 2**62
    shifted = h2 @ vectors_int
    residual_squared = 0
    for j, e in enumerate(energies):
        eigen_int = int(e * ENERGY_GRID)
        residual_squared += sum(
            (int(hv) * ENERGY_GRID - 2 * eigen_int * int(v)) ** 2
            for hv, v in zip(shifted[:, j], vectors_int[:, j], strict=True)
        )
    residual = sqrt_upper(
        Fraction(residual_squared, (2 * VECTOR_GRID * ENERGY_GRID) ** 2)
    )
    vectors_exact = vectors_int.astype(object)
    gram = vectors_exact.T @ vectors_exact
    for j in range(len(values)):
        gram[j, j] -= VECTOR_GRID**2
    gram_error = sqrt_upper(
        Fraction(sum(int(x) ** 2 for x in gram.flat), VECTOR_GRID**4)
    )
    assert gram_error < Fraction(1, 2)
    separation = high - max(energies)
    assert separation > 0
    projection_error = gram_error + 2 * residual / separation
    return vectors_exact, energies, residual, projection_error


def integer_phase(angle):
    """Approximate exp(-i angle) with a dyadic Gaussian integer and rational error."""
    degree = 20
    squarings = 0
    while abs(angle) > Fraction(1, 2):
        angle /= 2
        squarings += 1
    re = sum(
        (-1) ** j * angle ** (2 * j) / factorial(2 * j) for j in range(1 + degree // 2)
    )
    im = sum(
        -((-1) ** j) * angle ** (2 * j + 1) / factorial(2 * j + 1)
        for j in range(1 + (degree - 1) // 2)
    )
    a, b = round(re * PHASE_GRID), round(im * PHASE_GRID)
    error = 2 * Fraction(1, 2) ** (degree + 1) / factorial(degree + 1) + Fraction(
        1, PHASE_GRID
    )
    for _ in range(squarings):
        a, b = (
            round(Fraction(a * a - b * b, PHASE_GRID)),
            round(Fraction(2 * a * b, PHASE_GRID)),
        )
        error = error * (2 + error) + Fraction(1, PHASE_GRID)
        # Round the error upwards to prevent denominator growth.
        grid = 1 << 100
        error = Fraction(-(-(error.numerator * grid) // error.denominator), grid)
    return a, b, error


def integer_channel(
    vectors, energies, residual, projection_error, source, capture, time
):
    phases = [integer_phase(time * (e - energies[0])) for e in energies]
    re, im, phase_error = zip(*phases, strict=True)
    left, right = vectors[capture, :], vectors[source, :]
    real = (left * np.array(re, dtype=object)) @ right.T
    imag = (left * np.array(im, dtype=object)) @ right.T
    scale = VECTOR_GRID**2 * PHASE_GRID
    k = sp.Matrix(
        [
            [int(a) + sp.I * int(b) for a, b in zip(r, s, strict=True)]
            for r, s in zip(real, imag, strict=True)
        ]
    )
    unscaled = sp.zeros(4, 4)
    for mark in range(2):
        vector = k[mark::2, :].T.reshape(4, 1)
        unscaled += vector * vector.H
    raw = unscaled.applyfunc(sp.expand) / (8 * scale**2)
    # Propagate F itself with Duhamel, then add projector and phase errors.
    # Our exact Gram bound gives ||F|| < 2 and ||F||^2 < 2.
    delta = projection_error + 2 * time * residual + 2 * max(phase_error)
    local_norm = sqrt_upper(
        Fraction(
            sum(
                int(a) ** 2 + int(b) ** 2
                for a, b in zip(real.flat, imag.flat, strict=True)
            ),
            scale**2,
        )
    )
    branch_error = 2 * delta  # Restriction to two input columns and four output rows.
    error = branch_error * (2 * local_norm + branch_error) / 8
    return raw, error


def partial_transpose(raw):
    return sp.Matrix(4, 4, lambda i, j: raw[(j // 2) * 2 + i % 2, (i // 2) * 2 + j % 2])


def exact_witness(raw, zero_raw):
    pt = partial_transpose(raw)
    _, eigenvectors = np.linalg.eigh(np.asarray(pt.evalf(), dtype=complex))
    w = eigenvectors[:, 0]
    idx = np.argmax(abs(w))
    w = w * (w[idx].conjugate() / abs(w[idx]))
    probe_grid = 1 << 24
    ints = np.rint(probe_grid * np.column_stack((w.real, w.imag))).astype(int)
    vector = sp.Matrix([int(a) + sp.I * int(b) for a, b in ints])
    q = sp.expand((vector.H * vector)[0])
    at_time = Fraction(sp.expand((vector.H * pt * vector)[0] / q))
    at_zero = Fraction(
        sp.expand((vector.H * partial_transpose(zero_raw) * vector)[0] / q)
    )
    return ints.tolist(), at_time, at_zero


def check_witness(current, zero):
    raw, err = current
    raw0, err0 = zero
    integers, w, w0 = exact_witness(raw, raw0)
    return {
        "witness_vector_real_and_imag_integer": integers,
        "wait_witness_interval": [lower_float(w - err), upper_float(w + err)],
        "same_witness_at_zero_interval": [
            lower_float(w0 - err0),
            upper_float(w0 + err0),
        ],
        "wait_minus_zero_interval": [
            lower_float(w - w0 - err - err0),
            upper_float(w - w0 + err + err0),
        ],
        "joint_trace_error_upper": upper_float(err),
    }


def approximate_banks(v, e, source, capture, waits):
    admitted = ACCEPTANCE * float(np.linalg.norm(v[source, :], "fro") ** 2) / 2
    rows = {}
    baseline = None
    for wait in waits:
        block = (v[capture] * np.exp(-1j * float(wait) * (e - e[0]))) @ v[source].T
        joint = np.zeros((4, 4), dtype=complex)
        for mark in range(2):
            f = block[mark::2].T.reshape(-1)
            joint += np.outer(f, f.conj()) / 8
        pt = joint.reshape(2, 2, 2, 2).transpose(2, 1, 0, 3).reshape(4, 4)
        if baseline is None:
            baseline = joint
        arrival = float(np.trace(joint).real)
        assert arrival <= admitted + 1e-10
        rows[str(wait)] = {
            "write_refusal": 1 - ACCEPTANCE,
            "coarse_window_refusal": round(ACCEPTANCE - admitted, 13),
            "in_window_missing_port": round(admitted - arrival, 13),
            "signed_port_and_coin_click": round(arrival, 13),
            "partial_transpose_minimum_numeric": float(
                f"{np.linalg.eigvalsh(pt)[0]:.7g}"
            ),
            "local_joint_record_change_frobenius_numeric": float(
                f"{np.linalg.norm(joint - baseline):.7g}"
            ),
        }
    return rows


def finite_coarse_filter(matrices, low, high):
    """Coherent majority phase estimation and erase, using only a yes/no flag."""
    clock = Fraction(1, 64)
    bits, repeats = 22, 51
    grid = 1 << bits
    for h2 in matrices.values():
        diagonal = h2.diagonal()
        radius = np.asarray(abs(h2).sum(axis=1)).ravel() - abs(diagonal)
        assert min(diagonal - radius) >= 0 and max(diagonal + radius) <= 256
    margin = min(Fraction(1), (high - low) / 2)
    cutoff = (low + high) / 2
    pi0, pi1 = pi_interval()
    assert 128 * clock < pi0
    bins = floor(grid * margin * clock / (2 * pi1))
    assert bins >= 2
    first = ceil(-grid * clock / (2 * pi0))
    last = floor(grid * clock * cutoff / (2 * pi0))
    assert first == ceil(-grid * clock / (2 * pi1))
    assert last == floor(grid * clock * cutoff / (2 * pi1))
    probability_one = Fraction(1, 2 * (bins - 1))
    filter_error = 2**repeats * probability_one ** ((repeats + 1) // 2)
    assert filter_error < Fraction(1, 10**20)
    return {
        "signed_phase_read_bins_accepting_below_coarse_cut": [first, last],
        "maximum_wrong_bin_probability_one_read": str(probability_one),
        "coherent_phase_bits_each": bits,
        "coherent_read_repetitions": repeats,
        "controlled_native_H_exposure_elementary": str(clock),
        "maximum_clock_exposure_for_one_query": str((grid // 2) * clock),
        "total_absolute_controlled_exposure": str(2 * repeats * (grid - 1) * clock),
        "successful_operator_error_upper": upper_float(filter_error),
        "executed_clock_circuit": False,
    }, filter_error


def certify_two(matrices, low, high, e, v, source, capture, waits):
    exact, energies, residual, frame_bound = certify_frame(
        matrices[PREPARED], e, v, high
    )
    clock, filter_error = finite_coarse_filter(matrices, low, high)
    values = {
        t: integer_channel(exact, energies, residual, frame_bound, source, capture, t)
        for t in (Fraction(0), *waits)
    }
    window_tests = {
        str(t): check_witness(values[t], values[Fraction(0)]) for t in waits
    }
    for summary in window_tests.values():
        assert summary["wait_witness_interval"][1] < -filter_error
        assert summary["same_witness_at_zero_interval"][0] > filter_error
    zero_test = check_witness(values[Fraction(0)], values[Fraction(0)])
    assert zero_test["wait_witness_interval"][1] < 0
    halfbeat = waits[-1]
    lo, hi = pi_interval()
    angle = halfbeat * (energies[1] - energies[0])
    # On distinct disjoint intervals each eigenvalue lies within twice the
    # corresponding column residual, by the spectral theorem and ||F_j||>1/2.
    assert energies[1] - energies[0] > 4 * residual
    assert hi - Fraction(1, 1024) < angle - 4 * halfbeat * residual
    assert angle + 4 * halfbeat * residual < lo + Fraction(1, 1024)
    return {
        "prepared_energy_centers_dyadic": [str(x) for x in energies],
        "each_energy_error_upper": upper_float(2 * residual),
        "integer_spectral_residual_upper": upper_float(residual),
        "projector_operator_error_upper": upper_float(frame_bound),
        "zero_wait_own_negative_witness": zero_test,
        "fixed_positive_wait_witnesses_and_zero_clock_comparison": window_tests,
        "calibrated_near_half_beat_wait": str(halfbeat),
        "finite_coarse_filter_certificate": clock,
    }


def run():
    local_certificate()
    board, launch, receive, marked = geometry()
    capture = [(s + 2 * r) * board.n + receive for r, s in product(range(2), repeat=2)]
    banks = {}
    for sign in (-1, 1):
        source = source_slots(board, launch, receive, marked, sign)
        matrices = {p: native_generator(board, marked, sign, p) for p in PARITIES}
        spectra = {
            p: eigh(h.toarray() / 2, subset_by_index=(0, max(TARGETS) + 4))
            for p, h in matrices.items()
        }
        cuts = band_cuts(matrices, {p: data[0] for p, data in spectra.items()})
        energies, frame = spectra[PREPARED]
        halfbeat = Fraction(round(1024 * pi / (energies[1] - energies[0])), 1024)
        windows = {}
        for target, low, high, ranks in cuts:
            occupied = ranks[PREPARED]
            e, v = energies[:occupied], frame[:, :occupied]
            data = {
                "full_rank": sum(ranks.values()),
                "each_parity_rank": {str(p): r for p, r in ranks.items()},
                "empty_dyadic_interval": [str(low), str(high)],
                "binary_energy_cutoff_anywhere_in_gap": str((low + high) / 2),
                "readouts": approximate_banks(
                    v, e, source, capture, (*WAITS, halfbeat)
                ),
            }
            if target == 1:
                assert occupied == 1 and sum(ranks.values()) == 1
                data["local_reference_entanglement_and_wait_change_exact"] = False
            if target == 2:
                assert occupied == 2 and sum(ranks.values()) == 2
                data["certificate"] = certify_two(
                    matrices,
                    low,
                    high,
                    e,
                    v,
                    source,
                    capture,
                    (CERTIFICATE_WAIT, halfbeat),
                )
                data["individual_energy_record_comparator"] = (
                    "separable and stationary at every wait"
                )
            windows[str(target)] = data
        banks[str(sign)] = windows
    return {
        "schema": "oph.material_spin_band.low_window.v1",
        "predeclared_whole_space_rank_targets": TARGETS,
        "populated_parity": PREPARED,
        "charge_and_coin_read_address": receive,
        "same_wait_banks": WAITS,
        "marks": marked,
        "signs": banks,
    }


def main():
    if not __debug__:
        raise RuntimeError(
            "Run without optimization so certificate assertions execute."
        )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="store_true")
    target = Path(__file__).with_name("results.json")
    data = json.dumps(run(), indent=2, sort_keys=True) + "\n"
    if parser.parse_args().record:
        target.write_text(data)
    elif target.read_text() != data:
        raise RuntimeError("Independent replay differs from the archived receipt.")
    print("Verified global cut counts and lowest two-mode local-reference witness.")


if __name__ == "__main__":
    main()
