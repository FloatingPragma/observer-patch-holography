"""Cubic weak-probe contrast of opposite stationary first-character sources.

In the rotating frame of the common source energy E, the first outside
amplitude is -i integral(exp(-i(H_out-E)s), s=0..T) QVP. Its squared
effect is real on the source since H_out, V and Q are real. Conjugate
angular sources have exactly equal quadratic outside probabilities.
"""

import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import factorial, inf, nextafter, sqrt
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix, diags, eye, kron

from .....replay import compare
from ..measure import cohorts, matched_axis
from .apron import instrument
from .exact_band import certify, negative_pivots, point, square_root
from .wheel import LAND, RIM, preparation, rotation, wheel_basis

WAIT = 4
ORDER = 10
COEFFICIENT_SCALE = 8
STEP = Fraction(1, 16)
POLY = 26
RATES = tuple(map(Fraction, ("0.001", "0.01", "0.025", "0.05")))
HERE = Path(__file__).resolve().parent


def angular_input(basis, stock_h):
    energy, _ = certify()
    root5 = square_root(point(5))
    cosine = (root5 - point(1)) / point(4)
    sine = square_root((point(5) + root5) / point(8))
    z = float((cosine.lo + cosine.hi) / 2) - 1j * float((sine.lo + sine.hi) / 2)
    # The six disjoint orbits realize the positive angular character.
    seeds = np.array([i for i, (x, _) in enumerate(basis) if x == RIM[0]])
    frame = np.zeros((len(basis), len(seeds)), complex)
    rotate = rotation(basis)
    support = seeds
    for j in range(len(RIM)):
        frame[support, np.arange(len(seeds))] = (z**j) / sqrt(5)
        support = rotate[support]
    np.testing.assert_array_equal(support, seeds)
    _, vecs = np.linalg.eigh(frame.conj().T @ stock_h @ frame)
    source = frame @ vecs[:, 0]
    source /= np.linalg.norm(source)
    center = float((energy.lo + energy.hi) / 2)
    assert abs(np.linalg.norm(source) - 1) < 4e-15
    assert np.linalg.norm(stock_h @ source - center * source) < 1e-12
    # Exact Sturm count at 8.8 proves an in-character gap above 1/4.
    assert negative_pivots(Fraction("8.8"), point("15.25") + root5) == 5
    assert Fraction("8.8") - energy.hi > Fraction(1, 4)
    # At most 5 powers and one scale per frame entry. A 1e-13 frame
    # error plus the residual and gap bounds gives this outward
    # source-vector allowance, after aligning unobservable global phase.
    source_error = 1e-10
    return np.column_stack((source, source.conj())), source_error


def gamma(n):
    return Fraction(n, 2**53 - n)


def cascade(h0, interaction, source):
    """Solve the triangular Duhamel cascade with dyadic Taylor stepping."""
    n = h0.shape[0]
    assert n == 46 and ORDER * n < 512
    shift = diags(np.ones(ORDER), -1, shape=(ORDER + 1, ORDER + 1))
    h = (
        kron(eye(ORDER + 1), csr_matrix(h0))
        + kron(shift, csr_matrix(interaction / COEFFICIENT_SCALE))
    ).tocsr()
    h.eliminate_zeros()
    q = h * 8
    np.testing.assert_array_equal(q.data, np.rint(q.data))
    width = max(np.diff(q.indptr))
    assert width <= 8
    row_sum = int(max(np.asarray(abs(q).sum(axis=1)).ravel()))
    x = Fraction(row_sum, 8) * STEP
    assert x < 2
    # Each real sparse row acts on two components. Coefficient division,
    # phase, and accumulation use the same componentwise gamma estimate.
    weights = [x**j / factorial(j) for j in range(POLY + 1)]
    roundoff = (1 + gamma(POLY + 1)) * sum(
        weight * (1 + gamma(2 * int(width) + 8)) ** j
        for j, weight in enumerate(weights)
    ) - sum(weights)
    tail = 8 * x ** (POLY + 1) / factorial(POLY + 1)
    ticks = int(WAIT / STEP)
    err_per_tick = Fraction(23) * (roundoff + tail + Fraction(2) ** -990)
    # sqrt(506)<23 converts componentwise bounds to two-norm for each
    # complex column. The complex componentwise budget is doubled below.
    # The exact cascade has norm growth <=exp(2T/8)=exp(1)<3.
    vector_error = float(18 * ticks * err_per_tick)
    assert vector_error < 1e-8
    state = np.zeros(((ORDER + 1) * n, 2), complex)
    state[:n] = source
    for _ in range(ticks):
        term = state
        total = state.copy()
        for j in range(1, POLY + 1):
            term = (-1j * float(STEP) / j) * (h @ term)
            total += term
        state = total
    return state.reshape(ORDER + 1, n, 2), vector_error


def outward(center, pad):
    return [nextafter(float(center - pad), -inf), nextafter(float(center + pad), inf)]


def archived_reference(basis, route, src, rows):
    path = HERE.parents[1] / "word_reads.json"
    raw = path.read_bytes()
    case = next(c for c in json.loads(raw)["cases"] if c["fan_count"] == 4)
    stock = next(a for a in case["stock_arms"] if a["preparation"] == "clock_bypassed")
    late = next(
        f
        for f in stock["flights"]
        if f["appointment"] == 30.8125 and f["material_bridge_rate"] == 4
    )
    pair = next(
        p
        for p in late["pair_readings"]
        if p["port_face_and_words"] == [1045, 0, 334, 365]
    )
    pauli = np.asarray(pair["receiver_reference_pauli_moments"])
    reports = cohorts(pauli, matched_axis(pauli))
    restricted = src[[basis.index((LAND, w)) for w in route["wheel_input_words"]]]
    assert np.linalg.norm(restricted, ord=2) ** 2 < 0.24
    # Each 1/4-weighted Pauli term changes by at most (1+sqrt(3))*epsilon.
    # Trace norm on the two-word input is at most twice the spectral norm.
    # The exact projection map contracts it by at most 1/4, with extra
    # room for the eigensource bound and its finite matrix products.
    raw_error = 2 * (1 + sqrt(3)) * late["observable_absolute_error_upper_per_pair"]
    error = raw_error / 4 + 1e-11
    normalized = {}
    admission = {}
    normalized_error = 0.0
    for sign in ("plus", "minus"):
        rho = restricted.conj().T @ reports[sign][0] @ restricted
        passed = float(rho.trace().real)
        admission[sign] = outward(passed, error)
        assert passed > error
        normalized[sign] = rho / passed
        normalized_error += 2 * error / (passed - error)
    difference = normalized["plus"] - normalized["minus"]
    scale = float(difference[0, 0].real)
    assert 0.2 < scale < 0.21
    assert np.linalg.norm(difference - np.diag([scale, -scale]), ord="nuc") < 1e-12
    result = []
    for row in rows:
        # The outside effect has norm <=(2 eta T)^2 on the original words.
        # This bounds the induced probability error from both density errors.
        pad = row["source_disturbance_upper"] * (normalized_error + 1e-12)
        lo, hi = row["k_plus_one_minus_k_minus_one_outside_probability"]
        bounds = [nextafter(scale * lo - pad, -inf), nextafter(scale * hi + pad, inf)]
        if row["exterior_face_rate"] >= 0.01:
            assert bounds[1] < 0
        result.append(
            {
                "exterior_face_rate": row["exterior_face_rate"],
                "conditional_reference_plus_minus_outside_probability": bounds,
                "source_calibration_allowance": pad,
            }
        )
    return {
        "source_word_reads_sha256": sha256(raw).hexdigest(),
        "source_wait": 30.8125,
        "doublet_admission_per_original_start": admission,
        "nominal_k_conjugate_contrast_multiplier": scale,
        "source_density_trace_norm_allowance_in_contrast": normalized_error,
        "conditional_reference_runs": result,
    }


def run():
    basis, route = wheel_basis(), preparation()
    h0, embed, changed, _ = instrument(basis, route, 0)
    h1, check_embed, check_changed, _ = instrument(basis, route, 1)
    assert embed == check_embed and len(h0) == 46
    np.testing.assert_array_equal(changed, check_changed)
    np.testing.assert_array_equal(h0[np.ix_(changed, ~changed)], 0)
    v = h1 - h0
    np.testing.assert_array_equal(v, v.T)
    np.testing.assert_array_equal(v @ v, 2 * v)
    assert np.any(v) and np.all(np.sum(abs(v), axis=1) <= 2)
    src, source_error = angular_input(basis, h0[np.ix_(embed, embed)])
    initial = np.zeros((len(h0), 2), complex)
    initial[embed] = src
    stacked, tick_error = cascade(h0, v, initial)
    out = stacked[:, changed]
    pad = 1.0001 * tick_error + 3 * source_error
    assert np.linalg.norm(out[0]) < pad
    contrasts = {}
    probabilities = {}
    for order in range(2, ORDER + 1):
        coeff = sum(
            np.sum((out[j].conj() * out[order - j]).real, axis=0)
            for j in range(1, order)
        )
        # Exact scaled amplitudes have norm <=1/j! since T ||V||/8=1.
        radius = (
            sum(
                pad * (1 / factorial(j) + 1 / factorial(order - j)) + pad**2
                for j in range(1, order)
            )
            + 1e-14
        )
        probabilities[order] = [outward(c, radius) for c in coeff]
        contrasts[order] = outward(coeff[0] - coeff[1], 2 * radius)
    assert contrasts[2][0] < 0 < contrasts[2][1]
    contrasts[2] = [0.0, 0.0]
    cubic = [q * COEFFICIENT_SCALE**3 for q in contrasts[3]]
    assert -0.7711 < cubic[0] <= cubic[1] < -0.7709
    rows = []
    for eta in RATES:
        alpha = COEFFICIENT_SCALE * eta
        a = 2 * alpha
        # Sum both tails, using the ratio bound from degree ORDER+1.
        radius = 2 * a ** (ORDER + 1) / factorial(ORDER + 1) / (1 - a / (ORDER + 2))
        low, high = (
            sum(contrasts[m][i] * float(alpha**m) for m in range(3, ORDER + 1))
            for i in (0, 1)
        )
        enclosure = [
            nextafter(low - float(radius) - 1e-17, -inf),
            nextafter(high + float(radius) + 1e-17, inf),
        ]
        assert enclosure[1] < 0
        rows.append(
            {
                "exterior_face_rate": float(eta),
                "k_plus_one_minus_k_minus_one_outside_probability": enclosure,
                "absolute_tail_beyond_tenth_order": float(radius),
                "source_disturbance_upper": float(alpha**2),
            }
        )
    return {
        "enrolled_configurations": len(h0),
        "wait": WAIT,
        "coefficient_scale": COEFFICIENT_SCALE,
        "quadratic_outside_coefficient_per_character": [
            [q * COEFFICIENT_SCALE**2 for q in a] for a in probabilities[2]
        ],
        "k_plus_one_minus_k_minus_one_cubic_coefficient_in_face_rate": cubic,
        "polynomial_coefficients_in_eight_times_rate": [
            {"order": m, "contrast": contrasts[m]} for m in range(2, ORDER + 1)
        ],
        "conditional_angular_character_runs": rows,
        "archived_late_reference": archived_reference(basis, route, src, rows),
    }


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    args = cli.parse_args()
    rows = run()
    path = HERE / "weak_probe.json"
    if args.write:
        path.write_text(json.dumps(rows, indent=2, sort_keys=True) + "\n")
    else:
        compare(rows, json.loads(path.read_text()), "weak_probe")
        print(
            "Verified exact quadratic cancellation and certified cubic weak-probe contrast."
        )


if __name__ == "__main__":
    main()
