"""Rational enclosures of the native formed source and full reduced strain read."""

import argparse
import json
from math import factorial, isqrt

import numpy as np
import sympy as sp

from ....source_implication.isometric_collar import certify as prior
from . import measure

PRECISION = 2**51
DEGREE = 20
STEPS = 32
Q = sp.Rational


def formed_state(board, route, charge, weights):
    """Centered native H, Chebyshev columns in two integer limbs, t=4."""
    h = prior.integer_flight(board, route, charge)
    row = int(abs(h).sum(axis=1).max())
    assert row <= 64
    d = h.shape[0]
    hi = np.zeros(d, dtype=np.int64)
    left, right = measure.formation.starting_coins(route, charge)
    hi[4 * route[0] + 2 * left + right] = prior.UNIT // prior.BASE
    lo = np.zeros_like(hi)
    sums = np.zeros((2, d), dtype=object)
    old_hi = old_lo = None
    for k, c in enumerate(weights):
        val = hi.astype(object) * prior.BASE + lo
        sums[k % 2] += (1 if k % 4 in (0, 3) else -1) * c * val
        if k == prior.DEGREE:
            break
        assert row * (int(np.max(abs(hi))) + 1) < 2**62 and row * prior.BASE < 2**62
        nh, nl = prior.divided_product(h, hi, lo, 6 if k == 0 else 5)
        if k:
            nl -= old_lo
            nh += (nl >> prior.LIMB) - old_hi
            nl &= prior.MASK
        old_hi, old_lo, hi, lo = hi, lo, nh, nl
    re, im = (v.reshape(-1, 4) for v in sums)
    raw = sp.Matrix(re.T @ re + im.T @ im) + sp.I * sp.Matrix(im.T @ re - re.T @ im)
    rho = raw / sp.trace(raw)
    n = prior.DEGREE
    wave_error = Q(3 * (1 + isqrt(d)) * n * (n + 1) * (n + 2), 6 * prior.UNIT)
    wave_error += Q(3 * (n + 1), prior.FRAC) + Q(1, 10**59)
    # Normalizing an approximate unit vector at distance e costs at most 2e
    # in vector norm. Partial trace contracts trace distance of pure states.
    err = 4 * wave_error
    assert err < Q(1, 10**18)
    left = sp.Matrix(2, 2, lambda i, j: sum(rho[2 * i + k, 2 * j + k] for k in (0, 1)))
    right = sp.Matrix(2, 2, lambda i, j: sum(rho[2 * k + i, 2 * k + j] for k in (0, 1)))
    product = sp.kronecker_product(left, right)
    lower = Q(1, 10**9)
    grad, log_error = prior.entropy_coordinates(rho, err, lower, [product - rho])
    # ||log(true rho)|| < 21 follows from the certified eigenvalue cutoff.
    kerr = 2 * log_error + 3 * err * 21
    return (rho, product), err, (grad[0] - kerr, grad[0] + kerr)


def propagate(a, start):
    """Round-to-nearest Taylor replay, finite rational error per original column."""
    d, count = start.shape
    row = int(abs(a).sum(axis=1).max())
    length = Q(row, 40 * STEPS)
    assert length <= 1 and (a.T != a).nnz == 0
    rounding = Q(3 * DEGREE * (1 + isqrt(2 * d)), 2 * PRECISION)
    # The Hermitian scalar Taylor integral remainder costs L^(k+1)/(k+1)!.
    local = rounding + 2 * length ** (DEGREE + 1) / factorial(DEGREE + 1)
    assert 4 * STEPS * local < Q(1, 1000)
    wave = np.concatenate((start * (PRECISION // 2), np.zeros_like(start)), axis=1)
    limit = 2**63 - 1
    for step in range(1, 4 * STEPS + 1):
        term = wave.copy()
        for k in range(1, DEGREE + 1):
            denom = 40 * STEPS * k
            assert row * int(np.max(abs(term))) + denom < limit
            product = a @ term
            term[:, :count] = (product[:, count:] + denom // 2) // denom
            term[:, count:] = (-product[:, :count] + denom // 2) // denom
            assert int(np.max(abs(wave))) + int(np.max(abs(term))) < limit
            wave += term
        if step / STEPS in measure.TIMES:
            yield step / STEPS, wave, local * step


def exact_gram(w):
    """Integer dot products with row-bounded 18-bit limbs."""
    bits = 18
    low = (1 << bits) - 1
    digits = [w & low, (w >> bits) & low, w >> (2 * bits)]
    peak = max(int(abs(d).max()) for d in digits)
    assert w.shape[0] * peak**2 < 2**62
    return sum(
        (a.T @ b).astype(object) * (1 << (bits * (i + j)))
        for i, a in enumerate(digits)
        for j, b in enumerate(digits)
    )


def effects(wave):
    result = []
    for rail in (0, 1):
        # Four groups of head columns: real and imaginary parts on +Z, then -Z.
        r, i = np.split(wave, 2, axis=1)
        values = [
            u[:, 4 * rail : 4 * rail + 4]
            for q in np.split(np.column_stack((r, i)), 2)
            for u in np.split(q, 2, axis=1)
        ]
        g = sp.Matrix(exact_gram(np.column_stack(values))) / PRECISION**2

        def block(a, b):
            return g[4 * a : 4 * a + 4, 4 * b : 4 * b + 4]

        def cross(a, b):
            return (
                block(2 * a, 2 * b)
                + block(2 * a + 1, 2 * b + 1)
                + sp.I * (block(2 * a, 2 * b + 1) - block(2 * a + 1, 2 * b))
            )

        off = cross(1, 0)
        result.append((off + off.H, sp.I * (off - off.H), cross(0, 0) - cross(1, 1)))
    return result


def hermitian_trace(effect, state):
    """Hilbert–Schmidt pairing, all scalar products over rational reals."""
    total = sum(effect[i, i] * state[i, i] for i in range(4))
    for i in range(4):
        for j in range(i):
            a, b = effect[i, j].as_real_imag()
            c, d = state[i, j].as_real_imag()
            total += 2 * (a * c + b * d)
    return total


def enclosed(value, error):
    return [
        float(np.nextafter(float(value - error), -np.inf)),
        float(np.nextafter(float(value + error), np.inf)),
    ]


def certificate(nominal):
    weights = prior.weights()[7]
    source = []
    log_two_center = Q(float(np.log(2)))
    prior.log_bound(sp.Integer(2), log_two_center)
    checks = []
    for level in (0, 1):
        board = measure.formation.federation(level)
        route = measure.formation.region(board)
        states = {}
        for sign in (-1, 1):
            pair, state_error, (lo, hi) = formed_state(board, route, sign, weights)
            states[sign] = (pair, state_error, (lo, hi))
            nominal_rho = prior.sym_matrix(measure.form(board, route, sign)[0])
            float_state_error = state_error + 4 * prior.one_norm(nominal_rho - pair[0])
            source.append(
                dict(
                    cells=board.cells,
                    forming_charge=sign,
                    collar_modular_nats=[enclosed(lo, 0)[0], enclosed(hi, 0)[1]],
                    formed_trace_norm_error=enclosed(state_error, 0)[1],
                    nominal_formed_trace_norm_error=enclosed(float_state_error, 0)[1],
                )
            )
        for q in (-1, 1):
            a, _, _ = measure.system(board, route, q)
            start = measure.launch(board, route)
            floating = measure.fast_flow(a, start / 2)
            for time, wave, col_error in propagate(a, start):
                t, target = next(floating)
                assert t == time and int(abs(wave).max()) < 2**53
                real, imag = np.split(wave, 2, axis=1)
                peak = max(
                    abs(target.real - real / PRECISION).max(),
                    abs(target.imag - imag / PRECISION).max(),
                )
                # Integer centers convert exactly. Twice the observed component
                # error plus one subnormal bounds one binary64 subtraction.
                float_col_error = col_error + (1 + isqrt(2 * wave.shape[0])) * (
                    2 * Q(float(peak)) + Q(1, 2**1074)
                )
                ops = effects(wave)
                for s, (pair, state_error, (lo, hi)) in states.items():
                    row = next(
                        r
                        for r in nominal["runs"]
                        if (
                            r["level"],
                            r["forming_charge"],
                            r["pump_charge"],
                            r["time"],
                        )
                        == (level, s, q, time)
                    )
                    arm_bounds = {}
                    for name, numerators in (
                        ("reference", (2, 1)),
                        ("ordinary", (1, 2)),
                        ("decorrelated", (2, 1)),
                    ):
                        rho = pair[name == "decorrelated"]
                        # Four coherent head columns per rail; both rails mix classically.
                        delta = 2 * col_error
                        err = (
                            2 * delta
                            + delta**2
                            + (2 if name == "decorrelated" else 1) * state_error
                        )
                        moment_bounds = []
                        for j in range(3):
                            v = sum(
                                Q(w, 3) * hermitian_trace(ops[rail][j], rho)
                                for rail, w in enumerate(numerators)
                            )
                            assert sp.im(v) == 0
                            archived = Q(row["arms"][name]["bloch"][j])
                            moment_bounds.append(abs(v - archived) + err)
                        arm_bounds[name] = enclosed(max(moment_bounds), 0)[1]

                    def arm(name):
                        return [Q(v) for v in row["arms"][name]["bloch"]]

                    o = [u - v for u, v in zip(arm("ordinary"), arm("reference"))]
                    c = [u - v for u, v in zip(arm("decorrelated"), arm("reference"))]
                    oe = Q(arm_bounds["reference"]) + Q(arm_bounds["ordinary"])
                    ce = Q(arm_bounds["reference"]) + Q(arm_bounds["decorrelated"])
                    minors = [
                        (
                            abs(o[i] * c[j] - o[j] * c[i])
                            - (abs(o[i]) + abs(o[j])) * ce
                            - (abs(c[i]) + abs(c[j])) * oe
                            - 2 * oe * ce
                        )
                        for i, j in ((0, 1), (1, 2), (2, 0))
                    ]
                    assert max(minors) > 0
                    # ||o||_1 upper also bounds ||o||_2 in the distance formula.
                    distance = max(minors) / (sum(abs(v) for v in o) + 3 * oe)
                    mlow = 3 * lo / (log_two_center + Q(1, 10**12))
                    mhigh = 3 * hi / (log_two_center - Q(1, 10**12))
                    mid, spread = (mlow + mhigh) / 2, (mhigh - mlow) / 2
                    residual = max(
                        abs(ci - mid * oi) - ce - mhigh * oe - spread * abs(oi)
                        for ci, oi in zip(c, o)
                    )
                    assert residual > 0
                    checks.append(
                        dict(
                            level=level,
                            forming_charge=s,
                            pump_charge=q,
                            time=time,
                            moment_absolute_errors=arm_bounds,
                            fixed_point_column_error=enclosed(col_error, 0)[1],
                            nominal_column_error=enclosed(float_col_error, 0)[1],
                            collar_distance_from_ordinary_line_lower=enclosed(
                                distance, 0
                            )[0],
                            modular_residual_max_component_lower=enclosed(residual, 0)[
                                0
                            ],
                        )
                    )
    return dict(
        fractional_bits=51,
        steps_per_tick=STEPS,
        degree=DEGREE,
        source=source,
        checks=checks,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    data = certificate(json.loads((measure.HERE / "results.json").read_text()))
    result = json.dumps(data, indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "certificate.json").write_text(result)
    else:
        print(result, end="")
