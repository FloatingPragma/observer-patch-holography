"""Outward enclosures of time-dependent angular optics on the native stationary source."""

from fractions import Fraction
from math import comb, isqrt
import argparse
import json

import numpy as np
import sympy as sp
from mpmath import iv

from ..beam_response import measure as beam
from ..full_algebra import certificate as full
from ..full_algebra import native_precision as native
from ..null_optics.certificate import bounded_sources
from . import measure

prior = native.prior


def amplitude_jets(h, senders, receivers, coeff):
    """Evaluate the integer kernel on the rows needed for H^p U, p <= 2."""
    launches, ends = beam.indices(senders), beam.indices(receivers)
    row_ops = [None, h[ends], (h[ends] @ h).tocsr()]
    support = np.unique(np.r_[ends, row_ops[1].indices, row_ops[2].indices])
    dim, columns = h.shape[0], len(launches)
    row = int(abs(h).sum(axis=1).max())
    assert row < 64
    n = prior.DEGREE
    err = Fraction(3 * (1 + isqrt(dim)) * n * (n + 1) * (n + 2), 6 * prior.UNIT)
    err += Fraction(3 * (n + 1), prior.FRAC) + Fraction(1, 10**59)
    assert 4 * err < Fraction(1, 10**19)
    high = np.zeros((dim, columns), dtype=np.int64)
    low = high.copy()
    high[launches] = (prior.UNIT // prior.BASE) * np.eye(columns, dtype=np.int64)
    total = np.zeros((2, len(coeff), len(support), columns), dtype=object)
    prev_hi = prev_lo = None
    for k in range(n + 1):
        assert (
            row * (int(np.max(np.abs(high))) + 1) < 2**62 and row * prior.BASE < 2**62
        )
        val = high[support].astype(object) * prior.BASE + low[support]
        total[k % 2] += (
            (1 if k % 4 in (0, 3) else -1) * coeff[:, k, None, None] * val[None]
        )
        if k < n:
            nh, nl = prior.divided_product(h, high, low, 6 if k == 0 else 5)
            if k:
                nl -= prev_lo
                nh += (nl >> prior.LIMB) - prev_hi
                nl &= prior.MASK
            prev_hi, prev_lo, high, low = high, low, nh, nl
    u = total // prior.UNIT
    jets = []
    for p in range(3):
        if p == 0:
            a, b = u[:, :, np.searchsorted(support, ends)]
        else:
            matrix = row_ops[p][:, support].toarray().astype(object)
            a, b = matrix @ u[0], matrix @ u[1]
            if p == 1:
                a, b = b, -a
            else:
                a, b = -a, -b
        jets.append(
            [
                part.reshape(len(coeff), 5, 4, 5, 4).transpose(0, 1, 3, 2, 4)
                for part in (a, b)
            ]
        )
    return jets


def probability_jets(amplitudes):
    """Exact numerators. Order p has common denominator 2**p * FRAC**2."""
    result = []
    for p in range(3):
        re, im = 0, 0
        for k in range(p + 1):
            a, b = amplitudes[k]
            c, d = amplitudes[p - k]
            at, bt = a.swapaxes(-1, -2), b.swapaxes(-1, -2)
            re += comb(p, k) * (at @ c + bt @ d)
            im += comb(p, k) * (at @ d - bt @ c)
        result.append([re, im])
    return np.array(result, dtype=object)


def translated_basis(probs):
    # Integer translations do not change the central covariance or any time derivative.
    def axes(marginal):
        mean = marginal @ beam.BEARINGS / probs.sum()
        shift = np.rint(mean * beam.UNIT).astype("int64")
        return np.c_[np.full(5, beam.UNIT), beam.PIXELS - shift].astype(object)

    return axes(probs.sum(1)), axes(probs.sum(0))


def enclose_moments(operators, state, source_error, variations=None):
    # Shape: parity, derivative, real/imaginary, receiver, sender, final/initial coin.
    denominator = prior.FRAC**2
    nominal = np.array([np.array(b, dtype=complex) for b in state])
    center = (
        np.array(operators[:, 0, 0], dtype=float) / denominator
        + 1j * np.array(operators[:, 0, 1], dtype=float) / denominator
    )
    out, inp = translated_basis(beam.probabilities(center, nominal))
    data = []
    varied = {name: [] for name in ("ordinary", "unpair", "rotation")}
    for order in range(3):
        denom = 2**order * denominator * beam.UNIT**2
        entries = []
        for left in out.T:
            for right in inp.T:
                weights = np.outer(left, right)
                assert sum(abs(w) for w in weights.ravel()) <= 100 * beam.UNIT**2
                re, im = (
                    (operators[:, order] * weights[None, None, :, :, None, None])
                    .sum((2, 3))
                    .swapaxes(0, 1)
                )
                blocks = [
                    (sp.Matrix(a) + sp.I * sp.Matrix(b)) / denom for a, b in zip(re, im)
                ]
                center = sum(full.inner(h, b) for h, b in zip(blocks, state))
                # Each true derivative of an amplitude block has norm <=32**p.
                # The derivative error of a 4x4 block is <=4*32**p*1e-19
                # from row absolute sums and the pointwise kernel enclosure.
                # Binomial Leibniz and <=100 total absolute bearing weight give this bound.
                kernel_error = sp.Rational(100 * 10 * 64**order, 10**19)
                norm = max(prior.one_norm(h) for h in blocks)
                entries.append(full.ri_band(center, norm * source_error + kernel_error))
                if variations is not None:
                    _, targets, axis, rots, _ = variations

                    def expectation(d, err):
                        return full.ri_band(
                            sum(full.inner(h, b) for h, b in zip(blocks, d)),
                            norm * err + 2 * kernel_error,
                        )

                    varied["ordinary"].append(
                        expectation(targets["ordinary"], 2 * source_error)
                    )
                    varied["unpair"].append(
                        expectation(targets["unpair"], 3 * source_error)
                    )
                    varied["rotation"].append(
                        sum(
                            a * expectation(d, source_error) for a, d in zip(axis, rots)
                        )
                    )
        data.append(np.array(entries, dtype=object).reshape(3, 3))
    if variations is not None:
        return data, {
            name: np.array(vals, dtype=object).reshape(3, 3, 3)
            for name, vals in varied.items()
        }
    return data


def run():
    iv.dps = 32
    f = beam.front.optics.pairing.formation
    board = f.federation(1)
    route = f.region(board)
    rings = beam.screens(board, route)
    coeff = prior.weights()[[int(2 * t) - 1 for t in measure.TIMES]]
    state, rad = bounded_sources(sp.eye(4) / 4)["state"]
    data = {t: {} for t in measure.TIMES}
    for q in (-1, 1):
        effects = []
        for odd in (False, True):
            h = native.generator(board, route, q, odd)
            effects.append(probability_jets(amplitude_jets(h, *rings, coeff)))
        effects = np.array(effects, dtype=object)
        for index, t in enumerate(measure.TIMES):
            moments = enclose_moments(effects[:, :, :, index], state, rad)
            j = measure.matrix_jets(moments)
            det, b, c = measure.optical_matrices(j)
            data[t][q] = (det, b, c)
    assert data[4.0][-1][0].a > 0 and data[4.0][1][0].a > 0
    assert data[4.5][-1][0].b < 0 and data[4.5][1][0].a > 0
    assert data[5.0][-1][0].a > 0 and data[5.0][1][0].a > 0

    def describe(v):
        det, b, c = v
        return dict(
            signed_area=str(det),
            expansion=str(b.trace()),
            twist=str((b[1, 0] - b[0, 1]) / 2),
            shear_diagonal=str((b[0, 0] - b[1, 1]) / 2),
            shear_off_diagonal=str((b[0, 1] + b[1, 0]) / 2),
            trace_focusing_if_affine=str(-c.trace()),
        )

    out = {
        str(t): {str(q): describe(v) for q, v in vals.items()}
        for t, vals in data.items()
    }
    neg, pos = data[4.0][-1], data[4.0][1]
    delta_b, delta_c = neg[1] - pos[1], neg[2] - pos[2]
    ratios = [delta_c[i, j] / delta_b[i, j] for i, j in ((0, 0), (0, 1))]
    assert ratios[0].b < 0 and ratios[1].a > 0
    u_trace = delta_c.trace() / delta_b.trace()
    common_trace = -neg[2].trace() + u_trace * neg[1].trace()
    assert common_trace.b < 0
    return dict(
        clock_times=list(measure.TIMES),
        reference=out,
        trace_only_nonaffinity=str(u_trace),
        shared_trace_times_squared_affine_rate=str(common_trace),
        necessary_common_nonaffinity_11_and_12=[str(v) for v in ratios],
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    s = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if args.write:
        (measure.HERE / "certificate.json").write_text(s)
    else:
        print(s, end="")
