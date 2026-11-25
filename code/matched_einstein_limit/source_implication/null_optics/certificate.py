"""Rational charge-resolved exit-clock witnesses over the exact station enclosure."""

import argparse
from fractions import Fraction
import json

import numpy as np
import sympy as sp

from ..full_algebra import native_precision as precision
from .matched import certified_geometry
from .measure import HERE, PAIRS, SEND, RECEIVE, pairing, source


ABSOLUTE = (1, -300, 299, -171, 299, -128, -118, 289, -171, 289, -290, 1)
CHECKER = (1, -1, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1)
PRIOR = precision.prior


def clock_effect(kernel, i, j, cap):
    arrivals = []
    coefficient = np.zeros((2, 4, 4), dtype=object)
    for n in range(1, cap):
        first = kernel[:, n - 1, j, i].copy()
        for k in range(1, n):
            first -= PRIOR.multiply(kernel[:, n - k - 1, j, j], arrivals[k - 1])
        arrivals.append(first)
        coefficient += (n - cap) * PRIOR.multiply(PRIOR.adjoint(first), first)
    operator = sp.Matrix(coefficient[0]) + sp.I * sp.Matrix(coefficient[1])
    return sp.eye(4) * sp.Rational(cap, 2) + operator / (2 * PRIOR.FRAC)


def effect_error(cap):
    multiplication = Fraction(8, PRIOR.FRAC)
    kernel = Fraction(1, 10**19) + multiplication
    f = [Fraction(0)]
    error = Fraction(0)
    for n in range(1, cap):
        f.append(
            kernel
            + sum(kernel + (1 + kernel) * f[j] + multiplication for j in range(1, n))
        )
        error += Fraction(cap - n, 2) * (2 * f[n] + f[n] ** 2 + multiplication)
    assert error < Fraction(1, 10**13)
    return sp.Rational(error)


def marginals(x):
    return (
        sp.Matrix(2, 2, lambda i, j: sum(x[2 * i + k, 2 * j + k] for k in range(2))),
        sp.Matrix(2, 2, lambda i, j: sum(x[2 * k + i, 2 * k + j] for k in range(2))),
    )


def bounded_sources(formed):
    # Same rational, normalized center as the full-algebra existence certificate.
    mid = [PRIOR.sym_matrix(b) for b in source()["state"]]
    mid = [(b + b.H) / 2 for b in mid]
    mid = [b / sum(a.trace() for a in mid) for b in mid]
    right = marginals(sum(mid, sp.zeros(4)))[1]
    delta = [sp.kronecker_product(marginals(b)[0], right) - b for b in mid]
    single = formed - sp.kronecker_product(*marginals(formed))
    # sqrt(8)*2e-7 bounds the station's trace-distance error in trace norm.
    return dict(
        state=(mid, sp.Rational("6e-7")),
        unpair=(delta, sp.Rational("18e-7")),
        dephased_unpair=([sum(delta, sp.zeros(4)) / 2] * 2, sp.Rational("18e-7")),
        formed_pairing=([single, sp.zeros(4)], sp.Rational("1e-18")),
    )


def leading_pair_transfer(board, route):
    # The centered one-leg kernel has its first nonzero term at the third hop.
    targets = []
    for odd in (False, True):
        for q in (1, -1):
            h = precision.generator(board, route, q, odd)
            for i, j in ((0, 3), (3, 0)):
                launch = np.zeros((h.shape[0], 4), dtype=np.int64)
                launch[4 * route[i] : 4 * route[i] + 4] = np.eye(4, dtype=int)
                assert int(abs(h).sum(axis=1).max()) ** 6 < 2**63
                for _ in range(3):
                    assert not np.any(launch[4 * route[j] : 4 * route[j] + 4])
                    launch = h @ launch
                block = launch[4 * route[j] : 4 * route[j] + 4]
                gram = block.T @ block
                index = int(np.argmax(gram.diagonal()))
                assert index == (3 if (q == 1) == (i == 0) else 0)
                expected = np.zeros((4, 4), dtype=int)
                expected[index, index] = 5184 if odd else 4096
                np.testing.assert_array_equal(gram, expected)
                targets.append(
                    dict(
                        odd=odd,
                        charge=q,
                        sender=route[i],
                        receiver=route[j],
                        index=index,
                        coefficient=str(sp.Rational(int(gram[index, index]), 2304)),
                    )
                )
    return targets


def certify():
    board = pairing.formation.federation(1)
    route = pairing.formation.region(board)
    stations = np.array([4 * p + k for p in route for k in range(4)])
    co = PRIOR.weights()
    kernels = {}
    for odd in (False, True):
        for q in (1, -1):
            if not odd and q == 1:
                u, vector, error = PRIOR.certified_kernels(board, route, q, co)
                assert error < 1e-20
                a, b = (sp.Matrix(v.reshape(-1, 4)) for v in vector)
                rho = a.T * a + b.T * b + sp.I * (b.T * a - a.T * b)
                rho /= rho.trace()
            else:
                u = precision.kernel_panels(
                    precision.generator(board, route, q, odd), stations, co
                )
            kernels[odd, q] = u
    sources = bounded_sources(rho)
    covariance = sp.re(sources["formed_pairing"][0][0][0, 0])
    assert covariance - sp.Rational("1e-18") > sp.Rational(17, 1000)
    answer = []
    for name, weights, orientation, cap in (
        ("state", ABSOLUTE, 1, 12),
        ("unpair", CHECKER, 1, 12),
        ("dephased_unpair", CHECKER, 1, 12),
        ("unpair", CHECKER, 1, 8),
        ("unpair", CHECKER, 1, 4),
        ("formed_pairing", CHECKER, -1, 12),
    ):
        assert not np.any(SEND.T @ weights) and not np.any(RECEIVE.T @ weights)
        center, trace_error = sources[name]
        bound = effect_error(cap)
        rows = []
        for q in (1, -1):
            for i in range(4):
                indices = [k for k, (a, b) in enumerate(PAIRS) if a == i and weights[k]]
                if not indices:
                    continue
                val = 0
                for k in indices:
                    a, b = PAIRS[k]
                    for odd in (False, True):
                        op = clock_effect(kernels[odd, q], a, b, cap)
                        val += weights[k] * sp.re(sp.trace(op * center[odd]))
                signed = q * orientation * val
                # A physical state has norm 1, a difference of states has norm <=2.
                # Each one-way capped arrival observable lies between 0 and cap/2.
                uncertainty = sum(abs(weights[k]) for k in indices) * (
                    sp.Rational(cap, 2) * trace_error + 3 * bound
                )
                assert signed > uncertainty
                rows.append(
                    dict(
                        charge=q,
                        sender=route[i],
                        center=float(val),
                        error_upper=float(uncertainty),
                        signed_lower=float(signed - uncertainty),
                    )
                )
        answer.append(
            dict(
                source=name,
                horizon=cap / 2,
                weights=list(weights),
                orientation=orientation,
                rows=rows,
                minimum_signed_lower=min(row["signed_lower"] for row in rows),
            )
        )
    delta, delta_err = sources["unpair"]
    global_delta = sum(delta, sp.zeros(4))
    assert all(c == sp.zeros(2) for c in marginals(global_delta))
    assert global_delta[0, 0] == global_delta[3, 3]
    standardized = sp.Rational(145, 72) * sp.re(global_delta[0, 0])
    uncertainty = sp.Rational(145, 72) * delta_err
    assert standardized + uncertainty < -sp.Rational(232, 1000)
    return dict(
        preprocessed_leading_unpair_slope=[
            float(standardized - uncertainty),
            float(standardized + uncertainty),
        ],
        state_radius_HS="2e-7",
        matched_pairing=certified_geometry(
            kernels, sources, clock_effect, effect_error(12)
        ),
        readings=answer,
        shortest_head_projectors=leading_pair_transfer(board, route),
        formed_leading_common_population_slope=float(sp.Rational(16, 9) * covariance),
    )


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--write", action="store_true")
    opt = cli.parse_args()
    text = json.dumps(certify(), indent=2, allow_nan=False) + "\n"
    if opt.write:
        (HERE / "certificate.json").write_text(text)
    else:
        print(text, end="")
