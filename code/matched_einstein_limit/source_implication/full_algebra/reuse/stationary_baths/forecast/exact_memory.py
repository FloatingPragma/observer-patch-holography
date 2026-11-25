"""Certify equal native initial ranges and a different response to four radar cycles."""

import argparse
import json

import numpy as np
import sympy as sp

from .. import measure as bath
from ... import invisible
from . import measure

R = sp.Rational
iv = measure.certificate.c.iv
ri = measure.certificate.c.ri
SCALE = 2**46
P = [
    sp.Matrix(a.real.astype(int))
    if not np.any(a.imag)
    else sp.Matrix(a.real.astype(int)) + sp.I * sp.Matrix(a.imag.astype(int))
    for a in measure.PAULI
]


def dyadic(a, scale=SCALE):
    a = np.asarray(a)
    real = sp.Matrix(np.rint(a.real * scale).astype(np.int64))
    imag = sp.Matrix(np.rint(a.imag * scale).astype(np.int64))
    return (real + sp.I * imag) / scale


def fro(m):
    return iv.sqrt(ri(sum(sp.re(x) ** 2 + sp.im(x) ** 2 for x in m)))


def head_coordinates(channel):
    # Hermitian Pauli coordinates satisfy ||x||_2 = 2 ||rho||_HS.
    # Each coefficient sums sixteen signed entries of the binary64 channel.
    basis = np.concatenate(([np.eye(4)], measure.PAULI)).reshape(16, 16)
    return dyadic((basis.conj() @ channel @ basis.T).real / 4)


def study(channels, effects, frozen=None):
    iv.dps = 35
    # Rounding the symmetrized native static effect to 46 fractional bits
    # adds <4.1e-14 in operator norm, including conversion from its integer.
    er = ri(R("1e-13"))
    # Averaging already bounded binary64 channels and rounding entries to 46
    # fractional bits adds <2e-12 in induced Hermitian trace norm (16x16
    # Frobenius bound times two for 4x4 output trace norm).
    ce = ri(invisible.stacked_error() + R("2e-12"))
    gain = (1 + ce) ** 4
    effects = [[dyadic((e + e.conj().T) / 2) for e in row] for row in effects]
    # A real coordinate matrix represents the Hermitian projection of the
    # channel. Its multiplication and fourth power use exact rationals.
    ch = [head_coordinates(c) ** 4 for c in channels.mean(axis=1)]
    observables = [
        sp.Matrix(
            [
                [sp.re((row[b] * p).trace()) / 4 for p in [sp.eye(4), *P]]
                for row in effects
            ]
        )
        for b in range(2)
    ]
    future_blocks = [s * c for s, c in zip(observables, ch, strict=True)]
    b = observables[0][:, 1:].row_join(observables[1][:, 1:])
    future = future_blocks[0][:, 1:].row_join(future_blocks[1][:, 1:])
    approx = np.array(b, float)
    k = (
        dyadic(np.linalg.pinv(approx), 2**28)
        if frozen is None
        else sp.Matrix(frozen["right_inverse_numerators"])
        / int(frozen["right_inverse_scale"])
    )
    eps = iv.sqrt(ri(12 * 30)) * er
    defect = fro(b * k - sp.eye(12)) + eps * fro(k)
    if defect.b >= 1:
        raise ArithmeticError("right inverse unresolved")
    inverse = fro(k) / (1 - defect)
    center = future_blocks[0][:, 0] * R("0.8657") + future_blocks[1][:, 0] * R("0.1343")
    g = bath.station.range_gradient(np.array(center, float).ravel())[1]
    _, _, vh = np.linalg.svd(approx, full_matrices=True)
    kernel = vh[12:].T
    numeric = np.array(future, float)
    targets = [
        *[("route_" + str(i), r) for i, r in enumerate(numeric)],
        ("volume", g @ numeric),
    ]
    answers = {}
    # The future effect on a Hermitian input has additional error from
    # channel propagation and static-effect rounding. Each capped score <=3.
    per_unit = 3 * (gain - 1) + er * gain
    mid_range = [ri(x) + iv.mpf([-1, 1]) * per_unit for x in center]

    def volume_slope(ranges, response):
        det, weights = measure.determinant(ranges)
        if det.a <= 0:
            raise ArithmeticError("future tetrahedron unresolved")
        gradient = [
            sum(float(k) * w for k, w in zip(column, weights, strict=True) if k)
            / (12 * iv.sqrt(det) * r)
            for column, r in zip(measure.POWERS.T, ranges, strict=True)
        ]
        return sum(x * y for x, y in zip(response, gradient, strict=True))

    for t, (name, a) in enumerate(targets):
        v = kernel @ (kernel.T @ a)
        v = (
            dyadic((v / np.linalg.norm(v)).reshape(30, 1))
            if frozen is None
            else sp.Matrix(
                frozen["directions_on_exact_initial_range_fiber"][name][
                    "direction_numerators"
                ]
            )
            / int(frozen["scale"])
        )
        # Move from v to D = v - K(BK)^(-1)Bv for the *exact* static B.
        # It preserves both traces because its 30 Pauli coordinates do.
        correction = inverse * (fro(b * v) + eps * fro(v))
        bound = iv.sqrt(ri(2)) * (fro(v) * per_unit + 3 * correction)
        rr = [ri(x) + iv.mpf([-1, 1]) * bound for x in future * v]
        derivative = rr[t] if name != "volume" else volume_slope(mid_range, rr)
        answers[name] = dict(
            direction_numerators=[int(x * SCALE) for x in v],
            coordinate_correction_norm_upper=float(
                np.nextafter(float(correction.b), np.inf)
            ),
            derivative_interval=[
                float(np.nextafter(float(derivative.a), -np.inf)),
                float(np.nextafter(float(derivative.b), np.inf)),
            ],
        )
        if name == "volume":
            half = ri(R(1, 32))
            tube = [
                x + iv.mpf([-1, 1]) * half * y
                for x, y in zip(mid_range, rr, strict=True)
            ]
            separation = 2 * half * volume_slope(tube, rr)
            block_floor = ri(R("0.1343") / 4) - half * (fro(v) + correction) / 2
            if min(separation.a, block_floor.a) <= 0:
                raise ArithmeticError("exact finite volume contrast unresolved")
            answers[name].update(
                half_displacement=1 / 32,
                minimum_source_eigenvalue_lower=float(
                    np.nextafter(float(block_floor.a), -np.inf)
                ),
                finite_volume_difference_interval=[
                    float(np.nextafter(float(separation.a), -np.inf)),
                    float(np.nextafter(float(separation.b), np.inf)),
                ],
            )
    return dict(
        scale=SCALE,
        right_inverse_scale=2**28,
        right_inverse_numerators=[
            [int(k[i, j] * 2**28) for j in range(k.cols)] for i in range(k.rows)
        ],
        right_inverse_defect_upper=float(np.nextafter(float(defect.b), np.inf)),
        future_range_operator_error_upper=float(
            np.nextafter(float(per_unit.b), np.inf)
        ),
        directions_on_exact_initial_range_fiber=answers,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    path = measure.HERE / "exact_memory.json"
    if args.write:
        out = study(*bath.load_native())
        path.write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    else:
        saved = json.loads(path.read_text())
        out = study(*bath.load_native(), frozen=saved)
        if out != saved:
            raise ValueError(
                "exact memory certificate differs from archived dyadic proposals or bounds"
            )
        print(
            out["directions_on_exact_initial_range_fiber"]["volume"][
                "finite_volume_difference_interval"
            ]
        )
