"""A commuting source family at fixed volume, energy and both marginal entropies."""

import argparse
import json
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp

from ..common_front.certificate import inputs
from ..full_algebra import certificate as full
from ..full_algebra import native_precision as native
from ..null_optics.certificate import marginals
from .controls import nominal_experiment

HERE = Path(__file__).resolve().parent
iv = mp.iv


def det3(columns):
    a, b, c = columns
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - b[0] * (a[1] * c[2] - a[2] * c[1])
        + c[0] * (a[1] * b[2] - a[2] * b[1])
    )


def real_bound(mid, error):
    return full.ri(mid) + iv.mpf([-1, 1]) * error


def verify():
    iv.dps = 32
    x, _, _, _, radius = inputs()
    ir = full.ri(radius)
    eye = sp.eye(4)
    dirs = []
    errors = []
    for sector, (offset, scale) in enumerate(
        (
            (sp.Rational(1, 5), sp.Rational(1, 10)),
            (sp.Rational(1, 20), sp.Rational(3, 100)),
        )
    ):
        h = (x[sector] - offset * eye) / scale
        for positive in (eye + h, eye - h):
            assert all(positive[:k, :k].det() > 0 for k in range(1, 5))
        h_error = ir / full.ri(scale)
        for power, err in ((h, h_error), (h * h, (2 + h_error) * h_error)):
            dirs.append((power - eye * power.trace() / 4).applyfunc(sp.expand))
            errors.append(err)
    nums, unit, ge = native.exact_effects()
    effects = full.rational_effects(nums, unit)
    eps = full.ri(sp.Rational(ge))
    norms = [iv.sqrt(full.ri(full.inner(d, d))) for d in dirs]
    ranges = []
    response = []
    for e in effects:
        trace = sum(t.trace() for t in e) / 8
        norm = iv.sqrt(
            full.ri(sum(full.inner(t - trace * eye, t - trace * eye) for t in e))
        )
        ranges.append(
            real_bound(sum(full.inner(m, t) for m, t in zip(x, e)), eps + norm * ir)
        )
        row = []
        for j, d in enumerate(dirs):
            b = e[j // 2]
            centered = b - eye * b.trace() / 4
            width = iv.sqrt(full.ri(full.inner(centered, centered)))
            row.append(
                real_bound(
                    full.inner(b, d),
                    width * errors[j] + 2 * eps * (norms[j] + errors[j]),
                )
            )
        response.append(row)
    area, volume = full.enc.geometry_jet(ranges)
    response = np.array(response, dtype=object)
    a, v = area.g @ response, volume.g @ response
    left, right = [marginals(b)[0] for b in x], marginals(x[0] + x[1])[1]
    logs = [
        [sp.kronecker_product(b, sp.eye(2)) for b in left],
        [sp.kronecker_product(sp.eye(2), right)] * 2,
        x,
    ]
    bands = []
    for blocks in logs:
        values = []
        for j, d in enumerate(dirs):
            b = blocks[j // 2]
            # The log varies by at most ||delta||/.019. For a marginal the
            # operator error is at most the full trace-norm error, <3*radius.
            slopes, log_error = native.prior.entropy_coordinates(
                b, 3 * radius, sp.Rational(".019"), [d]
            )
            # d and its displacement are traceless: subtracting any scalar
            # from the near-reference log does not change the expectation.
            # Its centered norm is bounded by half its spectral spread.
            w = np.linalg.eigvalsh(np.array(b, dtype=complex))
            lo, hi = (
                sp.Rational(int(w[0] * 10000) - 2, 10000),
                sp.Rational(int(w[-1] * 10000) + 2, 10000),
            )
            for positive in (b - (lo + 3 * radius) * eye, (hi - 3 * radius) * eye - b):
                assert all(positive[:k, :k].det() > 0 for k in range(1, 5))
            spread = iv.ln(full.ri(hi / lo)) / 2
            values.append(
                real_bound(
                    slopes[0],
                    2 * norms[j] * full.ri(log_error) + 2 * spread * errors[j],
                )
            )
        bands.append(values)

    cols = list(zip(v, bands[0], bands[1]))
    # Fix the coefficient of the even quadratic shape to one. The other
    # three coefficients match volume and the two individual entropies.
    chosen = [cols[i] for i in (0, 2, 3)]
    det = det3(chosen)
    assert det.a > 0 or det.b < 0
    coeff = [None, iv.mpf(1), None, None]
    for where, j in enumerate((0, 2, 3)):
        changed = list(chosen)
        changed[where] = [-val for val in cols[1]]
        coeff[j] = det3(changed) / det

    # Evaluate the small held signal after subtracting constraint rows to
    # reduce interval dependency from cancellation between large columns.
    def constrained(row):
        mid = np.array([[float(r.mid) for r in c] for c in chosen])
        dual = np.linalg.solve(mid, np.array([float(row[k].mid) for k in (0, 2, 3)]))
        weights = [iv.mpf(float(t)) for t in dual]
        return sum(
            y * (r - sum(w * z for w, z in zip(weights, col)))
            for y, r, col in zip(coeff, row, cols)
        )

    ds, da = constrained(bands[2]), constrained(a)
    z = sp.kronecker_product(sp.eye(2), sp.diag(1, -1))
    rz = constrained(
        [real_bound(full.inner(z, d), 2 * er) for d, er in zip(dirs, errors)]
    )
    # Right marginal has trace one, and its squared eigenvalue gap is the
    # sum of squares of the three Pauli expectations.
    pauli = [
        sp.Matrix([[0, 1], [1, 0]]),
        sp.Matrix([[0, -sp.I], [sp.I, 0]]),
        sp.diag(1, -1),
    ]
    gap = iv.sqrt(sum(real_bound(full.inner(right, p), 3 * ir) ** 2 for p in pauli))
    frequency = iv.ln((1 + gap) / (1 - gap))
    assert ds.a > 0 and da.b < 0 and rz.a > 0
    s = iv.ln(2) / (2 * iv.pi)
    # cosh(log 2)=5/4, evaluated as an exact rational.
    two_sided_defect = full.ri(sp.Rational(5, 4)) - iv.cos(frequency * s)
    cosine = iv.cos(frequency * s)
    forward_defect = (5 - 4 * cosine) / (1 + 2 * cosine)
    return dict(
        true_entropy_slope=str(ds),
        true_area_slope=str(da),
        right_Z_slope=str(rz),
        coefficients=[str(c) for c in coeff],
        minor=str(det),
        right_modular_frequency=str(frequency),
        right_doubling_error_per_unit_response=str(two_sided_defect),
        right_forward_error_per_unit_response=str(forward_defect),
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    data = verify()
    data["rounded_source_control"] = nominal_experiment()
    result = json.dumps(data, indent=2, allow_nan=False) + "\n"
    if args.write:
        (HERE / "results.json").write_text(result)
    else:
        print(result, end="")
