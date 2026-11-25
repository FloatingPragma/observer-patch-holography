"""Certified retained modular responses on the common coupled eight-state region.

K is minus the natural logarithm of the reference joint density after the
four-tick pump. Every slope uses that same K. A rotation is a unit angular
displacement about the spectral axis of one incoming head marginal; both
marginals are invariant. Unitary conjugation leaves incoming entropy fixed,
while the mobile trace can give a signed modular response on the retained region.
"""

import argparse
import json

import numpy as np
import sympy as sp

from . import live_source, live_survey
from .source import certificate as cert
from .source import measure as source


def entry_bound(mat):
    """Entrywise nuclear-norm upper bound valid also for non-Hermitian readback."""
    return sum(abs(z) for side in mat.as_real_imag() for z in side)


def rounded(a):
    """Reduce denominators, with a computed trace-norm radius."""
    unit = 2**78
    re, im = a.as_real_imag()
    m = (
        re.applyfunc(lambda x: round(x * unit))
        + sp.I * im.applyfunc(lambda x: round(x * unit))
    ) / unit
    assert m == m.H
    return m, entry_bound(a - m)


def spectral_center(matrix):
    """Validate a floating orthonormal eigensystem proposal by polar decomposition."""
    values, vectors = np.linalg.eigh(np.array(matrix, dtype=complex))
    v = cert.prior.sym_matrix(vectors)
    w = sp.diag(*[sp.Rational(float(k)) for k in values])
    defect = cert.prior.one_norm(v.H * v - sp.eye(matrix.rows))
    assert defect < sp.Rational(1, 1000)
    residue = cert.prior.one_norm(matrix - v * w * v.H)
    radius = residue + 3 * defect * max(abs(k) for k in w.diagonal())
    return v, w, defect, radius


def logarithm(ref, error):
    v, w, defect, radius = spectral_center(ref)
    radius += error
    lower = sp.Rational(14, 10000)
    assert min(w.diagonal()) - radius > lower
    assert max(w.diagonal()) < 1
    logs = [sp.Rational(float(np.log(float(a)))) for a in w.diagonal()]
    for a, b in zip(w.diagonal(), logs):
        cert.prior.log_bound(a, b)
    log_approx = (-v * sp.diag(*logs) * v.H).applyfunc(sp.expand)
    # In the unitary polar basis, the resolvent formula bounds the matrix log
    # by radius/lower. Returning to V costs <= 3 ||log(D)|| defect.
    log_error = (
        radius / lower + sp.Rational(1, 10**12) + 3 * defect * max(abs(b) for b in logs)
    )
    numerical = np.array(log_approx, dtype=complex)
    num_round = cert.prior.one_norm(cert.prior.sym_matrix(numerical) - log_approx)
    return log_approx, log_error, numerical, log_error + num_round, lower


def rotation_centers(formed, ferr):
    """Head rotations use the actual marginal eigenvalue gap as angular unit."""
    left = sp.Matrix(
        2, 2, lambda i, j: sum(formed[2 * i + k, 2 * j + k] for k in range(2))
    )
    right = sp.Matrix(
        2, 2, lambda i, j: sum(formed[2 * k + i, 2 * k + j] for k in range(2))
    )
    result = {}
    for tag, m in (("left", left), ("right", right)):
        m = m - sp.eye(2) / 2
        gap2 = sp.re(2 * sp.trace(m * m))
        assert sp.Rational(1, 4) + 4 * ferr < gap2 < 1
        nu = sp.Rational(float(1 / np.sqrt(float(gap2))))
        # True and center gaps exceed 1/2 on this source. Weyl bounds their
        # difference by 2*ferr. Check the numerical reciprocal independently.
        reciprocal_error = 2 * abs(nu**2 * gap2 - 1)
        generator_error = 10 * ferr + reciprocal_error
        g = (
            sp.kronecker_product(nu * m, sp.eye(2))
            if tag == "left"
            else sp.kronecker_product(sp.eye(2), nu * m)
        )
        # The true generator has norm 1/2, so the input commutator error is
        # <= 2 ||G-G0|| ||rho_center||_1 + ||rho-rho_center||_1.
        t = (-sp.I * (g * formed - formed * g)).applyfunc(sp.expand)
        result[tag] = (t, 2 * generator_error + ferr)
    return result


def certify():
    (formed, product), ferr, _ = live_source.formation_center()
    _, _, colerr, grams = live_source.pump_center(True)
    native_error = 4 * colerr + 4 * colerr**2

    def forward(head, weights):
        return live_source.rational_read(grams, head, weights)

    ref = forward(formed, (2, 1))
    ref_center, ref_round = rounded(ref)
    log, logerr, numlog, numlogerr, cutoff = logarithm(
        ref_center, ferr + native_error + ref_round
    )
    center_diff = {
        "ordinary": (forward(formed, (1, 2)) - ref, 2 * (ferr + native_error)),
        "decorrelated": (forward(product, (2, 1)) - ref, 3 * ferr + 2 * native_error),
    }
    axes = live_source.qubit_hidden.correlation_axes()
    _, _, _, refined = live_source.exact_witness(grams, colerr, formed, ferr)
    tangent = sum(
        (
            sp.Rational(v) * axes[j]
            for v, j in zip(refined["coefficients"], (8, 2, 7, 6))
        ),
        sp.zeros(4),
    )
    input_error = sp.Rational(refined["l2_error"])
    center_diff["hidden"] = (
        forward(tangent, (2, 1)),
        input_error + (1 + input_error) * native_error,
    )
    for tag, (inp, err) in rotation_centers(formed, ferr).items():
        center_diff[tag] = forward(inp, (2, 1)), err + (1 + err) * native_error
    _, numerical, charges = live_survey.retained(True)
    slopes, output_errors, rotation_norms = {}, {}, {}
    for tag, (out, err) in center_diff.items():
        c, r = rounded(out)
        err += r
        check = err + entry_bound(c - cert.prior.sym_matrix(numerical[tag]))
        output_errors[tag] = cert.enclosed(check, 0)[1]
        val = sp.re(sp.trace(log * c))
        # Weyl bounds each of the eight eigenvalues. Use the much smaller
        # measured trace norm of the retained tangent in the log error budget.
        _, w, _, radius = spectral_center(c)
        trace_center = sum(abs(v) for v in w.diagonal())
        trace_upper = trace_center + 8 * radius
        # lambda_min > .0014 and lambda_max <=1 give ||K|| <7.
        slope_err = trace_upper * logerr + 7 * err
        slopes[tag] = cert.enclosed(val, slope_err)
        if tag in ("left", "right"):
            rotation_norms[tag] = cert.enclosed(trace_center, 8 * radius + err)
    return dict(
        cells=80,
        forming_charge=-1,
        pump_charge=-1,
        pump_ticks=4,
        region_order=["bridge strain", "original left head", "original right head"],
        reference_minimum_eigenvalue_lower=float(cutoff),
        reference_log_real=numlog.real.tolist(),
        reference_log_imag=numlog.imag.tolist(),
        reference_log_op_error=cert.enclosed(numlogerr, 0)[1],
        modular_slopes=slopes,
        survey_modular_centers=charges["retained"],
        rotation_output_trace_norms=rotation_norms,
        rotation_trace_norm_errors={
            tag: output_errors[tag] for tag in ("left", "right")
        },
        nominal_tangent_trace_norm_errors=output_errors,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    text = json.dumps(certify(), indent=2, allow_nan=False) + "\n"
    if parser.parse_args().write:
        (source.HERE.parent / "live_modular.json").write_text(text)
    else:
        print(text, end="")
