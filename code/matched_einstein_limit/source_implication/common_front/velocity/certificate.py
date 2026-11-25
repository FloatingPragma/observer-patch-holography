"""Uniform likelihood and source-derivative bounds for every positive head readout."""

import argparse
from itertools import product
import json

import numpy as np
import sympy as sp
from mpmath import iv
from scipy.linalg import eigvalsh

from .. import certificate as front
from .measure import HERE


def positive(a):
    return all(a[:k, :k].det() > 0 for k in range(1, a.rows + 1))


def run():
    iv.dps = 32
    mid, targets, axis, rotations, r = front.inputs()
    eye = sp.eye(4)
    assert all(positive(x - sp.Rational(1, 50) * eye) for x in mid)
    m = sp.Rational(19999, 10**6)
    assert sp.Rational(1, 50) - r > m
    angles = json.loads(
        (HERE.parent.parent / "null_optics/certificate.json").read_text()
    )["matched_pairing"]
    scores = {}
    for name, field in [
        ("ordinary", "ordinary_angle_interval"),
        ("unpair", "angle_interval"),
    ]:
        a, b = angles[field]
        alpha = iv.mpf([int(np.floor(a)), int(np.ceil(b))])
        # Outward integer corners contain every admissible phase adjustment
        # times normalized Bloch-axis component. The entire rectangular hull
        # is checked. Correlations of the three coefficients can only shrink it.
        coeff = []
        for x in axis:
            bounds = alpha * x
            low = int(np.floor(float(bounds.a))) - 1
            high = int(np.ceil(float(bounds.b))) + 1
            assert low < bounds.a <= bounds.b < high
            coeff.append([low, high])
        eps = sp.Rational("6e-7") * (2 if name == "ordinary" else 3)
        eps += r * sum(max(abs(lo), abs(hi)) for lo, hi in coeff)
        corners = []
        minimum = 0.0
        for beta in product(*coeff):
            for b in range(2):
                x = mid[b] - r * eye
                d = targets[name][b] + sum(
                    (c * rs[b] for c, rs in zip(beta, rotations)), sp.zeros(4)
                )
                for sign in (-1, 1):
                    h = sign * d + eps * eye
                    # Floating generalized eigenvalues only propose an integer.
                    minimum = max(
                        minimum,
                        float(
                            eigvalsh(
                                np.array(h, dtype=complex), np.array(x, dtype=complex)
                            )[-1]
                        ),
                    )
                    corners.append((x, h))
        factor = int(np.ceil(minimum)) + 1
        assert all(positive(factor * x - h) for x, h in corners)
        scores[name] = dict(
            absolute_log_score_upper=factor,
            rotation_product_intervals=coeff,
            rational_matrix_inequalities_checked=len(corners),
        )
    return dict(
        reference_minimum_eigenvalue_lower=str(m),
        finite_ordinary_or_pairing_likelihood_lower=str(8 * m * m),
        finite_ordinary_or_pairing_likelihood_upper=str(1 / m),
        score=scores,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    out = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if p.parse_args().write:
        (HERE / "certificate.json").write_text(out)
    else:
        print(out, end="")
