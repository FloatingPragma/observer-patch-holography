#!/usr/bin/env python3
"""Exact local-tick null and inherited N4 strictness certificate.

This never constructs or estimates the full N4 second local-tick moment.
Strict positivity follows from the proof in README.md and an independently
recomputed existing edge Gram certificate.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
import sys

import sympy as sp
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT / "code"))
from native_response.certificate import certify_fixed_angles
from native_response.source import edges


def verify_sources():
    pins = json.loads((HERE / "source_pins.json").read_text())
    for name, digest in pins.items():
        actual = hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
        if actual != digest:
            raise ValueError(f"source digest mismatch: {name}")
    return pins


def baseline(alpha):
    """Exact closed-mean comparator, per observer-incident tick."""
    alpha = F(alpha)
    if not 0 <= alpha <= 1:
        raise ValueError("closed comparator requires 0 <= alpha <= 1")
    c1 = 1 - 6 * alpha / 11
    mixing = 6 * alpha**2 * (19 * alpha + 30) / (
        121 * (11 * alpha**2 + 50 * alpha + 50))
    return c1, c1**2 + mixing


def corrected(c1, c2):
    """Exact rational statistic; no statistical error rule is implied."""
    c1, c2 = F(c1), F(c2)
    r = 1 - c1
    return c2 - c1**2 - r**2 * (209 * r + 180) / (
        1331 * r**2 + 3300 * r + 1800)


def mean_maps(alpha, observer):
    """Independent full 12-port compression, in column-mean convention."""
    alpha = sp.Rational(alpha)
    a, b = sp.zeros(12), sp.zeros(12)
    for left, right in edges():
        v = sp.zeros(12, 1)
        v[left], v[right] = 1, -1
        seam = sp.eye(12) - alpha * v * v.T / 2
        if observer in (left, right):
            b += seam / 30
        else:
            a += seam / 30
    return a, b, b * (sp.eye(12) - a).inv()


def mean_correlations(alpha, observer):
    _, _, kernel = mean_maps(alpha, observer)
    f = sp.ones(12, 1) * -sp.Rational(1, 12)
    f[observer] += 1
    norm = (f.T * f)[0]
    return tuple(F((f.T * kernel**n * f)[0] / norm) for n in (1, 2))


def symbolic_shell_identity():
    """Check the shell inverse before substituting either coupling."""
    a = sp.Symbol("alpha")
    remote = sp.Matrix([[0, 0, 0, 0], [0, 2, -2, 0],
                        [0, -2, 3, -1], [0, 0, -5, 5]])
    column = sp.Matrix([5, -1, 0, 0])
    value = (sp.Matrix([[5, -5, 0, 0]]) *
             (10 * sp.eye(4) + a * remote).inv() * column)[0]
    expected = 5 * (2*a + 3) * (3*a + 10) / (11*a*a + 50*a + 50)
    assert sp.cancel(value - expected) == 0
    c2 = 1 - sp.Rational(12, 11)*a + sp.Rational(6, 55)*a*a*value
    c1 = 1 - sp.Rational(6, 11)*a
    excess = 6*a*a*(19*a + 30) / (121*(11*a*a + 50*a + 50))
    assert sp.cancel(c2 - c1*c1 - excess) == 0
    r = sp.Symbol("r")
    assert sp.cancel(excess.subs(a, 11*r/6) -
                     r*r*(209*r + 180)/(1331*r*r + 3300*r + 1800)) == 0
    return True


def calculate():
    pins = verify_sources()
    symbolic_shell_identity()
    exact_rows = []
    for alpha in (F(0), F(1, 9), F(1, 2), F(1)):
        expected = baseline(alpha)
        for observer in range(12):
            assert mean_correlations(alpha, observer) == expected
        assert corrected(*expected) == 0
        exact_rows.append({"alpha": str(alpha), "c1": str(expected[0]),
                           "c2": str(expected[1]), "corrected": "0"})
    rows = []
    for certificate in certify_fixed_angles():
        interval = certificate["enclosures"]["gershgorin_margin"]
        m = F(interval["lower_numerator"], interval["denominator"])
        lower = F(3, 50) * m
        assert lower > 0
        rows.append({"coupling": f"1/{certificate['theta_denominator']}",
                     "certified_b_tick_lower_bound": str(lower),
                     "edge_gram_certificate": certificate})
    return {"schema": 1, "date": "2026-10-12", "source_sha256": bound_pins(pins),
            "tick_probability": "1/6", "occupation_variance": "5/36",
            "exact_comparator_checks": {"observers": 12, "rows": exact_rows},
            "symbolic_shell_identity": True, "n4_certified_lower_bounds": rows,
            "scope": "One/two incident-event endpoint instrument; full N4 c2 not numerically evaluated; preparation, detection and readout supplied."}


def main():
    parser = argparse.ArgumentParser()
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--write", type=Path)
    group.add_argument("--check", type=Path)
    args = parser.parse_args()
    result = json.dumps(calculate(), indent=2, sort_keys=True) + "\n"
    if args.check:
        if args.check.read_text() != result:
            raise SystemExit("Local-tick receipt mismatch")
        print("Local-tick receipt reproduced")
    elif args.write:
        with args.write.open("x") as stream:
            stream.write(result)
        print(args.write)
    else:
        print(result, end="")


if __name__ == "__main__":
    main()
