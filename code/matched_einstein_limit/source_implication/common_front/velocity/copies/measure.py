"""Replayed two-copy rail detector, complete binary outcomes, and collective lower bounds."""

import argparse
from fractions import Fraction
import json
from math import ceil
from pathlib import Path

import numpy as np
import sympy as sp
from mpmath import iv

from .. import certificate as previous
from .... import signed_transfers

HERE = Path(__file__).parent
front = previous.front
eye = sp.eye(2)
X = sp.Matrix([[0, 1], [1, 0]])
Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
Z = sp.diag(1, -1)
PAULIS = [X, Y, Z]


def rail_moments(mid, r, all_parities=False):
    def lim(center, error):
        return front.full.ri_band(sp.re(center), error)

    def inner(h, x):
        return sp.re((h * x).trace())

    def moment(a):
        value = lim(inner(sp.kronecker_product(Z, a), mid[0]), 2 * r)
        if all_parities:
            # At native dwell pi/8 the odd marked contact has frequency 9.
            # Its left read seen from the coin is diag(1,cos(9*pi/8)).
            c = iv.sqrt(2 + iv.sqrt(iv.mpf(2))) / 2
            odd_i = lim(inner(sp.kronecker_product(eye, a), mid[1]), 2 * r)
            odd_z = lim(inner(sp.kronecker_product(Z, a), mid[1]), 2 * r)
            value += ((1 - c) * odd_i + (1 + c) * odd_z) / 2
        return value

    g = [moment(a) for a in PAULIS]
    left = moment(eye)
    v = [
        lim(sum(inner(sp.kronecker_product(eye, a), x) for x in mid), 3 * r)
        for a in PAULIS
    ]
    return g, [left * t for t in v]


def source_data():
    mid, targets, _, _, r = front.inputs()
    g, prod = rail_moments(mid, r)
    odds = front.full.ri_band(mid[1].trace(), 2 * r)
    m = sp.Rational("0.019999")
    energy = {}
    for label, name, err in [
        ("ordinary", "ordinary", sp.Rational("12e-7")),
        ("collar", "unpair", sp.Rational("18e-7")),
    ]:
        d = targets[name]
        q0 = sum((a * (x - r * sp.eye(4)).inv() * a).trace() for a, x in zip(d, mid))
        # Triangle bound in the Hilbert-Schmidt norm after right multiplication
        # by the inverse square root of the perturbed reference density.
        q = (
            iv.sqrt(front.full.ri(sp.re(q0)))
            + front.full.ri(err) / iv.sqrt(front.full.ri(m))
        ) ** 2
        hs = iv.sqrt(
            front.full.ri(sum(sp.re((a * a).trace()) for a in d))
        ) + front.full.ri(err)
        energy[label] = (q, hs)
    return mid, targets, g, prod, odds, energy


def power(p0, delta, p0_nom, changed_nom, rng):
    mag = abs(delta).a
    scale = 10**18
    lower = Fraction(max(0, int(float(mag) * scale) - 1024), scale)
    exact = front.full.ri(sp.Rational(lower.numerator, lower.denominator))
    assert 0 < exact <= mag
    n = ceil(Fraction(19) / (lower * lower))
    # Two independent cohorts, a labelled reference and the blinded test.
    # Each product statistic is binary, so retain click and complement.
    draws = {}
    for label, p in [
        ("reference", p0_nom),
        ("unchanged", p0_nom),
        ("changed", changed_nom),
    ]:
        c = int(rng.binomial(n, p))
        draws[label] = [c, n - c]
    sign = 1 if delta.a > 0 else -1
    differences = {
        k: sign * (v[0] - draws["reference"][0]) / n
        for k, v in draws.items()
        if k != "reference"
    }
    assert differences["unchanged"] < float(lower) / 2 < differences["changed"]
    return dict(
        reference_probability=str(p0),
        probability_change=str(delta),
        minimum_gap_used=str(lower),
        samples_per_cohort=n,
        counts=draws,
        signed_differences=differences,
        decision_threshold=float(lower) / 2,
    )


def needed(q):
    eps = iv.mpf("0.01")
    requirement = iv.ln(1 + (1 - 2 * eps) ** 2 / (eps * (1 - eps))) / iv.ln(1 + q)
    lower = ceil(float(requirement.a))
    # Correct a floating proposal using the outward interval endpoint.
    while iv.mpf(lower - 1) >= requirement.a:
        lower -= 1
    while iv.mpf(lower) < requirement.a:
        lower += 1
    assert iv.mpf(lower - 1) < requirement.a <= iv.mpf(lower)
    return dict(
        single_copy_chi_square_upper=str(q.b),
        tested_copies_required_at_least=lower,
        real_lower_bound_interval=str(requirement),
    )


def run():
    iv.dps = 40
    mid, targets, g, h, p0, energy = source_data()
    nominal = np.array([np.array(x, dtype=complex) for x in mid])
    raw = np.array([np.array(x, dtype=complex) for x in targets["unpair"]])

    # Exact two-copy product statistic from the full parity+rail instrument.
    def correlations(x):
        z = np.array(Z, dtype=complex)
        corr = [
            np.trace(x[0] @ np.kron(z, np.array(a, dtype=complex))).real for a in PAULIS
        ]
        return np.array(corr)

    # Drive the large seeded binomials from the rational source center,
    # converting each final probability to binary64 exactly once.
    def exact_correlations(block):
        return sp.Matrix([sp.re((block * sp.kronecker_product(Z, a)).trace()) for a in PAULIS])

    g_nom = exact_correlations(mid[0])
    d_nom = exact_correlations(targets["unpair"][0])
    odds = sp.re(mid[1].trace())
    product_base = sp.S.Half + g_nom.dot(g_nom) / 6

    def probability(x):
        v = correlations(x)
        return 0.5 + v @ v / 6

    pg = 0.5 + sum(a * a for a in g) / 6
    dc = [t - s for t, s in zip(h, g)]
    derivative = sum(a * b for a, b in zip(g, dc)) / 3
    quadratic = sum(a * a for a in dc) / 6
    shots = np.random.default_rng(270034)
    saved = json.loads(
        (Path(signed_transfers.__file__).with_name("signed_transfers.json")).read_text()
    )
    curves = {row["source"]: row for row in saved["tests"]}
    panels = {}
    for t in [sp.S.One, sp.Rational(1, 8), sp.Rational(1, 5000000)]:
        s = front.full.ri(t)
        tf = float(t)
        panels[str(t)] = {}
        for name in ["ordinary", "collar"]:
            prob0, delta = (
                (p0, s * (0.5 - p0))
                if name == "ordinary"
                else (pg, s * (derivative + s * quadratic))
            )
            base, dnom = (
                (odds, t * (sp.S.Half - odds))
                if name == "ordinary"
                else (product_base, t * g_nom.dot(d_nom) / 3 + t**2 * d_nom.dot(d_nom) / 6)
            )
            if name == "collar":
                np.testing.assert_allclose(
                    probability(nominal + tf * raw) - float(base), float(dnom),
                    atol=2e-12, rtol=0,
                )
            assert float(base) in prob0 and float(dnom) in delta
            row = power(prob0, delta, float(base), float(base + dnom), shots)
            row["tested_heads_per_bernoulli"] = 1 if name == "ordinary" else 2
            chi2, hs = energy[name]
            row["arbitrary_joint_read_raw_lottery"] = needed(s * s * chi2)
            if t == sp.Rational(1, 5000000):
                a = curves[name]["finite"][0]["angle"]
                assert curves[name]["finite"][0]["probability"] == tf
                # The installed endpoint angle is fixed. The signed meter is
                # invariant under it; the optimal joint distinguishability can rise.
                angle = iv.mpf(str(abs(a))) + iv.mpf("1e-17")
                q = (s * hs + angle) ** 2 / iv.mpf("0.019999")
                row["archived_adjusted_angle"] = a
                row["arbitrary_joint_read_installed_angle"] = needed(q)
            panels[str(t)][name] = row
    g_all, h_all = rail_moments(mid, front.inputs()[-1], all_parities=True)
    d_all = [b - a for a, b in zip(g_all, h_all)]
    a_all = sum(a * d for a, d in zip(g_all, d_all)) / 3
    b_all = sum(d * d for d in d_all) / 6
    assert a_all.b < 0 and (a_all + b_all).b < 0
    return dict(
        virtual_count_seed=270034,
        unconditional_native_rail=dict(
            reference_probability=str(0.5 + sum(a * a for a in g_all) / 6),
            collar_first_derivative=str(a_all),
            collar_quadratic_coefficient=str(b_all),
            archived_endpoint_change=str((a_all + b_all / 5000000) / 5000000),
        ),
        coin_Z_right_random_axis_product_derivative=str(derivative),
        coin_product_quadratic=str(quadratic),
        parity_derivative=str(0.5 - p0),
        collar_parity_derivative=0,
        confidence_error_per_hypothesis=0.01,
        positive_source_resource_history=panels,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--write", action="store_true")
    data = json.dumps(run(), indent=2, allow_nan=False) + "\n"
    if p.parse_args().write:
        (HERE / "results.json").write_text(data)
    else:
        print(data, end="")
