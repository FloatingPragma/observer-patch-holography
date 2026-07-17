#!/usr/bin/env python3
"""Machine witness for the issue #320 collar Poisson theorem.

The witness enumerates exact opportunity-count laws on declared finite
collar families and checks the exact total-variation distance to the
matching Poisson law against the explicit bounds of the quantitative
Poisson theorem in cosmology/oph_dark_matter_paper.tex:

  independent families:  d_TV <= sum_i p_i^2 + |mu_r - mu|   (Le Cam)
  one-dependent family:  d_TV <= b_1 + b_2                   (Chen-Stein,
                                                              b_3 = 0)

Count laws are computed exactly in Fractions.  Poisson probability masses
are evaluated in Decimal arithmetic at DECIMAL_PRECISION digits, so the
reported total-variation values carry an absolute numerical error far
below NUMERIC_TOLERANCE, itself far below every bound margin checked.
"""

from __future__ import annotations

import argparse
import json
import pathlib
from decimal import Decimal, getcontext
from fractions import Fraction
from typing import Any

DECIMAL_PRECISION = 60
NUMERIC_TOLERANCE = Decimal("1e-40")

# Pixel fixed point consumed by the Z6 presence branch, exact decimal.
P_PIXEL = Fraction("1.630968209403959")
LAMBDA_UNIT = Fraction(1)
LAMBDA_Z6_PRESENCE = 1 - P_PIXEL / 24

getcontext().prec = DECIMAL_PRECISION


def frac_str(x: Fraction) -> str:
    return f"{x.numerator}/{x.denominator}"


def frac_to_decimal(x: Fraction) -> Decimal:
    return Decimal(x.numerator) / Decimal(x.denominator)


def bernoulli_sum_law(ps: list[Fraction]) -> list[Fraction]:
    """Exact law of a sum of independent Bernoulli(p_i) by convolution."""
    law = [Fraction(1)]
    for p in ps:
        q = 1 - p
        nxt = [Fraction(0)] * (len(law) + 1)
        for n, mass in enumerate(law):
            nxt[n] += mass * q
            nxt[n + 1] += mass * p
        law = nxt
    return law


def one_dependent_law(m: int, q: int) -> list[Fraction]:
    """Exact law of N = sum_{i=1..m} X_i with X_i = [U_i=0][U_{i+1}=0].

    U_1..U_{m+1} are i.i.d. uniform on {0..q-1}.  X_i depends only on
    (U_i, U_{i+1}), so X_i is independent of {X_j : |i-j| >= 2} and the
    Chen-Stein mixing term b_3 vanishes for neighborhoods
    B_i = {i-1, i, i+1}.
    """
    p_zero = Fraction(1, q)
    p_nonzero = 1 - p_zero
    # state: (count, U_last == 0) -> probability
    states: dict[tuple[int, bool], Fraction] = {
        (0, True): p_zero,
        (0, False): p_nonzero,
    }
    for _ in range(m):
        nxt: dict[tuple[int, bool], Fraction] = {}
        for (count, last_zero), mass in states.items():
            for new_zero, step in ((True, p_zero), (False, p_nonzero)):
                hit = 1 if (last_zero and new_zero) else 0
                key = (count + hit, new_zero)
                nxt[key] = nxt.get(key, Fraction(0)) + mass * step
        states = nxt
    law = [Fraction(0)] * (m + 1)
    for (count, _), mass in states.items():
        law[count] += mass
    return law


def poisson_pmf(mu: Fraction, n_max: int) -> list[Decimal]:
    mu_dec = frac_to_decimal(mu)
    base = (-mu_dec).exp()
    pmf = []
    term = base
    pmf.append(term)
    for n in range(1, n_max + 1):
        term = term * mu_dec / Decimal(n)
        pmf.append(term)
    return pmf


def total_variation(law: list[Fraction], mu: Fraction) -> Decimal:
    """Exact-to-precision d_TV between a finitely supported law and Poisson(mu)."""
    pois = poisson_pmf(mu, len(law) - 1)
    l1 = Decimal(0)
    pois_mass = Decimal(0)
    for mass, pm in zip(law, pois):
        l1 += abs(frac_to_decimal(mass) - pm)
        pois_mass += pm
    tail = Decimal(1) - pois_mass
    return (l1 + tail) / 2


def independent_family(
    family_id: str,
    branch: str,
    lam: Fraction,
    sqrt_x: Fraction,
    m: int,
    weights: list[int] | None,
    drift: Fraction,
) -> dict[str, Any]:
    """Evaluate one independent Bernoulli family against the Le Cam bound."""
    mu = lam * sqrt_x
    mu_r = mu * (1 + drift)
    if weights is None:
        ps = [mu_r / m] * m
    else:
        total = sum(weights)
        ps = [mu_r * w / total for w in weights]
    assert all(0 < p < 1 for p in ps)
    assert sum(ps) == mu_r
    law = bernoulli_sum_law(ps)
    tv = total_variation(law, mu)
    sum_p_sq = sum(p * p for p in ps)
    delta = abs(mu_r - mu)
    bound = sum_p_sq + delta
    p_max = max(ps)
    coarse = mu_r * p_max + delta
    return {
        "family_id": family_id,
        "model": "independent_bernoulli",
        "branch": branch,
        "inputs": {
            "lambda_collar": frac_str(lam),
            "sqrt_x": frac_str(sqrt_x),
            "m_cells": m,
            "weights": weights if weights is not None else "uniform",
            "mean_drift": frac_str(drift),
        },
        "mu_target": frac_str(mu),
        "mu_r": frac_str(mu_r),
        "p_max": frac_str(p_max),
        "bound_kind": "lecam_sum_p_squared_plus_delta",
        "bound_exact": frac_str(bound),
        "bound_decimal": str(frac_to_decimal(bound)),
        "coarse_bound_exact": frac_str(coarse),
        "exact_tv": str(+tv),
        "tv_le_bound": tv + NUMERIC_TOLERANCE <= frac_to_decimal(bound),
        "bound_le_coarse": bound <= coarse,
        "margin_decimal": str(frac_to_decimal(bound) - tv),
    }


def dependent_family(family_id: str, m: int, q: int) -> dict[str, Any]:
    """Evaluate the one-dependent family against the Chen-Stein b1+b2 bound."""
    p = Fraction(1, q * q)
    mu_r = m * p
    law = one_dependent_law(m, q)
    assert sum(law) == 1
    tv = total_variation(law, mu_r)
    # Neighborhoods B_i = {i-1, i, i+1} intersected with {1..m}.
    b1 = Fraction(0)
    for i in range(1, m + 1):
        nbhd = [j for j in (i - 1, i, i + 1) if 1 <= j <= m]
        b1 += sum(p * p for _ in nbhd)
    # E[X_i X_{i+1}] = P(U_i = U_{i+1} = U_{i+2} = 0) = q^-3.
    pair = Fraction(1, q**3)
    b2 = 2 * (m - 1) * pair
    bound = b1 + b2
    return {
        "family_id": family_id,
        "model": "one_dependent_indicators",
        "branch": "chen_stein_finite_range",
        "inputs": {"m_cells": m, "q_alphabet": q, "dependence_range_cells": 1},
        "mu_r": frac_str(mu_r),
        "b1_exact": frac_str(b1),
        "b2_exact": frac_str(b2),
        "b3_exact": "0/1",
        "bound_kind": "chen_stein_b1_plus_b2",
        "bound_exact": frac_str(bound),
        "bound_decimal": str(frac_to_decimal(bound)),
        "exact_tv": str(+tv),
        "tv_le_bound": tv + NUMERIC_TOLERANCE <= frac_to_decimal(bound),
        "margin_decimal": str(frac_to_decimal(bound) - tv),
    }


def compute() -> dict[str, Any]:
    families: list[dict[str, Any]] = []
    for branch, lam in (
        ("unit", LAMBDA_UNIT),
        ("z6_presence", LAMBDA_Z6_PRESENCE),
    ):
        for sqrt_x in (Fraction(1, 2), Fraction(1), Fraction(2)):
            for m in (8, 64):
                families.append(
                    independent_family(
                        f"independent_uniform_{branch}_sqrtx_{sqrt_x.numerator}_{sqrt_x.denominator}_m{m}",
                        branch,
                        lam,
                        sqrt_x,
                        m,
                        None,
                        Fraction(0),
                    )
                )
    for m in (8, 64):
        families.append(
            independent_family(
                f"independent_nonuniform_z6_presence_m{m}",
                "z6_presence",
                LAMBDA_Z6_PRESENCE,
                Fraction(1),
                m,
                list(range(1, m + 1)),
                Fraction(0),
            )
        )
    for m in (8, 64):
        families.append(
            independent_family(
                f"mean_drift_unit_m{m}",
                "unit",
                LAMBDA_UNIT,
                Fraction(1),
                m,
                None,
                Fraction(1, m),
            )
        )
    for m in (8, 64):
        families.append(dependent_family(f"one_dependent_q6_m{m}", m, 6))

    all_pass = all(f["tv_le_bound"] for f in families)
    return {
        "issue": 320,
        "certificate_id": "issue-320-collar-poisson-witness-v1",
        "claim": (
            "On refinement-stable collar opportunity models, the exact "
            "opportunity-count law sits within the explicit Le Cam or "
            "Chen-Stein total-variation bound of the matching Poisson law "
            "on every declared finite family."
        ),
        "theorem": {
            "statement": (
                "d_TV(L(N_r(J)), Poisson(mu(J))) <= sum_i p_{r,i}^2 + delta_r "
                "<= mu_r(J) p_max(r) + delta_r under disjoint-collar "
                "independence, rare local activation, and refinement "
                "stability; b_1 + b_2 replaces the Le Cam term under "
                "finite-range dependence with b_3 = 0."
            ),
            "premises": [
                "disjoint-collar independence (or finite dependence range R_r with exact neighborhood independence)",
                "rare local activation: p_max(r) -> 0",
                "refinement stability: |mu_r(J) - mu(J)| <= delta_r -> 0",
            ],
            "paper_surface": "cosmology/oph_dark_matter_paper.tex, thm:collar-poisson, cor:collar-poisson-process, eq:lecambound",
            "lambda_collar_identification": (
                "lambda_collar is the refinement-stable active-opportunity "
                "intensity per unit dimensionless cut coordinate sqrt(x); "
                "its numerical value is branch data (unit branch 1, Z6 "
                "presence branch 1 - P/24), not an output of the theorem."
            ),
        },
        "numerics": {
            "count_laws": "exact rational (Fractions convolution / transfer states)",
            "poisson_masses": f"Decimal, {DECIMAL_PRECISION} significant digits",
            "tolerance_added_to_tv_before_comparison": str(NUMERIC_TOLERANCE),
        },
        "families": families,
        "all_families_within_bound": all_pass,
        "claim_boundary": {
            "counting_step": "theorem grade on the declared model class",
            "physical_realization": (
                "conditional: that settled-galaxy repair opportunities "
                "realize a refinement-stable collar opportunity model on "
                "codimension-one cut support stays a physical premise"
            ),
            "flux_recovery_closure": "conditional, unchanged by this witness",
            "lambda_collar_value": "branch data, not certified here",
        },
        "consumers": {
            "paper": "cosmology/oph_dark_matter_paper.tex (Poisson Activation, conditional theorem 2)",
            "simulator_galaxy_proxy": (
                "static_galaxy_measurement_report.json / "
                "STATIC_GALAXY_RAR_BTFR_RECEIPT consume the activation law "
                "p(x) = 1 - exp(-lambda_collar sqrt(x))"
            ),
            "sparc_comparison": (
                "code/dark_matter/scripts/sparc_rar_compare.py and the SPARC "
                "likelihood scripts consume the same activation law; the "
                "counting step adds no freedom beyond lambda_collar"
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    default_out = (
        pathlib.Path(__file__).resolve().parents[1]
        / "certificates"
        / "collar_poisson_witness_certificate.json"
    )
    parser.add_argument("--out", type=pathlib.Path, default=default_out)
    args = parser.parse_args()
    payload = compute()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    worst = min(Decimal(f["margin_decimal"]) for f in payload["families"])
    print(f"families: {len(payload['families'])}")
    print(f"all_families_within_bound: {payload['all_families_within_bound']}")
    print(f"smallest margin (bound - exact TV): {worst}")
    print(f"certificate: {args.out}")
    return 0 if payload["all_families_within_bound"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
