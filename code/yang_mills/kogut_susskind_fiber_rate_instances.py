#!/usr/bin/env python3
"""Exact-rational instance tables for the Kogut--Susskind fiber-rate module.

This producer ties the Lean module
``Lean/Screen/KogutSusskindFiberRateComparison.lean`` to the committed
finite transfer receipt's conventions.  It recomputes, in exact
``fractions.Fraction`` arithmetic, the two concrete instances proved in
Lean (two fibers with ratios ``(2, 3)`` at ``lambda = 1`` and three fibers
with ratios ``(1, 2, 3)`` at ``lambda = 1/2``): per-fiber rates
``c = lambda * (r + 1/r)``, flip invariance ``c(1/r) = c(r)``, the
heat-bath off-diagonal identity ``c / (1 + r^2) = lambda / r`` under the
``pi = Omega^2`` convention, the induced fiber weights, the full
subset-sum eigenvalue tables of the independent-product generator, the
exact spectral gaps ``min_l c_l``, and the ``2 lambda`` floor with its
equality case at ``r = 1``.  It then cross-checks the committed
``z2_finite_transfer_receipt.json``: the rate-identity string and the
exact ``2 lambda`` floor values of every Kogut--Susskind run.

Scope boundary.  The tables certify structure of the idealized
independent-product model only.  They are not a gauge-orbit quotient
result, not a production or custody record, not a physical receipt, and
not evidence for OPH; the producer digest in the output is a self-digest,
not authenticated custody.  The committed receipt's own Kogut--Susskind
rates live on the non-product orbit space and vary over orbits; only the
rate identity, the floor value, and the ``pi = Omega^2`` heat-bath
convention are shared and checked here.

Replay: ``python3 kogut_susskind_fiber_rate_instances.py`` from this
directory rewrites ``receipts/kogut_susskind_fiber_rate_instances.json``
deterministically.
"""

from __future__ import annotations

import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path

SCHEMA = "oph.yang_mills.kogut_susskind_fiber_rate_instances.v1"

HERE = Path(__file__).resolve().parent
COMMITTED_RECEIPT = HERE / "receipts" / "z2_finite_transfer_receipt.json"
OUTPUT = HERE / "receipts" / "kogut_susskind_fiber_rate_instances.json"

LEAN_MODULE = "Lean/Screen/KogutSusskindFiberRateComparison.lean"
LEAN_THEOREMS = [
    "fiberRate_flip",
    "two_mul_le_fiberRate",
    "fiberRate_eq_two_mul_iff",
    "fiberRate_heatBath_offdiagonal",
    "heatBath_offdiagonal_unique_rate",
    "dirichlet_domination_single",
    "dirichlet_domination_family",
    "twoPoint_kogutSusskind_domination",
    "generator_walsh",
    "dirichlet_constant_iff",
    "bestConstant_isGreatest_min_rate",
    "kogutSusskind_uniform_dirichlet_floor",
    "kogutSusskind_ksWeight_floor",
    "instanceA_rate_values",
    "instanceA_eigenvalue_table",
    "instanceA_gap",
    "instanceB_rate_values",
    "instanceB_top_eigenvalue",
    "instanceB_gap",
    "instanceB_attains_floor",
]


def fiber_rate(lam: Fraction, r: Fraction) -> Fraction:
    """Kogut--Susskind Doob collar rate ``lambda * (r + 1/r)``."""
    return lam * (r + 1 / r)


def frac(x: Fraction) -> str:
    return str(Fraction(x))


def instance_table(name: str, lam: Fraction, ratios: list[Fraction]) -> dict:
    n = len(ratios)
    rates = [fiber_rate(lam, r) for r in ratios]
    floor = 2 * lam
    # Scalar identities per fiber.
    for r, c in zip(ratios, rates):
        assert r > 0
        assert c == lam * (r * r + 1) / r, "rate identity failed"
        assert fiber_rate(lam, 1 / r) == c, "flip invariance failed"
        # Heat-bath off-diagonal identity under pi = Omega^2.
        assert c * (1 / (1 + r * r)) == lam / r, "off-diagonal identity failed"
        # Induced fiber weights.
        w_true = r * r / (1 + r * r)
        w_false = 1 / (1 + r * r)
        assert w_true + w_false == 1
        assert w_true / w_false == r * r
        # AM-GM floor and equality case.
        assert c >= floor
        assert (c == floor) == (r == 1)
    # Subset-sum eigenvalue table of the independent-product generator.
    eigenvalues = {}
    for k in range(n + 1):
        for subset in itertools.combinations(range(n), k):
            eigenvalues[",".join(map(str, subset)) or "empty"] = frac(
                sum((rates[l] for l in subset), Fraction(0))
            )
    gap = min(rates)
    nonempty_min = min(
        sum((rates[l] for l in subset), Fraction(0))
        for k in range(1, n + 1)
        for subset in itertools.combinations(range(n), k)
    )
    assert nonempty_min == gap, "gap is not the minimum fiber rate"
    assert gap >= floor
    attains = gap == floor
    assert attains == any(r == 1 for r in ratios)
    return {
        "instance": name,
        "n_fibers": n,
        "lambda": frac(lam),
        "ratios": [frac(r) for r in ratios],
        "rates": [frac(c) for c in rates],
        "fiber_weights_true_false": [
            [frac(r * r / (1 + r * r)), frac(1 / (1 + r * r))] for r in ratios
        ],
        "subset_sum_eigenvalues": eigenvalues,
        "spectral_gap": frac(gap),
        "floor_two_lambda": frac(floor),
        "gap_attains_floor": attains,
        "equality_ratio_present": any(r == 1 for r in ratios),
    }


def committed_cross_check() -> dict:
    data = json.loads(COMMITTED_RECEIPT.read_text(encoding="utf-8"))
    assert data["schema"] == "oph.yang_mills.z2_finite_transfer_receipt.v2"
    assert data["physical_clay_receipt"] is False
    assert data["grid_scope"]["universal_no_go"] is False
    runs = []
    for run in data["runs"]:
        if run["transfer"] != "kogut_susskind":
            continue
        floor_block = run["variable_rate_floor"]
        assert floor_block["identity"] == "c_l(o) = lambda * (r_l(o) + 1/r_l(o))"
        lam = Fraction(str(run["parameters"]["lam"]))
        assert Fraction(str(floor_block["lower_bound_value"])) == 2 * lam
        assert floor_block["numerical_min_respects_bound"] is True
        runs.append(
            {
                "L": run["L"],
                "lambda": frac(lam),
                "floor_two_lambda": frac(2 * lam),
            }
        )
    assert runs, "no Kogut--Susskind runs in the committed receipt"
    return {
        "committed_receipt": "receipts/z2_finite_transfer_receipt.json",
        "committed_receipt_sha256_at_read": hashlib.sha256(
            COMMITTED_RECEIPT.read_bytes()
        ).hexdigest(),
        "shared_conventions_checked": [
            "rate identity string c_l(o) = lambda * (r_l(o) + 1/r_l(o))",
            "floor value 2 * lambda exact for every Kogut--Susskind run",
            "numerical minimum respects the floor in every run",
            "pi = Omega^2 heat-bath conditional-weight convention",
        ],
        "kogut_susskind_runs": runs,
    }


def build() -> dict:
    instances = [
        instance_table("A_two_fiber", Fraction(1), [Fraction(2), Fraction(3)]),
        instance_table(
            "B_three_fiber",
            Fraction(1, 2),
            [Fraction(1), Fraction(2), Fraction(3)],
        ),
    ]
    receipt = {
        "schema": SCHEMA,
        "verdict": "STATIC_EXACT_RATIONAL_TABLE_CONFORMANT",
        "scope": (
            "independent-product model of two-point fibers with one constant "
            "rate per fiber; exact rational arithmetic; structure conformance "
            "only.  Not a gauge-orbit quotient result, not a production, "
            "provenance, or custody record, not a physical receipt, and not "
            "evidence for OPH."
        ),
        "producer_self_digest_sha256": hashlib.sha256(
            Path(__file__).read_bytes()
        ).hexdigest(),
        "self_digest_note": (
            "producer_self_digest_sha256 is a self-digest for replay "
            "convenience; it is not authenticated custody"
        ),
        "lean_module": LEAN_MODULE,
        "lean_theorems": LEAN_THEOREMS,
        "instances": instances,
        "committed_receipt_cross_check": committed_cross_check(),
        "open_step": (
            "variable-rate approximate-tensorization comparison for the "
            "committed non-product L=2,3 orbit laws, then the anisotropic "
            "Wilson-to-Hamiltonian scan.  GaugeOrbitQuotientGap formalizes "
            "invariant-subspace descent only for the independent product "
            "under moved-link symmetry; it does not instantiate these laws."
        ),
    }
    return receipt


def main() -> None:
    receipt = build()
    OUTPUT.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {OUTPUT}")
    for inst in receipt["instances"]:
        print(
            f"{inst['instance']}: lambda={inst['lambda']} "
            f"ratios={inst['ratios']} rates={inst['rates']} "
            f"gap={inst['spectral_gap']} floor={inst['floor_two_lambda']} "
            f"attains={inst['gap_attains_floor']}"
        )


if __name__ == "__main__":
    main()
