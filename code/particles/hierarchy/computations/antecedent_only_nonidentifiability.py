#!/usr/bin/env python3
"""Antecedent-only non-identifiability witnesses for hierarchy issue #518.

The existing ``N_CRC^EW`` row fixes its output by imposing the bridge residual
whose solution it reports.  The existing ``epsilon_H=0`` row accepts the two
commuting-square defects as inputs.  Neither is a producer.

This module removes those target conditions and asks what the remaining typed
antecedents determine.  Exact finite countermodels show:

* the structural Banach/readback conditions admit distinct capacity fixed
  points on one unchanged source packet; and
* the scalar hierarchy packet plus map type signatures admits both commuting
  and noncommuting comparison-square completions.

The result is a no-go theorem about the present premise set.  It does not emit
a physical capacity or a Higgs naturality value.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Mapping


HERE = Path(__file__).resolve().parent
DEFAULT_OUT = (
    HERE.parent
    / "certificates"
    / "antecedent_only_nonidentifiability_receipt.json"
)

SOURCE_PACKET = {
    "branch": "source_audit_pixel_branch",
    "P": (
        "1.63097209585889737696451390350695562847912625483895268486516"
    ),
    "alpha_U": (
        "0.04112424744181668514088993388965971943770774203135879"
    ),
    "beta_EW": 4,
    "m_rep": 24,
    "capacity_structural_axioms": [
        "finite log-coordinate interval",
        "total monotone self-map",
        "strict contraction",
        "unique fixed point",
    ],
    "naturality_type_signatures": {
        "rho_sH": "Q_s -> Q_H",
        "n_s": "Q_s -> Q_s",
        "n_H": "Q_H -> Q_H",
        "h_s": "Q_s -> O_s",
        "chi_sH": "O_s -> O_H",
        "h_H": "Q_H -> O_H",
    },
    "excluded_target_conditions": [
        "B_EW(P,N)=0",
        "Pi_EW(P,N)=4P",
        "N_CRC^EW",
        "epsilon_n",
        "epsilon_h",
        "epsilon_H",
        "measured weak scale",
        "measured Lambda or Planck capacity",
    ],
}


def canonical_hash(payload: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _q(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def capacity_completion(
    *,
    completion_id: str,
    fixed_log_capacity: Fraction,
    domain_lo: Fraction = Fraction(0),
    domain_hi: Fraction = Fraction(4),
    averaging_weight: Fraction = Fraction(1, 2),
) -> dict[str, Any]:
    """Build ``F_c(x)=(1-lambda)x+lambda*c`` in exact arithmetic."""
    if not domain_lo < fixed_log_capacity < domain_hi:
        raise ValueError("fixed point must lie in the interior of the domain")
    if not 0 < averaging_weight <= 1:
        raise ValueError("averaging weight must lie in (0,1]")
    one_minus = 1 - averaging_weight
    image_lo = one_minus * domain_lo + averaging_weight * fixed_log_capacity
    image_hi = one_minus * domain_hi + averaging_weight * fixed_log_capacity
    return {
        "completion_id": completion_id,
        "antecedent_fingerprint": canonical_hash(SOURCE_PACKET),
        "map_family": "F_c(x)=(1-lambda)*x+lambda*c",
        "domain_log_capacity": [_q(domain_lo), _q(domain_hi)],
        "fixed_log_capacity": _q(fixed_log_capacity),
        "averaging_weight_lambda": _q(averaging_weight),
        "image_log_capacity": [_q(image_lo), _q(image_hi)],
        "derivative": _q(one_minus),
        "monotone": one_minus >= 0,
        "strict_contraction": abs(one_minus) < 1,
        "self_map": domain_lo <= image_lo <= image_hi <= domain_hi,
        "fixed_point_residual": "0",
        "unique_fixed_point": True,
    }


def _compose(
    outer: Mapping[str, str],
    inner: Mapping[str, str],
    key: str,
) -> str:
    return outer[inner[key]]


def evaluate_naturality_completion(completion: Mapping[str, Any]) -> dict[str, int]:
    """Evaluate both finite commuting-square defects without supplied flags."""
    q_source = list(completion["sets"]["Q_s"])
    rho = completion["maps"]["rho_sH"]
    n_s = completion["maps"]["n_s"]
    n_h = completion["maps"]["n_H"]
    h_s = completion["maps"]["h_s"]
    chi = completion["maps"]["chi_sH"]
    h_h = completion["maps"]["h_H"]

    normal_defects = [
        int(_compose(rho, n_s, x) != _compose(n_h, rho, x))
        for x in q_source
    ]
    obstruction_defects = [
        int(_compose(chi, h_s, x) != _compose(h_h, rho, x))
        for x in q_source
    ]
    epsilon_n = max(normal_defects, default=0)
    epsilon_h = max(obstruction_defects, default=0)
    return {
        "epsilon_n": epsilon_n,
        "epsilon_h": epsilon_h,
        "epsilon_H": max(epsilon_n, epsilon_h),
    }


def naturality_completion(
    *,
    completion_id: str,
    n_h: Mapping[str, str],
    h_h: Mapping[str, str],
) -> dict[str, Any]:
    """Build and evaluate one exact two-point comparison-square completion."""
    completion: dict[str, Any] = {
        "completion_id": completion_id,
        "antecedent_fingerprint": canonical_hash(SOURCE_PACKET),
        "metric": "discrete metric on finite codomains",
        "sets": {
            "Q_s": ["s0", "s1"],
            "Q_H": ["h0", "h1"],
            "O_s": ["os0", "os1"],
            "O_H": ["oh0", "oh1"],
        },
        "maps": {
            "rho_sH": {"s0": "h0", "s1": "h1"},
            "n_s": {"s0": "s0", "s1": "s1"},
            "n_H": dict(n_h),
            "h_s": {"s0": "os0", "s1": "os1"},
            "chi_sH": {"os0": "oh0", "os1": "oh1"},
            "h_H": dict(h_h),
        },
    }
    completion["evaluated_defects"] = evaluate_naturality_completion(completion)
    return completion


def build_receipt() -> dict[str, Any]:
    source_fingerprint = canonical_hash(SOURCE_PACKET)
    capacity_models = [
        capacity_completion(
            completion_id="capacity_completion_c1",
            fixed_log_capacity=Fraction(1),
        ),
        capacity_completion(
            completion_id="capacity_completion_c3",
            fixed_log_capacity=Fraction(3),
        ),
    ]
    naturality_models = [
        naturality_completion(
            completion_id="commuting_completion",
            n_h={"h0": "h0", "h1": "h1"},
            h_h={"h0": "oh0", "h1": "oh1"},
        ),
        naturality_completion(
            completion_id="normal_form_noncommuting_completion",
            n_h={"h0": "h1", "h1": "h0"},
            h_h={"h0": "oh0", "h1": "oh1"},
        ),
        naturality_completion(
            completion_id="obstruction_noncommuting_completion",
            n_h={"h0": "h0", "h1": "h1"},
            h_h={"h0": "oh1", "h1": "oh0"},
        ),
    ]
    return {
        "artifact": "oph_issue_518_antecedent_only_nonidentifiability",
        "issue": 518,
        "arithmetic": "fractions.Fraction and exhaustive finite-map evaluation",
        "source_packet": SOURCE_PACKET,
        "source_packet_sha256": source_fingerprint,
        "target_artifacts_consumed": [],
        "capacity_nonidentifiability": {
            "theorem": (
                "The target-free scalar packet and generic Banach/readback "
                "conditions do not identify a physical capacity."
            ),
            "countermodels": capacity_models,
            "same_antecedents": True,
            "distinct_fixed_log_capacities": True,
            "candidate_family_cardinality": "continuum: every c in the domain interior",
            "result": "NO_UNIQUE_CAPACITY_PRODUCER_FROM_DECLARED_ANTECEDENTS",
            "minimal_missing_object": (
                "a source-derived public-record/readback map F with physical "
                "semantics and an independently computed fixed point; a "
                "separate theorem may then test, rather than impose, B_EW=0"
            ),
            "not_claimed": [
                "either finite countermodel is the physical capacity",
                "the source-derived readback map cannot exist",
                "the EW bridge residual is false",
            ],
        },
        "naturality_nonidentifiability": {
            "theorem": (
                "The scalar hierarchy packet and map type signatures do not "
                "identify epsilon_H without executable source-derived maps."
            ),
            "countermodels": naturality_models,
            "same_antecedents": True,
            "evaluated_epsilon_H_values": [0, 1, 1],
            "result": "NO_UNIQUE_NATURALITY_DEFECT_FROM_DECLARED_ANTECEDENTS",
            "minimal_missing_object": (
                "source-derived finite or symbolic definitions of rho_sH, n_s, "
                "n_H, h_s, chi_sH, and h_H on a common domain, followed by an "
                "independent evaluation of both commuting-square residuals"
            ),
            "not_claimed": [
                "the physical comparison square has nonzero defect",
                "an exact naturality theorem is impossible after the maps exist",
                "the two-point completions are physical Higgs models",
            ],
        },
        "promotion_boundary": {
            "promotable_result": (
                "non-identifiability of the two targets from the present "
                "target-free premise packet"
            ),
            "not_promotable": [
                "N_CRC^EW as a physical capacity",
                "epsilon_H=0 as a physical naturality result",
            ],
        },
        "pass": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--print-json", action="store_true")
    args = parser.parse_args()
    payload = build_receipt()
    text = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.check:
        tracked = json.loads(Path(args.output).read_text(encoding="utf-8"))
        if tracked != payload:
            print("tracked non-identifiability receipt is stale")
            return 1
    else:
        Path(args.output).write_text(text, encoding="utf-8")
    if args.print_json:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
