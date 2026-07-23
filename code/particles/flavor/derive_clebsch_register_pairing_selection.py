#!/usr/bin/env python3
"""Clebsch register selection: pairing theorem and frozen weight-set scan.

The down-type Clebsch lane (``derive_down_type_register_clebsch_lane.py``)
consumes the register pattern ``(b/tau, s/mu, d/e) = (1, 1/3, 3)`` under the
open premise gate ``CLEBSCH_REGISTER_SELECTION_THEOREM``.  This builder
splits that premise into its three parts and decides each on the corpus:

1. Pairing (theorem, conditional on the declared exterior matter contract of
   the super-Tannakian lift, issue #314).  The certified Yukawa sector
   couples ``d_c`` and ``e_c`` through the same scalar ``Sbar`` while ``u_c``
   couples through ``S``, and the forbidden channel ``[Q, S, d_c]`` carries a
   certified zero invariant line.  In the one-scalar package every register
   relation at ``mu_U`` that ties a quark Yukawa to a lepton Yukawa through a
   shared scalar therefore pairs down-type with charged leptons, and no
   register relation can pair down-type through the up scalar.

2. Weight set (conditional selection under two frozen constraints).  The
   weight alphabet is the multiplicative closure ``{1/3, 1, 3}`` of the two
   scalars the corpus emits for the color factor: ``N_c`` read from the
   certified charge spectrum and its inverse, the invariant probability
   measure of the transitive C3 color action.  The constraints are frozen in
   this file before any enumeration and no measured mass or angle enters:

   F1 (measure balance): the three register weights multiply to one - the
      counting measure ``N_c``, the invariant probability measure ``1/N_c``,
      and the unit total measure compose to the identity on the orbit.
   F2 (register faithfulness): the weight map is injective on the three
      generation registers - a collapsed map factors through a proper
      quotient of the transitive orbit and carries no register.

   The scan enumerates all 27 ordered triples.  Exactly one unordered
   multiset survives F1 and F2, and the theorem is emitted conditional on F1
   and F2; if the count were ever different, the builder emits the exclusion
   record instead and the premise stays open.

3. Generation-order assignment: open.  Which weight attaches to which
   generation requires the family attachment (#569); no source argument for
   the third-generation-weight-one assignment exists in the corpus, so the
   assignment stays a declared premise of the down-type lane.

Consequence: the open premise of the down-type lane shrinks from the full
selection theorem to ``GENERATION_REGISTER_ORDER`` (plus the unchanged MCPR
architecture condition).  The lane keeps tier T2; nothing here promotes.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
A5 = ROOT / "a5_closure"
MATTER_RECEIPT_JSON = A5 / "receipts" / "super_tannakian_matter_reference.receipt.json"
MATTER_MANIFEST_JSON = A5 / "manifests" / "super_tannakian_matter_reference.json"
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "flavor" / "clebsch_register_pairing_selection.json"
)

UP_CHANNEL = ("Q", "S", "u_c")
DOWN_CHANNEL = ("Q", "Sbar", "d_c")
LEPTON_CHANNEL = ("L", "Sbar", "e_c")
FORBIDDEN_CHANNEL = ("Q", "S", "d_c")

# Frozen before evaluation; recorded verbatim in the artifact.
FROZEN_CONSTRAINTS = {
    "F1_measure_balance": (
        "the three register weights multiply to one: counting measure N_c, "
        "invariant probability measure 1/N_c, and unit total measure compose "
        "to the identity on the transitive C3 orbit"
    ),
    "F2_register_faithfulness": (
        "the weight map is injective on the three generation registers; a "
        "collapsed map factors through a proper quotient of the transitive "
        "orbit and carries no register"
    ),
}

FORBIDDEN_SOLVE_PATH_INPUTS = [
    "m_d_measured",
    "m_s_measured",
    "m_b_measured",
    "cabibbo_measured",
    "PDG",
    "CODATA",
    "down_type_lane_predictions",
]


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def derive_weight_generators(receipt: dict[str, Any]) -> dict[str, Any]:
    """Read N_c from the certified charge spectrum; no declared constant."""

    spectrum = receipt.get("realized_package", {}).get("charge_spectrum")
    if not isinstance(spectrum, dict):
        raise SystemExit("fail closed: certified charge spectrum missing from receipt")
    n_c = spectrum.get("-2/3")
    checks = {
        "up_conjugate_multiplicity_is_n_c": n_c == 3,
        "down_conjugate_multiplicity_matches": spectrum.get("1/3") == n_c,
        "quark_doublet_multiplicity_is_two_n_c": spectrum.get("1/6") == 2 * n_c,
    }
    if not all(checks.values()):
        raise SystemExit(
            "fail closed: charge spectrum does not carry the color multiplicity "
            f"structure: {checks}"
        )
    return {
        "charge_spectrum": spectrum,
        "N_c": n_c,
        "checks": checks,
        "generators": {
            "counting_measure": str(Fraction(n_c)),
            "invariant_probability_measure": str(Fraction(1, n_c)),
            "unit_total_measure": "1",
        },
        "alphabet": [str(Fraction(1, n_c)), "1", str(Fraction(n_c))],
        "provenance": (
            "N_c is the multiplicity of the up-conjugate charge class in the "
            "certified one-generation charge spectrum of the matter lift; 1/N_c "
            "is the invariant probability measure of the transitive C3 color "
            "action; no numeric constant is declared by this builder"
        ),
    }


def pairing_theorem(receipt: dict[str, Any], manifest: dict[str, Any]) -> dict[str, Any]:
    sector = receipt.get("yukawa_sector", {})
    channels = {
        tuple(row.get("channel", [])): row.get("invariant_dimension")
        for row in sector.get("channels", [])
    }
    control = sector.get("forbidden_channel_control", {})
    contract = manifest.get("exterior_matter_contract", {})
    declared_channels = [tuple(ch) for ch in contract.get("yukawa_channels", [])]

    checks = {
        "declared_channel_list_matches_receipt": sorted(declared_channels)
        == sorted(channels.keys()),
        "one_scalar_package": contract.get("extra_scalars") == []
        and contract.get("one_scalar") == "weak_block",
        "down_couples_through_Sbar": channels.get(DOWN_CHANNEL) == 1,
        "lepton_couples_through_Sbar": channels.get(LEPTON_CHANNEL) == 1,
        "up_couples_through_S": channels.get(UP_CHANNEL) == 1,
        "down_through_up_scalar_has_no_line": (
            tuple(control.get("channel", [])) == FORBIDDEN_CHANNEL
            and control.get("invariant_dimension") == 0
        ),
    }
    passed = all(checks.values())
    return {
        "statement": (
            "in the declared one-scalar exterior package, d_c and e_c couple "
            "through the same scalar Sbar and u_c couples through S, while the "
            "channel [Q, S, d_c] carries a certified zero invariant line; "
            "therefore every register relation at mu_U tying a quark Yukawa to "
            "a lepton Yukawa through a shared scalar pairs down-type with "
            "charged leptons, and no register relation pairs down-type through "
            "the up scalar"
        ),
        "conditional_on": (
            "the declared exterior matter contract and branch premises of the "
            "conditional super-Tannakian matter lift (#314)"
        ),
        "channel_invariant_dimensions": {"/".join(k): v for k, v in channels.items()},
        "forbidden_channel_control": control,
        "checks": checks,
        "passed": passed,
    }


def weight_set_scan(generators: dict[str, Any]) -> dict[str, Any]:
    alphabet = [Fraction(value) for value in generators["alphabet"]]
    candidates = list(itertools.product(alphabet, repeat=3))
    survivors_f1 = [triple for triple in candidates if triple[0] * triple[1] * triple[2] == 1]
    survivors_f2 = [triple for triple in survivors_f1 if len(set(triple)) == 3]
    surviving_multisets = sorted(
        {tuple(sorted(triple)) for triple in survivors_f2}
    )
    unique = len(surviving_multisets) == 1

    return {
        "frozen_constraints": FROZEN_CONSTRAINTS,
        "constraint_order": ["F1_measure_balance", "F2_register_faithfulness"],
        "alphabet": [str(v) for v in alphabet],
        "candidate_count": len(candidates),
        "survivors_after_F1": {
            "count": len(survivors_f1),
            "unordered_multisets": sorted(
                {tuple(str(v) for v in sorted(t)) for t in survivors_f1}
            ),
        },
        "survivors_after_F2": {
            "count": len(survivors_f2),
            "unordered_multisets": [
                [str(v) for v in multiset] for multiset in surviving_multisets
            ],
        },
        "unique_unordered_survivor": unique,
        "surviving_weight_set": (
            [str(v) for v in surviving_multisets[0]] if unique else None
        ),
        "order_assignment": {
            "status": "open",
            "reason": (
                "which weight attaches to which generation requires the family "
                "attachment (#569); no source argument for the third-generation "
                "weight-one assignment exists in the corpus"
            ),
            "remaining_premise": "GENERATION_REGISTER_ORDER",
        },
    }


def build_artifact(
    receipt: dict[str, Any],
    manifest: dict[str, Any],
    input_hashes: dict[str, str],
) -> dict[str, Any]:
    generators = derive_weight_generators(receipt)
    pairing = pairing_theorem(receipt, manifest)
    scan = weight_set_scan(generators)

    pairing_closed = pairing["passed"]
    weight_set_closed = scan["unique_unordered_survivor"]

    if pairing_closed and weight_set_closed:
        status = "pairing_and_weight_set_selected_order_open"
        proof_status = "closed_conditional_selection_theorem"
    elif pairing_closed:
        status = "pairing_selected_weight_set_exclusion_record"
        proof_status = "partial_selection_exclusion_record"
    else:
        status = "selection_failed_premise_unchanged"
        proof_status = "exclusion_record_only"

    return {
        "artifact": "oph_clebsch_register_pairing_selection",
        "generated_utc": _timestamp(),
        "github_issues": [591, 546],
        "status": status,
        "proof_status": proof_status,
        "promotion_allowed": False,
        "tier_consequence": (
            "the down-type lane stays T2: the MCPR architecture premise alone "
            "forbids promotion; this artifact shrinks the lane's register "
            "premise, it does not lift the lane"
        ),
        "theorem_statement": (
            "Conditional on the declared exterior matter contract (#314) and "
            "the two frozen constraints F1 and F2, the register relation at "
            "mu_U pairs the down-type quarks with the charged leptons and its "
            "unordered weight set is exactly {1/3, 1, 3}. The generation-order "
            "assignment of the weights is not selected and remains the open "
            "premise GENERATION_REGISTER_ORDER."
            if pairing_closed and weight_set_closed
            else "The corpus does not select the declared register pattern; "
            "the exclusion record below is the deliverable and the full "
            "selection premise stays open."
        ),
        "pairing_theorem": pairing,
        "weight_generators": generators,
        "weight_set_scan": scan,
        "premise_bookkeeping": {
            "replaced_premise": "CLEBSCH_REGISTER_SELECTION_THEOREM",
            "closed_parts": (
                ["pairing (conditional on #314 contract)", "unordered weight set (conditional on F1, F2)"]
                if pairing_closed and weight_set_closed
                else ["pairing (conditional on #314 contract)"]
                if pairing_closed
                else []
            ),
            "remaining_open_parts": ["GENERATION_REGISTER_ORDER"]
            + ([] if weight_set_closed else ["CLEBSCH_WEIGHT_SET"])
            + ([] if pairing_closed else ["CLEBSCH_PAIRING"]),
            "unchanged_conditions": [
                "MCPR eight-register architecture (declared model input)",
                "conditionality of the #314 matter lift on its branch premises",
            ],
        },
        "solve_path_discipline": {
            "forbidden_inputs": FORBIDDEN_SOLVE_PATH_INPUTS,
            "reads_measured_masses": False,
            "reads_measured_ratios": False,
            "reads_down_type_lane_outputs": False,
            "numeric_inputs": (
                "N_c and 1/N_c, both read from the certified charge spectrum"
            ),
        },
        "dependency_audit": {
            "direct_input_paths": [
                "a5_closure/receipts/super_tannakian_matter_reference.receipt.json",
                "a5_closure/manifests/super_tannakian_matter_reference.json",
            ],
            "input_sha256": input_hashes,
            "no_target_leak": True,
        },
        "notes": [
            "F1 and F2 are declared constraints with stated source rationale; the selection is conditional on them and says so.",
            "The scan enumerates the full 27-triple candidate space and records survivors at each stage.",
            "Exactly one unordered multiset surviving is the landing condition; any other count downgrades this artifact to an exclusion record on the next run.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Clebsch register pairing theorem and frozen weight-set scan."
    )
    parser.add_argument("--receipt", default=str(MATTER_RECEIPT_JSON))
    parser.add_argument("--manifest", default=str(MATTER_MANIFEST_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    receipt_path = Path(args.receipt)
    manifest_path = Path(args.manifest)
    payload = build_artifact(
        _load_json(receipt_path),
        _load_json(manifest_path),
        {
            "matter_receipt": _sha256(receipt_path),
            "matter_manifest": _sha256(manifest_path),
        },
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {output}")
    print(f"status: {payload['status']}")
    if payload["status"] == "selection_failed_premise_unchanged":
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
