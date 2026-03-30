#!/usr/bin/env python3
"""Audit simple compare-only normalizer candidates beneath the neutrino attachment route."""

from __future__ import annotations

import argparse
import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[2]
READBACK_JSON = ROOT / "particles" / "runs" / "neutrino" / "realized_same_label_gap_defect_readback.json"
THEOREM_OBJECT_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_weighted_cycle_theorem_object.json"
BRIDGE_CANDIDATE_JSON = ROOT / "particles" / "runs" / "neutrino" / "neutrino_lambda_nu_bridge_candidate.json"
DEFAULT_OUT = ROOT / "particles" / "runs" / "neutrino" / "neutrino_attachment_normalizer_candidate_audit.json"


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _monomial_value(values: dict[str, float], keys: tuple[str, ...], exponents: tuple[float, ...]) -> float:
    value = 1.0
    for key, exponent in zip(keys, exponents):
        value *= values[key] ** exponent
    return value


def _formula_string(keys: tuple[str, ...], exponents: tuple[float, ...]) -> str:
    pieces = []
    for key, exponent in zip(keys, exponents):
        if exponent == 1:
            pieces.append(key)
        else:
            pieces.append(f"{key}^{exponent:g}")
    return " * ".join(pieces)


def build_payload(
    *,
    readback: dict[str, Any],
    theorem_object: dict[str, Any],
    bridge_candidate: dict[str, Any],
) -> dict[str, Any]:
    q = list((readback.get("q_e") or {}).values())
    gap = list((readback.get("gap_e") or {}).values())
    defect = list((readback.get("defect_e") or {}).values())
    mu = list((readback.get("mu_e") or {}).values())
    eta = list((readback.get("eta_e") or {}).values())
    target = float((bridge_candidate.get("compare_only_bridge_factor") or {}).get("F_nu_star"))
    symmetric_values = {
        "gamma": float(theorem_object["live_inputs"]["gamma"]),
        "eps": float(theorem_object["live_inputs"]["eps"]),
        "ratio_hat": float(theorem_object["live_outputs"]["dimensionless_ratio_dm21_over_dm32"]),
        "D_nu": float(theorem_object["live_inputs"]["D_nu"]),
        "p_nu": float(theorem_object["live_inputs"]["p_nu"]),
        "sum_q": float(sum(q)),
        "sum_gap": float(sum(gap)),
        "sum_defect": float(sum(defect)),
        "sum_mu": float(sum(mu)),
        "prod_q": float(math.prod(q)),
        "prod_gap": float(math.prod(gap)),
        "prod_defect": float(math.prod(defect)),
        "sum_mu_q": float(sum(a * b for a, b in zip(mu, q))),
        "sum_mu_gap": float(sum(a * b for a, b in zip(mu, gap))),
        "sum_mu_defect": float(sum(a * b for a, b in zip(mu, defect))),
        "sum_mu_eta2": float(sum(a * (b * b) for a, b in zip(mu, eta))),
    }
    exponent_choices = (-3.0, -2.0, -1.0, -0.5, 0.5, 1.0, 2.0, 3.0)

    candidates: list[dict[str, Any]] = []
    for width in (1, 2):
        for keys in itertools.combinations(symmetric_values.keys(), width):
            for exponents in itertools.product(exponent_choices, repeat=width):
                value = _monomial_value(symmetric_values, keys, exponents)
                rel_error = abs(value - target) / target
                candidates.append(
                    {
                        "formula": _formula_string(keys, exponents),
                        "keys": list(keys),
                        "exponents": list(exponents),
                        "value": value,
                        "relative_error": rel_error,
                        "complexity": width,
                    }
                )

    candidates.sort(key=lambda item: (item["complexity"], item["relative_error"]))
    single_factor = [item for item in candidates if item["complexity"] == 1][:10]
    two_factor = [item for item in candidates if item["complexity"] == 2][:10]
    best_simple = min(candidates, key=lambda item: (item["relative_error"], item["complexity"]))
    best_overall = min(candidates, key=lambda item: item["relative_error"])

    return {
        "artifact": "oph_neutrino_attachment_normalizer_candidate_audit",
        "generated_utc": _timestamp(),
        "status": "compare_only_normalizer_search",
        "public_promotion_allowed": False,
        "target_bridge_factor": {
            "name": "F_nu_star",
            "value": target,
            "source": "weighted_least_squares_compare_only_bridge_factor",
        },
        "invariant_pool": symmetric_values,
        "search_family": {
            "kind": "symmetric_one_or_two_factor_monomials",
            "exponents": list(exponent_choices),
            "objective": "minimize relative error to F_nu_star",
        },
        "best_simple_symmetric_candidate": best_simple,
        "best_overall_symmetric_candidate": best_overall,
        "top_single_factor_candidates": single_factor,
        "top_two_factor_candidates": two_factor,
        "working_observation": (
            "The live attachment-normalizer surface has a simple symmetric clue: the best two-factor "
            "candidate on the current branch is close to an inverse defect-sum times inverse square-root "
            "hierarchy ratio. This is compare-only research output and not a theorem-grade lambda_nu law."
        ),
        "hard_guard": {
            "status": "do_not_promote",
            "reason": (
                "Any candidate normalizer extracted here is a local search clue beneath the attachment route. "
                "It does not bypass the exact positive-rescaling no-go or emit lambda_nu without the missing "
                "attachment theorem."
            ),
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Build the neutrino attachment normalizer candidate audit.")
    parser.add_argument("--readback", default=str(READBACK_JSON))
    parser.add_argument("--theorem-object", default=str(THEOREM_OBJECT_JSON))
    parser.add_argument("--bridge-candidate", default=str(BRIDGE_CANDIDATE_JSON))
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    payload = build_payload(
        readback=_load_json(Path(args.readback)),
        theorem_object=_load_json(Path(args.theorem_object)),
        bridge_candidate=_load_json(Path(args.bridge_candidate)),
    )
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
