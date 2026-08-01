#!/usr/bin/env python3
"""Independent exact replay of the issue-649 attachment receipt."""

from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
RECEIPT = HERE / "runtime" / "global_capacity_attachment_receipt.json"
SOURCE = ROOT / "code" / "cosmology" / "manifests" / "edge_center_clock_certificate.json"
PROJECTION = HERE / "manifests" / "global_capacity_attachment_source_projection.json"


def canonical(value: Any) -> bytes:
    return (
        json.dumps(value, ensure_ascii=True, separators=(",", ":"), sort_keys=True)
        + "\n"
    ).encode("ascii")


def tagged(raw: bytes) -> str:
    return "sha256:" + hashlib.sha256(raw).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise SystemExit(message)


def verify(path: Path = RECEIPT) -> None:
    value = json.loads(path.read_text(encoding="ascii"))
    require(value.get("schema") == "oph.global_capacity_attachment.v1", "schema drift")
    require(value.get("issue") == 649, "issue drift")
    require(
        value.get("status") == "BOUNDED_NONIDENTIFIABLE_GLOBAL_COMPOSITION",
        "status drift",
    )
    source_raw = SOURCE.read_bytes()
    source = json.loads(source_raw)
    projection_raw = PROJECTION.read_bytes()
    projection = json.loads(projection_raw)
    require(
        source.get("schema") == "oph.edge_center_clock_certificate.v3"
        and source.get("status")
        == "conditional_edge_center_arithmetic_with_open_source_and_clock_gates",
        "finite source status drift",
    )
    require(
        source.get("generator", {}).get("operational_clock_bound") is False,
        "finite source clock promotion",
    )
    reserve_input = source.get("antecedents", {}).get(
        "reserve_trace_branch_input", {}
    )
    require(
        reserve_input.get("binding_theorem_present") is False
        and reserve_input.get("derived_from_collar_counts") is False,
        "finite source branch boundary drift",
    )
    p_interval = source.get("inputs", {}).get("P_certified_enclosure", {})
    p_lo = Fraction(p_interval.get("lo"))
    p_hi = Fraction(p_interval.get("hi"))
    require(0 < p_lo <= p_hi < 4, "finite source q domain drift")
    require(
        projection
        == {
            "P_certified_enclosure": p_interval,
            "finite_source_schema": source["schema"],
            "finite_source_status": source["status"],
            "operational_clock_bound": source["generator"]["operational_clock_bound"],
            "parent_path": SOURCE.relative_to(ROOT).as_posix(),
            "presence_statement": source["antecedents"]["presence_reading"]["statement"],
            "reserve_trace_branch": {
                "binding_theorem_present": reserve_input["binding_theorem_present"],
                "derived_from_collar_counts": reserve_input["derived_from_collar_counts"],
                "statement": reserve_input["statement"],
            },
            "schema": "oph.global_capacity_attachment_source_projection.v1",
            "scope": (
                "Target-free projection of the conditional local survival datum. "
                "The full parent receipt and its negative controls are checked only "
                "by independent replay."
            ),
        },
        "source projection drift",
    )
    require(
        value.get("source_pin")
        == {
            "path": PROJECTION.relative_to(ROOT).as_posix(),
            "bytes": len(projection_raw),
            "sha256": tagged(projection_raw),
        },
        "source pin drift",
    )
    cone = value.get("producer_cone", {})
    require(
        cone.get("cosmological_target_payload_used_in_countermodel") is False,
        "target used",
    )
    require(
        cone.get("desired_capacity_used_in_countermodel") is False,
        "capacity target used",
    )
    require(cone.get("completion_selected_from_target") is False, "target selection")
    require(cone.get("full_parent_receipt_not_in_producer_cone") is True, "parent cone")

    rows = value.get("exact_witnesses")
    require(isinstance(rows, list) and len(rows) == 3, "witness set drift")
    for row in rows:
        q = Fraction(row["q"])
        capacity = Fraction(row["capacity"])
        survival = 1 - q
        require(0 < q < Fraction(1, 6), "q outside exact witness domain")
        require(Fraction(row["survival"]) == survival, "survival drift")
        completions = {item["completion"]: item for item in row["completions"]}
        require(set(completions) == {"neutral", "multiplicative"}, "completion set drift")
        neutral_expected = [capacity for _ in range(4)]
        multiplicative_expected = [survival**m * capacity for m in range(4)]
        require(
            [Fraction(x) for x in completions["neutral"]["values_m_0_to_3"]]
            == neutral_expected,
            "neutral values drift",
        )
        require(
            [Fraction(x) for x in completions["multiplicative"]["values_m_0_to_3"]]
            == multiplicative_expected,
            "multiplicative values drift",
        )
        for name, action in (
            ("neutral", lambda m, c: c),
            ("multiplicative", lambda m, c: survival**m * c),
        ):
            for left in range(4):
                for right in range(4 - left):
                    require(
                        action(left + right, capacity)
                        == action(left, action(right, capacity)),
                        f"{name} composition failed",
                    )
            checks = completions[name].get("refinement_partition_checks")
            require(isinstance(checks, list) and len(checks) == 4, "refinement checks drift")
            for check in checks:
                pieces = check.get("pieces")
                require(
                    isinstance(pieces, list)
                    and all(isinstance(piece, int) and piece > 0 for piece in pieces),
                    "refinement partition drift",
                )
                refined = capacity
                for piece in reversed(pieces):
                    refined = action(piece, refined)
                coarse = action(sum(pieces), capacity)
                require(refined == coarse, f"{name} refinement failed")
                require(
                    Fraction(check["coarse_value"]) == coarse
                    and Fraction(check["iterated_refined_value"]) == refined
                    and check.get("commutes") is True,
                    f"{name} refinement receipt drift",
                )
        require(neutral_expected[1] != multiplicative_expected[1], "countermodels agree")
        factors = row["blocked_event_factors"]
        expected_factors = {Fraction(1), 1 - q, 1 - 6 * q}
        require({Fraction(x) for x in factors.values()} == expected_factors, "factor menu drift")
        require(len(expected_factors) == 3, "factor menu collapsed")
        blocked = {
            item["completion"]: item for item in row["blocked_event_completions"]
        }
        require(
            set(blocked)
            == {
                "no_capacity_action",
                "one_class_projector",
                "six_class_total_projector",
            },
            "blocked-event completion set drift",
        )
        blocked_factors = {
            "no_capacity_action": Fraction(1),
            "one_class_projector": 1 - q,
            "six_class_total_projector": 1 - 6 * q,
        }
        one_cut_values = set()
        for name, factor in blocked_factors.items():
            item = blocked[name]
            require(Fraction(item["factor"]) == factor, f"{name} factor drift")
            expected_values = [factor**m * capacity for m in range(4)]
            require(
                [Fraction(x) for x in item["values_m_0_to_3"]] == expected_values,
                f"{name} action values drift",
            )
            require(
                item.get("identity_at_zero_cuts") is True
                and item.get("positive") is True
                and item.get("permutation_invariant") is True
                and item.get("disconnected_cut_composition") is True
                and item.get("refinement_regrouping_invariant") is True,
                f"{name} action law drift",
            )
            for left in range(4):
                for right in range(4 - left):
                    require(
                        factor ** (left + right) * capacity
                        == factor**left * (factor**right * capacity),
                        f"{name} composition failed",
                    )
            one_cut_values.add(expected_values[1])
        require(len(one_cut_values) == 3, "blocked-event actions agree")
        require(
            row.get("blocked_event_actions_pairwise_distinct_at_one_cut") is True,
            "blocked-event action verdict drift",
        )

    boundary = value.get("comparison_boundary")
    require(isinstance(boundary, dict) and all(v is False for v in boundary.values()), "promotion")
    blocked = value.get("blocked_event_verdict")
    require(
        blocked.get("one_class_selected") is False
        and blocked.get("six_class_total_selected") is False
        and blocked.get("no_capacity_action_selected") is False,
        "blocked-event selection drift",
    )
    require(
        blocked.get("status")
        == "SOURCE_DOES_NOT_SELECT_AMONG_DECLARED_BLOCKED_EVENT_ACTIONS",
        "blocked-event verdict drift",
    )
    expected_hash = value.get("receipt_sha256")
    unhashed = dict(value)
    unhashed.pop("receipt_sha256", None)
    require(expected_hash == tagged(canonical(unhashed)), "receipt hash drift")


if __name__ == "__main__":
    verify()
    print("GLOBAL_CAPACITY_ATTACHMENT_INDEPENDENT_VALID")
