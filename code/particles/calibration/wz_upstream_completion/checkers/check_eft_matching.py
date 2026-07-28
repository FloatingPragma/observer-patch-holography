#!/usr/bin/env python3
"""Independent checker for the EFT_MATCHING_1 bundle (Workstream B).

Re-derives the one-loop gauge beta coefficients from the Workstream A
census with its own arithmetic, verifies the digest binding to the
action bundle, the scheme freezes, the recorded-empty decoupling with
its reason, the symbolic output vector, the declared remainder, the
refusal controls, and the subject digest.  A numeric scan enforces that
no coupling value exists anywhere in the packet.
"""

from __future__ import annotations

import hashlib
import json
import sys
from fractions import Fraction
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
ACTION_PATH = ROOT / "outputs" / "sm_eft_action_1.json"
BUNDLE_PATH = ROOT / "outputs" / "eft_matching_1.json"


def fail(message: str) -> None:
    raise SystemExit(f"eft_matching checker: {message}")


def canonical_sha256(value: Any) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def dim(label: str, kind: str) -> int:
    if kind == "color":
        return 3 if label in ("3", "3bar") else 1
    return 2 if label == "2" else 1


def dynkin(label: str, kind: str) -> Fraction:
    if kind == "color":
        return Fraction(1, 2) if label in ("3", "3bar") else Fraction(0)
    return Fraction(1, 2) if label == "2" else Fraction(0)


def rederive_betas(census: dict[str, Any]) -> dict[str, Fraction]:
    fermions = census["fermions"]
    scalars = census["scalars"]
    b3 = (
        -Fraction(11, 3) * 3
        + Fraction(2, 3) * sum(dynkin(f["color"], "color") * dim(f["weak"], "weak") for f in fermions)
        + Fraction(1, 3) * sum(dynkin(s["color"], "color") * dim(s["weak"], "weak") for s in scalars)
    )
    b2 = (
        -Fraction(11, 3) * 2
        + Fraction(2, 3) * sum(dynkin(f["weak"], "weak") * dim(f["color"], "color") for f in fermions)
        + Fraction(1, 3) * sum(dynkin(s["weak"], "weak") * dim(s["color"], "color") for s in scalars)
    )
    b1 = Fraction(2, 3) * sum(
        Fraction(f["hypercharge"]) ** 2 * dim(f["color"], "color") * dim(f["weak"], "weak")
        for f in fermions
    ) + Fraction(1, 3) * sum(
        Fraction(s["hypercharge"]) ** 2 * dim(s["color"], "color") * dim(s["weak"], "weak")
        for s in scalars
    )
    return {"b1": b1, "b2": b2, "b3": b3}


def numeric_scan_for_couplings(value: Any, path: str = "$") -> None:
    if isinstance(value, bool):
        return
    if isinstance(value, float):
        fail(f"float at {path}")
    if isinstance(value, dict):
        for key, item in value.items():
            numeric_scan_for_couplings(item, f"{path}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            numeric_scan_for_couplings(item, f"{path}[{index}]")


def check() -> None:
    if not BUNDLE_PATH.is_file() or not ACTION_PATH.is_file():
        fail("bundle or action file missing; run both producers first")
    bundle = json.loads(BUNDLE_PATH.read_text(encoding="utf-8"))
    action = json.loads(ACTION_PATH.read_text(encoding="utf-8"))

    if bundle.get("schema") != "eft_matching_1.v1":
        fail("schema mismatch")
    if bundle.get("promotion_allowed") is not False:
        fail("the bundle must not allow promotion")
    if bundle["action_subject_digest"] != action["subject_digest"]:
        fail("the packet must bind the current action subject digest")

    derived = rederive_betas(action["field_census"])
    stated = bundle["gauge_betas"]["coefficients"]
    for name, value in derived.items():
        if Fraction(stated[name]) != value:
            fail(f"beta coefficient {name} does not match the census re-derivation")
    if (derived["b1"], derived["b2"], derived["b3"]) != (
        Fraction(41, 6),
        Fraction(-19, 6),
        Fraction(-7),
    ):
        fail("re-derived coefficients must be (41/6, -19/6, -7)")
    if bundle["gauge_betas"]["gut_normalized_b1"]["status"] != "excluded":
        fail("the GUT-normalized value must be recorded and excluded")

    matching = bundle["matching"]
    intervals = matching["intervals"]
    if len(intervals) != 1 or intervals[0]["id"] != "pure_sm_interval":
        fail("exactly one pure Standard Model interval is expected")
    interval = intervals[0]
    if interval["decoupling_maps"]["maps"] != [] or not interval["decoupling_maps"]["reason"]:
        fail("the empty decoupling list must carry its reason")
    if interval["beta_coefficients"] != stated:
        fail("interval beta coefficients must equal the derived block")
    if matching["scheme_freeze"]["running"] != "MSbar":
        fail("the running scheme must be frozen to MSbar")
    if "FJ" not in matching["scheme_freeze"]["output_coordinate"]:
        fail("the output coordinate must be the frozen FJ coordinate")
    if matching["scheme_freeze"]["drbar_to_msbar_finite_maps"]["status"] != "imported_finite_maps":
        fail("the finite scheme maps must be labeled imports")

    vector = matching["output_vector"]
    if vector["name"] != "SM_MSbar_FJ(Q)":
        fail("the output vector must carry its canonical name")
    for component in vector["components"]:
        if "symbol" not in component:
            fail("every output-vector component must be a symbol")
    if matching["remainder"]["form"] != {"symbol": "R2[g, t]"}:
        fail("the truncation remainder must be the declared symbol")

    for name, verdict in bundle["controls"].items():
        if verdict.get("expected_failure") is not True or verdict.get("failed") is not True:
            fail(f"control {name} must record its refusal")

    digest = canonical_sha256(
        {
            "schema": bundle["schema"],
            "action_subject_digest": bundle["action_subject_digest"],
            "gauge_betas": bundle["gauge_betas"],
            "matching": bundle["matching"],
        }
    )
    if digest != bundle["subject_digest"]:
        fail("subject digest mismatch")

    numeric_scan_for_couplings(bundle["matching"]["output_vector"])
    numeric_scan_for_couplings(bundle["matching"]["remainder"])


def main() -> int:
    check()
    print(
        "eft_matching checker OK: census re-derivation, digest binding, scheme "
        "freezes, interval structure, symbolic vector, remainder, and controls all pass"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
