#!/usr/bin/env python3
"""Issue #589 horizon-record attachment verdict.

The bridge identity ``N_star = log D_star = A_dS/(4 ell_star^2)`` consumes a
unique-zero direct capacity closure from issue #505: a unique source-selected
``D_star``. The direct packet is not evaluable because universal all-rung
membership in the complete source contract and its executable-to-Lean bridge
remain open. The horizon-record identification therefore has no capacity-side
carrier to identify. This producer consumes that bounded verdict, fails closed
on any drift in it, and emits the typed negative exit declared by the issue:
``NOT_EVALUABLE_NO_HORIZON_RECORD_ATTACHMENT``.

The verdict does not touch the independent de Sitter and Einstein branch
campaigns: the finite capacity identities, the shock normalization, and the
issue #503 tower work stand at their own scope. The bridge equation
``Lambda ell_star^2 = 3 pi / N_star`` remains a downstream contract with no
evaluable left-hand input, and no cosmological payload is read or compared.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
DIRECT_N_VERDICT_PATH = RUNTIME / "direct_n_closure_verdict.json"
COMPLETE_LIFT_RECEIPT_PATH = RUNTIME / "complete_packet_capacity_lift_receipt.json"
OUTPUT_PATH = RUNTIME / "horizon_record_attachment_verdict.json"

SCHEMA = "oph.horizon_record_attachment_verdict.v2"
STATUS = "NOT_EVALUABLE_NO_HORIZON_RECORD_ATTACHMENT"
REQUIRED_DIRECT_N_SCHEMA = "oph.direct_n_closure_verdict.v3"
REQUIRED_DIRECT_N_STATUS = "NOT_EVALUABLE_INCOMPLETE_CAPACITY_SOURCE_ANTECEDENT"
REQUIRED_LIFT_SCHEMA = "oph.complete_packet_capacity_lift.v2"
REQUIRED_LIFT_VERDICT = (
    "BOUNDED_GENERATION_REGISTER_COUNTERMODEL__UNIVERSAL_MEMBERSHIP_OPEN"
)


def canonical_json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("ascii")


def tagged_sha256(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def _load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _pin(path: Path) -> dict[str, Any]:
    payload = path.read_bytes()
    return {
        "path": path.relative_to(HERE.parent.parent).as_posix(),
        "bytes": len(payload),
        "sha256": tagged_sha256(payload),
    }


def _verify_self_hash(payload: dict[str, Any], key: str, label: str) -> None:
    unhashed = dict(payload)
    committed = unhashed.pop(key, None)
    expected = tagged_sha256(canonical_json_bytes(unhashed))
    if committed != expected:
        raise ValueError(f"{label} self-pin mismatch")


def build_verdict() -> dict[str, Any]:
    direct_n = _load(DIRECT_N_VERDICT_PATH)
    if direct_n.get("schema") != REQUIRED_DIRECT_N_SCHEMA:
        raise ValueError("direct N verdict schema drift")
    if direct_n.get("status") != REQUIRED_DIRECT_N_STATUS:
        raise ValueError("direct N verdict status drift")
    _verify_self_hash(direct_n, "verdict_sha256", "direct N verdict")
    boundary = direct_n.get("comparison_boundary", {})
    if boundary.get("direct_numeric_N_emitted") is not False:
        raise ValueError("direct N verdict emitted a numeric N")
    if boundary.get("cosmological_comparison_permitted") is not False:
        raise ValueError("direct N verdict permitted a comparison")
    bounded = direct_n.get("bounded_generation_register_result", {})
    if bounded.get("universal_all_rung_membership_proved") is not False:
        raise ValueError("direct N verdict promoted universal membership")
    if bounded.get("executable_lean_membership_bridge_proved") is not False:
        raise ValueError("direct N verdict promoted executable-Lean bridge")
    if (
        bounded.get("complete_a1_a3_source_class_nonidentifiability_proved")
        is not False
    ):
        raise ValueError("direct N verdict promoted complete-source no-go")
    lift = _load(COMPLETE_LIFT_RECEIPT_PATH)
    if lift.get("schema") != REQUIRED_LIFT_SCHEMA:
        raise ValueError("bounded lift schema drift")
    if lift.get("scientific_verdict") != REQUIRED_LIFT_VERDICT:
        raise ValueError("bounded lift verdict drift")
    _verify_self_hash(lift, "receipt_sha256", "bounded lift receipt")
    source_status = lift.get("source_contract_status", {})
    if source_status.get("universal_all_rung_membership_proved") is not False:
        raise ValueError("bounded lift promoted universal membership")
    if source_status.get("direct_n_status") != REQUIRED_DIRECT_N_STATUS:
        raise ValueError("bounded lift direct N status drift")

    verdict = {
        "schema": SCHEMA,
        "issue": 589,
        "status": STATUS,
        "consumed_verdicts": {
            "direct_n_closure": REQUIRED_DIRECT_N_STATUS,
            "bounded_generation_register_packet": lift["scientific_verdict"],
        },
        "missing_antecedent": (
            "the bridge identity requires a unique source-selected D_star "
            "from a unique-zero direct closure; universal all-rung membership "
            "in the complete A1-A3 capacity source contract and the "
            "executable-to-Lean membership bridge remain open, so no "
            "source-selected capacity carrier exists for the horizon identification"
        ),
        "unaffected_results": [
            "finite de Sitter capacity identities and shock normalization",
            "issue #503 inhabited-carrier tower campaign at its own scope",
            "the exact bounded generation-register capacity arithmetic",
        ],
        "bridge_equation_status": (
            "Lambda ell_star^2 = 3 pi / N_star stays a downstream contract "
            "with no evaluable left-hand input; it is not evidence in either "
            "direction"
        ),
        "comparison_boundary": {
            "cosmological_payload_read": False,
            "lambda_comparison_permitted": False,
            "forecast_entry_permitted": False,
        },
        "reopen_condition": (
            "a universal all-rung source-membership theorem, an "
            "executable-to-Lean membership bridge, a positive direct closure "
            "selecting one D_star, and an inhabited issue #503 carrier on one "
            "common source tower"
        ),
        "parent_pins": [
            _pin(DIRECT_N_VERDICT_PATH),
            _pin(COMPLETE_LIFT_RECEIPT_PATH),
        ],
    }
    verdict["verdict_sha256"] = tagged_sha256(canonical_json_bytes(verdict))
    return verdict


def write_runtime() -> None:
    OUTPUT_PATH.write_bytes(canonical_json_bytes(build_verdict()))


def verify_runtime() -> None:
    if OUTPUT_PATH.read_bytes() != canonical_json_bytes(build_verdict()):
        raise SystemExit("horizon-record attachment verdict is stale")
    if _load(OUTPUT_PATH)["status"] != STATUS:
        raise SystemExit("horizon-record attachment verdict status drift")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write:
        write_runtime()
    if args.verify:
        verify_runtime()
    if not args.write and not args.verify:
        print(canonical_json_bytes(build_verdict()).decode("ascii"), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
