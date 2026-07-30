#!/usr/bin/env python3
"""Issue #589 horizon-record attachment verdict.

The bridge identity ``N_star = log D_star = A_dS/(4 ell_star^2)`` consumes a
unique-zero direct capacity closure from issue #505: a unique source-selected
``D_star``. The completed capacity source class carries the locked
non-identifiability verdict, so no unique ``D_star`` exists on the declared
class and the horizon-record identification has no capacity-side carrier to
identify. This producer consumes the locked verdict, fails closed on any
drift in it, and emits the typed negative exit declared by the issue:
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

SCHEMA = "oph.horizon_record_attachment_verdict.v1"
STATUS = "NOT_EVALUABLE_NO_HORIZON_RECORD_ATTACHMENT"
REQUIRED_DIRECT_N_SCHEMA = "oph.direct_n_closure_verdict.v2"
REQUIRED_DIRECT_N_STATUS = (
    "LOCKED_NONIDENTIFIABILITY_COMPLETED_CAPACITY_SOURCE_CLASS"
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


def build_verdict() -> dict[str, Any]:
    direct_n = _load(DIRECT_N_VERDICT_PATH)
    if direct_n.get("schema") != REQUIRED_DIRECT_N_SCHEMA:
        raise ValueError("direct N verdict schema drift")
    if direct_n.get("status") != REQUIRED_DIRECT_N_STATUS:
        raise ValueError("direct N verdict status drift")
    boundary = direct_n.get("comparison_boundary", {})
    if boundary.get("direct_numeric_N_emitted") is not False:
        raise ValueError("direct N verdict emitted a numeric N")
    if boundary.get("cosmological_comparison_permitted") is not False:
        raise ValueError("direct N verdict permitted a comparison")
    lift = _load(COMPLETE_LIFT_RECEIPT_PATH)
    if lift.get("scientific_verdict") != "COMPLETE_SOURCE_CLASS_NO_UNIQUE_SLACK_ZERO":
        raise ValueError("complete lift verdict drift")

    verdict = {
        "schema": SCHEMA,
        "issue": 589,
        "status": STATUS,
        "consumed_verdicts": {
            "direct_n_closure": REQUIRED_DIRECT_N_STATUS,
            "complete_lift": lift["scientific_verdict"],
        },
        "missing_antecedent": (
            "the bridge identity requires a unique source-selected D_star "
            "from a unique-zero direct closure; the completed declared capacity "
            "source class does not entail a unique slack zero, so no "
            "capacity-side carrier exists for the horizon identification"
        ),
        "unaffected_results": [
            "finite de Sitter capacity identities and shock normalization",
            "issue #503 inhabited-carrier tower campaign at its own scope",
            "the locked capacity non-identifiability theorem itself",
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
            "an additional named source law selecting a unique slack zero, "
            "a positive completed-class direct closure packet, and an "
            "inhabited issue #503 carrier on one common source tower"
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
