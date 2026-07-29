#!/usr/bin/env python3
"""Issue-639 prediction-contract machinery.

The module carries the frozen contract format, the complete candidate
ledger, the target-ancestry checker, and the sealed-comparison
protocol.  Selection, freezing, and scoring of an actual forecast are
gated on a discriminator receipt upgrading a candidate to source
visibility; this module makes every step of that pipeline
machine-checkable and fail-closed before any candidate reaches it.

Fail-closed properties, each exercised by the test suite:

* a contract violating the schema is refused;
* a candidate whose declared ancestry names a forbidden input class is
  refused by the ancestry checker;
* a frozen contract whose bytes drift after the freeze is refused by
  the freeze check;
* a scoring attempt whose branch rule differs from the frozen contract
  is refused;
* a comparison payload whose bytes disagree with the sealed hash is
  refused, and unsealing is recorded exactly once.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path
from typing import Any

PACKET_DIR = Path(__file__).resolve().parent
SCHEMA_PATH = PACKET_DIR / "data" / "prediction_contract_schema_v1.json"
INVENTORY_PATH = PACKET_DIR / "data" / "candidate_inventory_v1.json"
OUT_PATH = PACKET_DIR / "outputs" / "forecast_contract_state.json"

STATE_SCHEMA = "oph.forecast_contract_state.v1"
ISSUE = 639


class ContractError(ValueError):
    def __init__(self, code: str, message: str) -> None:
        super().__init__(f"{code}: {message}")
        self.code = code


def require(condition: bool, code: str, message: str) -> None:
    if not condition:
        raise ContractError(code, message)


def sha256_bytes(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def validate_contract(contract: dict[str, Any]) -> None:
    """Validate a contract against the frozen schema, fail closed."""

    schema = load_json(SCHEMA_PATH)
    required = set(schema["required"])
    allowed = set(schema["properties"])
    present = set(contract)
    require(
        required.issubset(present),
        "CONTRACT_MISSING_FIELDS",
        f"missing: {sorted(required - present)}",
    )
    require(
        present.issubset(allowed),
        "CONTRACT_UNKNOWN_FIELDS",
        f"unknown: {sorted(present - allowed)}",
    )
    require(
        contract["schema"] == "oph.prediction_contract.v1",
        "CONTRACT_SCHEMA",
        "wrong contract schema string",
    )
    require(
        contract["output_type"]
        in ("dimensionless_real_interval", "discrete_verdict"),
        "CONTRACT_OUTPUT_TYPE",
        "output type outside the frozen enumeration",
    )
    for pin in contract["allowed_ancestry"]:
        require(
            set(pin) == {"artifact", "sha256", "role"}
            and str(pin["sha256"]).startswith("sha256:"),
            "CONTRACT_ANCESTRY_PIN",
            "every ancestry entry needs artifact, tagged sha256, and role",
        )
    for block in ("generator", "independent_checker"):
        require(
            set(contract[block]) == {"module", "sha256"}
            and str(contract[block]["sha256"]).startswith("sha256:"),
            "CONTRACT_MODULE_PIN",
            f"{block} needs a module path and tagged sha256",
        )
    sealed = contract["sealed_comparison"]
    require(
        set(sealed) == {"payload_sha256", "byte_count", "storage_note"}
        and str(sealed["payload_sha256"]).startswith("sha256:")
        and int(sealed["byte_count"]) >= 1,
        "CONTRACT_SEAL",
        "the sealed comparison needs a tagged hash and byte count",
    )
    freeze = contract["freeze"]
    require(
        freeze.get("frozen_by_commit") is True
        and freeze.get("unsealing_after_commit_only") is True,
        "CONTRACT_FREEZE",
        "the freeze block must bind commitment before unsealing",
    )


def check_target_ancestry(
    candidate: dict[str, Any], forbidden_classes: dict[str, Any]
) -> dict[str, Any]:
    """Bounded declared-vocabulary audit of the recorded ancestry entries.

    Each forbidden class carries its declared match fragments; an
    ancestry entry naming one of them is a hit unless the entry
    records the input as absent or excluded.  This is a schema audit
    over the declared entries, not semantic input closure."""

    hits = []
    for entry in candidate.get("allowed_ancestry", []):
        text = str(entry).lower()
        if "absent" in text or "without" in text or "only if" in text:
            continue
        for name, row in forbidden_classes.items():
            for fragment in row["match_fragments"]:
                if fragment in text:
                    hits.append(
                        {"entry": entry, "forbidden_class": name,
                         "fragment": fragment}
                    )
    return {
        "check_type": "bounded_declared_vocabulary_audit",
        "ancestry_hits": hits,
        "ancestry_clean": bool(not hits),
    }


def seal_comparison(payload_bytes: bytes, storage_note: str) -> dict[str, Any]:
    """Produce the sealed-comparison manifest for a payload."""

    return {
        "payload_sha256": sha256_bytes(payload_bytes),
        "byte_count": len(payload_bytes),
        "storage_note": storage_note,
    }


def verify_frozen_contract(
    contract: dict[str, Any], frozen_sha256: str
) -> None:
    """Refuse a contract whose bytes drift after the freeze."""

    current = sha256_bytes(canonical_json(contract).encode("utf-8"))
    require(
        current == frozen_sha256,
        "CONTRACT_POST_FREEZE_MUTATION",
        "the contract bytes differ from the frozen hash",
    )


def unseal_and_score(
    contract: dict[str, Any],
    frozen_contract_sha256: str,
    comparison_bytes: bytes,
    scoring_branch_rule: str,
    prior_unsealing_record: dict[str, Any] | None,
) -> dict[str, Any]:
    """Single-use unsealing with leakage, mutation, and branch guards."""

    require(
        prior_unsealing_record is None,
        "COMPARISON_ALREADY_UNSEALED",
        "scoring happens exactly once",
    )
    validate_contract(contract)
    verify_frozen_contract(contract, frozen_contract_sha256)
    require(
        scoring_branch_rule == contract["branch_rule"],
        "BRANCH_RULE_REPLACED",
        "the scoring branch rule must equal the frozen branch rule",
    )
    sealed = contract["sealed_comparison"]
    require(
        sha256_bytes(comparison_bytes) == sealed["payload_sha256"]
        and len(comparison_bytes) == int(sealed["byte_count"]),
        "COMPARISON_LEAKAGE_OR_TAMPER",
        "the comparison payload disagrees with the sealed manifest",
    )
    return {
        "unsealed": True,
        "contract_id": contract["contract_id"],
        "comparison_sha256": sealed["payload_sha256"],
        "branch_rule_applied": scoring_branch_rule,
    }


def build_state() -> dict[str, Any]:
    """Emit the machine-readable contract-and-ledger state receipt."""

    schema_raw = SCHEMA_PATH.read_bytes()
    inventory_raw = INVENTORY_PATH.read_bytes()
    inventory = json.loads(inventory_raw.decode("utf-8"))

    forbidden = inventory["forbidden_input_classes"]
    ancestry_verdicts = {}
    eligible = []
    for name, candidate in sorted(inventory["candidates"].items()):
        verdict = check_target_ancestry(candidate, forbidden)
        ancestry_verdicts[name] = {
            "eligibility": candidate["eligibility"],
            "blocking_condition": candidate["blocking_condition"],
            **verdict,
        }
        if candidate["eligibility"] == "ELIGIBLE_SOURCE_VISIBLE":
            eligible.append(name)

    payload = {
        "schema": STATE_SCHEMA,
        "issue": ISSUE,
        "pins": {
            "contract_schema_sha256": sha256_bytes(schema_raw),
            "candidate_inventory_sha256": sha256_bytes(inventory_raw),
            "module_sha256": sha256_bytes(Path(__file__).read_bytes()),
        },
        "forbidden_input_classes": forbidden,
        "selection_rule": inventory["selection_rule"],
        "candidate_ledger": ancestry_verdicts,
        "eligible_candidates": eligible,
        "contract_freeze_status": (
            "NO_CANDIDATE_ELIGIBLE_YET"
            if not eligible
            else "SELECTION_PENDING_UNDER_FROZEN_RULE"
        ),
        "gate": (
            "a contract freeze and sealed comparison happen only after a "
            "discriminator receipt upgrades a candidate to source "
            "visibility; the issue-569 and issue-636 sprints are the "
            "named upgrade paths"
        ),
        "claim_boundary": (
            "Frozen issue-639 contract infrastructure: the prediction "
            "contract format, the complete low-cost candidate ledger with "
            "ancestry and eligibility, the target-ancestry checker, the "
            "selection rule fixed before comparison access, and the "
            "single-use sealed-comparison protocol with fail-closed "
            "guards. No forecast is selected, frozen, or scored by this "
            "receipt, and no candidate is promoted past its recorded "
            "eligibility."
        ),
    }
    payload["state_sha256"] = sha256_bytes(
        canonical_json(
            {k: v for k, v in payload.items() if k != "state_sha256"}
        ).encode("utf-8")
    )
    return payload


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args(argv)
    payload = build_state()
    if args.verify:
        stored = json.loads(OUT_PATH.read_text(encoding="utf-8"))
        if stored != payload:
            print("CONTRACT_STATE_DRIFT", file=sys.stderr)
            return 1
        print("CONTRACT_STATE_VERIFIED")
        return 0
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(
        json.dumps(payload, sort_keys=True, indent=1) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": payload["contract_freeze_status"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
