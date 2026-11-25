#!/usr/bin/env python3
"""Validate the OPH completion ledger, docs/registers/OBSERVABLES.json.

The ledger lists every observable OPH must reproduce (kind "observable") and every
premise of its finite theorems (kind "premise", ids PR-nn). The public status page
reads the file directly, so this check keeps it loadable and its references live:

- schema, required keys and value types per row kind;
- one of the four statuses, whose definitions the file carries;
- unique ids; premise references resolve to premise rows;
- a Conditional observable names at least one premise the model does not demonstrate;
- every cited paper, code or evidence path exists in the repository;
- every claim id exists in claims/claim_registry.yaml;
- every PR id that the claim registry or the scorecard's theorem criteria cite is a
  premise row;
- prose rules: no em-dash, no forms of "honest", no link to the private OpenWave repository.

Usage: python tools/check_observables.py
"""

from __future__ import annotations

import json
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "docs" / "registers" / "OBSERVABLES.json"
CLAIMS = ROOT / "claims" / "claim_registry.yaml"
THEOREM_CRITERIA = ROOT / "model" / "evaluation" / "theorem_criteria.json"

SCHEMA = "oph.observables.v1"
STATUSES = ("Demonstrated", "Conditional", "Partial", "Open")
SCALES = {"Micro", "Quantum", "Spacetime", "Matter", "Cosmic"}
PREMISE_KINDS = {"declared", "derivable", "calibration", "mathematics"}
OBSERVABLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9+._-]*$")
PREMISE_ID = re.compile(r"^PR-\d{2,3}$")
COMMON = ("id", "kind", "name", "category", "status", "summary", "premises", "claims", "papers", "code", "tags")
OBSERVABLE_KEYS = COMMON + ("scale", "why", "value", "observed")
PREMISE_KEYS = COMMON + ("statement", "premise_kind", "used_by")
LISTS = ("premises", "claims", "papers", "code", "tags", "used_by")
PROSE = ("name", "summary", "why", "statement")
BANNED = (
    (re.compile("—"), "em-dash"),
    (re.compile(r"\bhonest", re.I), "the word 'honest'"),
    (re.compile(r"github\.com/muellerberndt/openwave", re.I), "a link to the private OpenWave repository"),
)


def is_path(value: str) -> bool:
    return "/" in value and " " not in value and not value.startswith(("http://", "https://"))


def main() -> int:
    errors: list[str] = []
    data = json.loads(LEDGER.read_text(encoding="utf-8"))

    if data.get("schema") != SCHEMA:
        errors.append(f"schema must be {SCHEMA!r}")
    statuses = data.get("statuses", {})
    if set(statuses) != set(STATUSES) or not all(isinstance(v, str) and v for v in statuses.values()):
        errors.append(f"statuses must define exactly {', '.join(STATUSES)} with non-empty text")
    rows = data.get("rows")
    if not isinstance(rows, list) or not rows:
        print("check_observables: rows must be a non-empty list", file=sys.stderr)
        return 1

    claim_ids = {c["claim_id"] for c in json.loads(CLAIMS.read_text(encoding="utf-8"))["claims"]}
    ids: Counter[str] = Counter(r.get("id") for r in rows)
    errors += [f"duplicate id {i}" for i, n in ids.items() if n > 1]
    premise_status = {r["id"]: r.get("status") for r in rows if r.get("kind") == "premise" and "id" in r}

    for r in rows:
        rid = r.get("id", "<missing id>")
        kind = r.get("kind")
        if kind == "observable":
            keys, pattern = OBSERVABLE_KEYS, OBSERVABLE_ID
        elif kind == "premise":
            keys, pattern = PREMISE_KEYS, PREMISE_ID
        else:
            errors.append(f"{rid}: kind must be 'observable' or 'premise'")
            continue
        if list(r) != list(keys):
            errors.append(f"{rid}: keys must be {list(keys)} in this order")
            continue
        if not pattern.match(rid):
            errors.append(f"{rid}: malformed {kind} id")
        if r["status"] not in STATUSES:
            errors.append(f"{rid}: status {r['status']!r} is not one of {STATUSES}")
        for key in ("name", "category", "summary") + (("statement",) if kind == "premise" else ("why",)):
            if not isinstance(r[key], str) or not r[key].strip():
                errors.append(f"{rid}: {key} must be non-empty text")
        for key in LISTS:
            if key in r and not (isinstance(r[key], list) and all(isinstance(x, str) and x for x in r[key])):
                errors.append(f"{rid}: {key} must be a list of non-empty strings")
        if kind == "observable":
            if r["scale"] not in SCALES:
                errors.append(f"{rid}: scale must be one of {sorted(SCALES)}")
            for key in ("value", "observed"):
                if r[key] is not None and not isinstance(r[key], str):
                    errors.append(f"{rid}: {key} must be text or null")
        else:
            if r["premise_kind"] not in PREMISE_KINDS:
                errors.append(f"{rid}: premise_kind must be one of {sorted(PREMISE_KINDS)}")
        for p in r["premises"]:
            if p not in premise_status:
                errors.append(f"{rid}: premise {p} is not a premise row")
        if kind == "observable" and r["status"] == "Conditional":
            if not any(premise_status.get(p) not in (None, "Demonstrated") for p in r["premises"]):
                errors.append(f"{rid}: a Conditional row names at least one premise the model does not demonstrate")
        for c in r["claims"]:
            if c not in claim_ids:
                errors.append(f"{rid}: claim {c} is not in claims/claim_registry.yaml")
        paths = r["papers"] + r["code"] + [u for u in r.get("used_by", []) if is_path(u)]
        for path in paths:
            if path.startswith(("http://", "https://")):
                continue
            if not (ROOT / path.split("#", 1)[0]).exists():
                errors.append(f"{rid}: path {path} does not exist")
        for key in PROSE:
            text = r.get(key) or ""
            for pattern_, label in BANNED:
                if pattern_.search(text):
                    errors.append(f"{rid}: {key} contains {label}")

    cited: set[str] = set()
    for claim in json.loads(CLAIMS.read_text(encoding="utf-8"))["claims"]:
        deps = claim.get("premise_dependencies") or {}
        for key in ("consumed", "open", "boundary"):
            cited.update(x for x in deps.get(key, []) if isinstance(x, str) and x.startswith("PR-"))
    if THEOREM_CRITERIA.exists():
        for group in json.loads(THEOREM_CRITERIA.read_text(encoding="utf-8"))["groups"]:
            for criterion in group["criteria"]:
                cited.update(criterion.get("premises", []))
    for p in sorted(cited - set(premise_status)):
        errors.append(f"premise {p} is cited by the claim registry or the theorem criteria but has no row")

    if errors:
        for e in errors:
            print(f"check_observables: {e}", file=sys.stderr)
        print(f"check_observables: {len(errors)} problem(s)", file=sys.stderr)
        return 1
    for kind in ("observable", "premise"):
        tally = Counter(r["status"] for r in rows if r["kind"] == kind)
        print(f"{kind}s: {sum(tally.values())} ({', '.join(f'{s} {tally[s]}' for s in STATUSES)})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
