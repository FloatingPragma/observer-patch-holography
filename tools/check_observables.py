#!/usr/bin/env python3
"""Validate the OPH completion ledger, docs/registers/OBSERVABLES.json.

The ledger lists every observable OPH must reproduce (kind "observable") and every
premise of its finite theorems (kind "premise", ids PR-nn). The status page
status.floatingpragma.ai reads the file directly. docs/registers/README.md documents
the format and the rules. This check keeps the file loadable and its references live:

- schema, key order, value types; one of the four statuses, defined in the file;
- short texts: `what` and `why` are one sentence each;
- unique ids; premise references resolve to premise rows;
- a Conditional observable names at least one premise the model does not demonstrate;
- a Demonstrated observable names no premise the model does not demonstrate;
- `lean` holds Lean files, `model` holds model/ and code/ files, `papers` holds the rest,
  and every path exists;
- every PR id that the claim registry or the scorecard's theorem criteria cite is a premise row;
- prose rules: no em-dash, no forms of "honest", no link to the private OpenWave repository.

It prints the counts and the completion score the status page shows.

Usage: python tools/check_observables.py
"""

from __future__ import annotations

import re
import sys
from collections import Counter
from pathlib import Path

if __package__:
    from .strict_json import load
else:
    from strict_json import load


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "docs" / "registers" / "OBSERVABLES.json"
CLAIMS = ROOT / "claims" / "claim_registry.yaml"
THEOREM_CRITERIA = ROOT / "model" / "evaluation" / "theorem_criteria.json"

SCHEMA = "oph.observables.v1"
STATUSES = ("Demonstrated", "Conditional", "Partial", "Open")
BASE_SCORE = {"Demonstrated": 1.0, "Conditional": 0.5, "Partial": 0.25, "Open": 0.0}
PREMISE_KINDS = {"declared", "derivable", "calibration", "mathematics"}
OBSERVABLE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9+._-]*$")
PREMISE_ID = re.compile(r"^PR-\d{2,3}$")
COMMON = ("id", "kind", "name", "category", "status", "what", "why", "premises", "lean", "model", "papers")
KEYS = {"observable": COMMON, "premise": COMMON + ("statement", "premise_kind")}
MAX_SHORT = 280
BANNED = (
    (re.compile("—"), "an em-dash"),
    (re.compile(r"\bhonest", re.I), "the word 'honest'"),
    (re.compile(r"github\.com/muellerberndt/openwave", re.I), "a link to the private OpenWave repository"),
)


def link_group(path: str) -> str:
    if path.endswith(".lean") or path.startswith("Lean/"):
        return "lean"
    if path.startswith(("model/", "code/")):
        return "model"
    return "papers"


def completion(rows: list[dict]) -> float:
    status = {r["id"]: r["status"] for r in rows}
    total = 0.0
    observables = [r for r in rows if r["kind"] == "observable"]
    for r in observables:
        score = BASE_SCORE[r["status"]]
        if r["status"] == "Conditional" and r["premises"]:
            score = 0.5 + 0.5 * sum(BASE_SCORE[status[p]] for p in r["premises"]) / len(r["premises"])
        total += score
    return total / max(len(observables), 1)


def main() -> int:
    errors: list[str] = []
    data = load(LEDGER)
    if data.get("schema") != SCHEMA:
        errors.append(f"schema must be {SCHEMA!r}")
    statuses = data.get("statuses", {})
    if list(statuses) != list(STATUSES) or not all(isinstance(v, str) and v for v in statuses.values()):
        errors.append(f"statuses must define {', '.join(STATUSES)} in this order with non-empty text")
    rows = data.get("rows")
    if not isinstance(rows, list) or not rows:
        print("check_observables: rows must be a non-empty list", file=sys.stderr)
        return 1

    if any(not isinstance(r, dict) or not isinstance(r.get("id"), str) for r in rows):
        print("check_observables: each row must be an object with a text id", file=sys.stderr)
        return 1
    ids = Counter(r["id"] for r in rows)
    errors += [f"duplicate id {i}" for i, n in ids.items() if n > 1]
    premise_status = {r["id"]: r.get("status") for r in rows if r.get("kind") == "premise" and "id" in r}

    for r in rows:
        rid = r.get("id", "<missing id>")
        kind = r.get("kind")
        if kind not in KEYS:
            errors.append(f"{rid}: kind must be 'observable' or 'premise'")
            continue
        if list(r) != list(KEYS[kind]):
            errors.append(f"{rid}: keys must be {list(KEYS[kind])} in this order")
            continue
        if not (OBSERVABLE_ID if kind == "observable" else PREMISE_ID).match(rid):
            errors.append(f"{rid}: malformed {kind} id")
        if r["status"] not in STATUSES:
            errors.append(f"{rid}: status {r['status']!r} is not one of {STATUSES}")
        for key in ("name", "category", "what", "why") + (("statement",) if kind == "premise" else ()):
            if not isinstance(r[key], str) or not r[key].strip():
                errors.append(f"{rid}: {key} must be non-empty text")
        for key in ("what", "why"):
            if isinstance(r[key], str) and len(r[key]) > MAX_SHORT:
                errors.append(f"{rid}: {key} is longer than {MAX_SHORT} characters; keep it to one sentence")
        if kind == "premise" and r["premise_kind"] not in PREMISE_KINDS:
            errors.append(f"{rid}: premise_kind must be one of {sorted(PREMISE_KINDS)}")
        valid_lists = True
        for key in ("premises", "lean", "model", "papers"):
            if not (isinstance(r[key], list) and all(isinstance(x, str) and x for x in r[key])):
                errors.append(f"{rid}: {key} must be a list of non-empty strings")
                valid_lists = False
            elif len(set(r[key])) != len(r[key]):
                errors.append(f"{rid}: {key} contains repeated entries")
        if not valid_lists:
            continue
        for p in r["premises"]:
            if p not in premise_status:
                errors.append(f"{rid}: premise {p} is not a premise row")
        if kind == "observable" and r["status"] == "Conditional":
            if not any(premise_status.get(p) not in (None, "Demonstrated") for p in r["premises"]):
                errors.append(f"{rid}: a Conditional row names at least one premise the model does not demonstrate")
        if kind == "observable" and r["status"] == "Demonstrated":
            waiting = [p for p in r["premises"] if premise_status.get(p) not in (None, "Demonstrated")]
            if waiting:
                errors.append(f"{rid}: a Demonstrated row lists premises the model does not demonstrate: {waiting}")
        for key in ("lean", "model", "papers"):
            for path in r[key]:
                if path.startswith(("http://", "https://")):
                    continue
                if link_group(path) != key:
                    errors.append(f"{rid}: {path} belongs in {link_group(path)}, not {key}")
                if not (ROOT / path.split("#", 1)[0]).is_file():
                    errors.append(f"{rid}: evidence path {path} must name a file")
        for key in ("name", "what", "why") + (("statement",) if kind == "premise" else ()):
            for pattern, label in BANNED:
                if pattern.search(r.get(key) or ""):
                    errors.append(f"{rid}: {key} contains {label}")

    cited: set[str] = set()
    for claim in load(CLAIMS)["claims"]:
        deps = claim.get("premise_dependencies") or {}
        for key in ("consumed", "open", "boundary"):
            cited.update(x for x in deps.get(key, []) if isinstance(x, str) and x.startswith("PR-"))
    if THEOREM_CRITERIA.exists():
        for group in load(THEOREM_CRITERIA)["groups"]:
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
    print(f"completion: {100 * completion(rows):.1f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
