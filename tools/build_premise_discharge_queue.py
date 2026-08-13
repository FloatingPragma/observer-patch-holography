#!/usr/bin/env python3
"""Build and validate the V3 premise-discharge queue (issue #739).

The queue is an immutable-origin worklist, not a second premise register.  Its
AV-0 statements and evidence are recovered from the anchored AV-0 Git commit;
the current register must still equal that origin inventory.  The architecture
version currently being evaluated is recorded separately, so moving to AV-1
cannot silently relabel the origin of an old obligation.

Reverse edges are reconstructed exactly from the canonical premise register,
observation ledger, claim registry, architecture decisions, and the explicit
non-observation surface map maintained by the observation-ledger validator.
Every input is content-pinned and every evidence path is pinned at the AV-0
origin commit.  Any mismatch fails closed.
"""

from __future__ import annotations

import argparse
import copy
import functools
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

import build_architecture_versions as architecture_versions
import build_observation_ledger as observation_ledger
import build_premise_register as premise_register
import strict_json


ROOT = Path(__file__).resolve().parents[1]
PREMISE_REGISTER_PATH = ROOT / "tracking" / "premise_register.json"
ARCHITECTURE_REGISTER_PATH = ROOT / "tracking" / "architecture_versions.json"
ISSUE_SNAPSHOT_PATH = ROOT / "tracking" / "open_issues" / "open_problem_ledger.json"
OBSERVATION_LEDGER_PATH = ROOT / "tracking" / "observation_ledger.json"
CLAIM_REGISTRY_PATH = ROOT / "claims" / "claim_registry.yaml"
OBSERVATION_LOADER_PATH = ROOT / "tools" / "build_observation_ledger.py"
QUEUE_PATH = ROOT / "tracking" / "premise_discharge_queue.json"
SURFACE_PATH = ROOT / "docs" / "PREMISE_DISCHARGE_QUEUE_V3.md"

SCHEMA = "oph.premise_discharge_queue.v2"
ISSUE = 739
AUDIT_ISSUE = 738
ARCHITECTURE_ISSUE = 741
ORIGIN_ARCHITECTURE_VERSION = "AV-0"
ISSUE_URL = "https://github.com/FloatingPragma/observer-patch-holography/issues"
QUEUED_DISPOSITIONS = ("remove", "axiomatize")
STATE = "deferred_open"
VERSION_RE = re.compile(r"^AV-(0|[1-9][0-9]*)$")
REVISION_RE = re.compile(r"^[0-9a-f]{40}$")
DIGEST_RE = re.compile(r"^[0-9a-f]{64}$")
OBJECT_RE = re.compile(r"^[0-9a-f]{40}$")
PREMISE_TOKEN_RE = re.compile(r"\bPR-[0-9]{2}\b")

TOP_KEYS = {
    "schema",
    "issue",
    "audit_custody_issue",
    "decision_custody",
    "policy",
    "input_custody",
    "reverse_edge_coverage",
    "summary",
    "path_contracts",
    "inherited_v2_residuals",
    "excluded_imports",
    "items",
}
ITEM_KEYS = {
    "triage_rank",
    "queue_id",
    "premise_id",
    "premise_name",
    "premise_type",
    "disposition",
    "state",
    "origin_architecture_version",
    "evaluation_architecture_version",
    "origin_premise_record_sha256",
    "fanout_band",
    "consumer_edges",
    "premise_statement",
    "next_action_or_decision",
    "origin_evidence",
    "available_exit_paths",
    "decision_custody",
}
CONSUMER_KEYS = {
    "lane_issues",
    "observation_rows",
    "claim_ids",
    "architecture_decisions",
    "non_observation_surfaces",
    "evidence_references",
    "unreferenced_evidence_paths",
}
CONTENT_PIN_KEYS = {"path", "byte_count", "sha256", "git_blob_sha1"}
EVIDENCE_PIN_KEYS = {
    "path",
    "origin_revision",
    "git_object_type",
    "git_object_sha1",
    "content_sha256",
    "byte_count",
    "descendant_count",
}

POLICY = (
    "Every remove or axiomatize premise has one visible deferred queue item. "
    "Generation is not discharge: a bounded no-go remains evidence and an open "
    "item. A derivation or axiomatization is recorded in the immutable "
    "architecture/audit registers established by issues #741 and #738 and is "
    "followed by replay of every exact reverse consumer. The issue may close "
    "after bootstrapping that durable custody; its live state is not decision "
    "semantics. Named-premise adequacy may proceed while an item is open; "
    "axioms-only necessity may not claim that item closed."
)

DECISION_CUSTODY = {
    "architecture_register": "tracking/architecture_versions.json",
    "audit_register": "tracking/audit_custody.json",
    "process": "append an immutable successor architecture version, replay every exact registered reverse edge, and attach an independent audit record before promotion",
    "bootstrap_issue": ARCHITECTURE_ISSUE,
}

PATH_CONTRACTS = {
    "derive": {
        "owner_issue": ISSUE,
        "closure_requirements": [
            "A target-clean source construction proves the exact registered premise statement from the declared architecture inputs.",
            "A committed producer, independent verifier, semantic mutation controls, and independent audit pin the result and its scope.",
            "The architecture decision register records the derivation against an immutable architecture version and enumerates every affected consumer replay.",
            "The premise register and all affected generated surfaces are updated together; no conditional claim is promoted merely because the queue item closes.",
        ],
        "non_closing_outputs": [
            "a bounded no-go without a positive source derivation",
            "a side-by-side theorem re-export without one typed composition",
            "a post-hoc or target-informed fit",
        ],
    },
    "axiomatize": {
        "owner_registry": "tracking/architecture_versions.json",
        "audit_registry": "tracking/audit_custody.json",
        "bootstrap_issue": ARCHITECTURE_ISSUE,
        "eligible_disposition": "axiomatize",
        "closure_requirements": [
            "The decision record states the exact architectural clause, rationale, alternatives, and retained source-selection no-gos.",
            "A basis-wide dependency, countermodel, and empirical-input audit identifies every affected premise and observation row.",
            "A successor architecture version freezes the changed normative bytes; no historical version is rewritten in place.",
            "Every conditional consumer is replayed before any promotion, with underlying audit artifacts retained under standing issue #738.",
        ],
        "non_closing_outputs": [
            "editing a premise statement without a version event",
            "calling a declared premise derived",
            "treating a version identifier as evidence of physical emergence",
        ],
    },
}

# Historical custody explicitly absorbed by the live #739 contract.  This is
# queue-level provenance, not an invented one-to-one attribution to newer
# atomic premises.
INHERITED_V2_RESIDUALS = [
    {"code": "B7", "issue": 683, "residual": "reference and real-enrichment selection"},
    {"code": "B9", "issue": 685, "residual": "spectral-information acceptance"},
    {"code": "E1", "issue": 692, "residual": "nonconstant source realization"},
    {"code": "B13", "issue": 702, "residual": "affinity-principle discharge"},
    {"code": "E5", "issue": 703, "residual": "clock derivation"},
    {"code": "B14", "issue": 705, "residual": "selection-rule derivation"},
    {"code": "B15", "issue": 706, "residual": "source matter selection"},
    {"code": "B16", "issue": 707, "residual": "character-category and global-form identities"},
    {"code": "E9", "issue": 716, "residual": "bracket from dynamics"},
    {"code": "B20", "issue": 725, "residual": "enriched-export decision"},
]


def fail(message: str) -> None:
    raise SystemExit(f"premise discharge queue: {message}")


def _display_path(path: Path) -> str:
    try:
        return path.relative_to(ROOT).as_posix()
    except ValueError:
        return str(path)


def load_json(path: Path) -> dict:
    try:
        return strict_json.load(path)
    except FileNotFoundError:
        fail(f"missing input {_display_path(path)}")
    except (json.JSONDecodeError, strict_json.DuplicateKeyError) as error:
        fail(f"invalid JSON in {_display_path(path)}: {error}")
    raise AssertionError("unreachable")


def _git_output(*args: str) -> bytes:
    result = subprocess.run(
        ["git", *args],
        cwd=ROOT,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode:
        detail = result.stderr.decode("utf-8", errors="replace").strip()
        fail(f"git {' '.join(args)} failed: {detail}")
    return result.stdout


def _git_blob_sha1(payload: bytes) -> str:
    header = f"blob {len(payload)}\0".encode("ascii")
    return hashlib.sha1(header + payload).hexdigest()


def _content_pin(path: Path) -> dict:
    if not path.is_file():
        fail(f"custody input must be a file: {_display_path(path)}")
    payload = path.read_bytes()
    return {
        "path": path.relative_to(ROOT).as_posix(),
        "byte_count": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "git_blob_sha1": _git_blob_sha1(payload),
    }


def _canonical_sha256(value: object) -> str:
    payload = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":")
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _strict_json_at_revision(revision: str, path: Path) -> dict:
    relative = path.relative_to(ROOT).as_posix()
    payload = _git_output("show", f"{revision}:{relative}")
    try:
        value = strict_json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError, strict_json.DuplicateKeyError) as error:
        fail(f"invalid pinned JSON at {revision}:{relative}: {error}")
    if not isinstance(value, dict):
        fail(f"pinned JSON at {revision}:{relative} must be an object")
    return value


@functools.lru_cache(maxsize=None)
def _evidence_pin(revision: str, relative: str) -> dict:
    if not REVISION_RE.fullmatch(revision):
        fail("evidence origin must be a full Git commit")
    path = ROOT / relative
    if path.is_absolute() and (".." in Path(relative).parts or relative.startswith("/")):
        fail(f"evidence path must be repository-relative: {relative}")
    if not path.exists():
        fail(f"current evidence path does not resolve: {relative}")

    object_spec = f"{revision}:{relative.rstrip('/')}"
    object_sha1 = _git_output("rev-parse", object_spec).decode("ascii").strip()
    object_type = _git_output("cat-file", "-t", object_sha1).decode("ascii").strip()
    if not OBJECT_RE.fullmatch(object_sha1) or object_type not in {"blob", "tree"}:
        fail(f"unsupported pinned evidence object {object_spec}")

    if object_type == "blob":
        payload = _git_output("show", object_spec)
        descendant_count = 1
    else:
        # A tree pin includes a canonical recursive manifest, so changing any
        # descendant changes both the tree object and the SHA-256 below.
        payload = _git_output("ls-tree", "-r", "-z", revision, "--", relative.rstrip("/"))
        descendant_count = payload.count(b"\0")
        if descendant_count == 0:
            fail(f"pinned evidence tree is empty: {object_spec}")

    return {
        "path": relative,
        "origin_revision": revision,
        "git_object_type": object_type,
        "git_object_sha1": object_sha1,
        "content_sha256": hashlib.sha256(payload).hexdigest(),
        "byte_count": len(payload),
        "descendant_count": descendant_count,
    }


def _issue_numbers(snapshot: dict) -> set[int]:
    rows = snapshot.get("rows") if isinstance(snapshot, dict) else None
    if not isinstance(rows, list):
        fail("open-issue snapshot rows must be a list")
    numbers: list[int] = []
    for row in rows:
        number = row.get("number") if isinstance(row, dict) else None
        if not isinstance(number, int) or isinstance(number, bool):
            fail("open-issue snapshot contains a malformed issue number")
        numbers.append(number)
    if len(numbers) != len(set(numbers)):
        fail("open-issue snapshot contains duplicate issue numbers")
    if snapshot.get("open_issue_count") != len(rows):
        fail("open-issue snapshot count does not match its row inventory")
    # #741 may close after its durable register exists.  Its live issue state
    # is deliberately not used as architecture decision semantics.
    for required in (AUDIT_ISSUE, ISSUE):
        if required not in numbers:
            fail(f"required live owner issue #{required} is absent from the snapshot")
    return set(numbers)


def _claim_rows(registry: dict, known_premises: set[str]) -> list[dict]:
    claims = registry.get("claims") if isinstance(registry, dict) else None
    if not isinstance(claims, list) or not claims:
        fail("claim registry claims must be a nonempty list")
    seen: set[str] = set()
    for claim in claims:
        claim_id = claim.get("claim_id") if isinstance(claim, dict) else None
        evidence = claim.get("evidence") if isinstance(claim, dict) else None
        if not isinstance(claim_id, str) or not claim_id or claim_id in seen:
            fail("claim registry contains a malformed or duplicate claim id")
        if not isinstance(evidence, list) or not all(
            isinstance(path, str) and path for path in evidence
        ):
            fail(f"{claim_id}: evidence must be a list of paths")
        unknown = sorted(
            set(PREMISE_TOKEN_RE.findall(json.dumps(claim, sort_keys=True)))
            - known_premises
        )
        if unknown:
            fail(f"{claim_id}: references unknown premise ids {unknown}")
        seen.add(claim_id)
    return claims


def _origin_anchor(architecture: dict) -> tuple[dict, dict]:
    anchors = architecture["version_anchors"]
    versions = architecture["versions"]
    anchor = next(
        (row for row in anchors if row["id"] == ORIGIN_ARCHITECTURE_VERSION), None
    )
    version = next(
        (row for row in versions if row["id"] == ORIGIN_ARCHITECTURE_VERSION), None
    )
    if anchor is None or version is None:
        fail("architecture register does not anchor AV-0")
    return anchor, version


def load_sources() -> dict:
    current_register = premise_register.load_json(PREMISE_REGISTER_PATH)
    premise_register.validate(current_register)

    architecture = architecture_versions.load_json(ARCHITECTURE_REGISTER_PATH)
    architecture_versions.validate(architecture)
    current_version = architecture.get("current_version")
    if not isinstance(current_version, str) or not VERSION_RE.fullmatch(current_version):
        fail("architecture register has a malformed current_version")
    origin_anchor, origin_version = _origin_anchor(architecture)
    origin_revision = origin_anchor["origin_revision"]

    origin_register = _strict_json_at_revision(origin_revision, PREMISE_REGISTER_PATH)
    origin_rows = premise_register.validate(origin_register)
    if current_register != origin_register:
        fail(
            "live premise register differs from the anchored AV-0 origin; "
            "append an explicit queue decision/event schema before changing an origin row"
        )

    issue_snapshot = load_json(ISSUE_SNAPSHOT_PATH)
    open_issues = _issue_numbers(issue_snapshot)

    observations = load_json(OBSERVATION_LEDGER_PATH)
    observation_rows = observation_ledger.validate(observations)
    claims_registry = load_json(CLAIM_REGISTRY_PATH)
    known_premises = {row["id"] for row in origin_rows}
    claims = _claim_rows(claims_registry, known_premises)

    return {
        "rows": origin_rows,
        "current_version": current_version,
        "open_issues": open_issues,
        "architecture": architecture,
        "origin_anchor": origin_anchor,
        "origin_version": origin_version,
        "origin_revision": origin_revision,
        "observation_rows": observation_rows,
        "claims": claims,
        "input_custody": {
            "premise_register": {
                **_content_pin(PREMISE_REGISTER_PATH),
                "origin_revision": origin_revision,
                "origin_git_blob_sha1": _git_output(
                    "rev-parse",
                    f"{origin_revision}:{PREMISE_REGISTER_PATH.relative_to(ROOT).as_posix()}",
                )
                .decode("ascii")
                .strip(),
            },
            "architecture_register": _content_pin(ARCHITECTURE_REGISTER_PATH),
            "origin_architecture_snapshot": {
                "id": ORIGIN_ARCHITECTURE_VERSION,
                "origin_revision": origin_revision,
                "snapshot_sha256": origin_version["snapshot_sha256"],
                "record_sha256": origin_anchor["record_sha256"],
            },
            "open_issue_snapshot": _content_pin(ISSUE_SNAPSHOT_PATH),
            "reverse_edge_inputs": {
                "observation_ledger": _content_pin(OBSERVATION_LEDGER_PATH),
                "claim_registry": _content_pin(CLAIM_REGISTRY_PATH),
                "non_observation_surface_map": _content_pin(OBSERVATION_LOADER_PATH),
            },
        },
    }


def _fanout_band(consumer_count: int) -> str:
    if consumer_count >= 4:
        return "high_fanout"
    if consumer_count >= 2:
        return "cross_lane"
    return "single_lane"


def _premise_number(row: dict) -> int:
    return int(row["id"].split("-", 1)[1])


def _queue_rows(rows: list[dict]) -> list[dict]:
    return sorted(
        (row for row in rows if row["disposition"] in QUEUED_DISPOSITIONS),
        key=lambda row: (-len(row["consuming_lanes"]), _premise_number(row)),
    )


def _consumer_edges(row: dict, sources: dict) -> dict:
    premise_id = row["id"]
    observations = []
    for observation in sources["observation_rows"]:
        if premise_id in observation["premises"]:
            observations.append({"id": observation["id"], "role": "consumed_premise"})
        if premise_id in observation["open_premises"]:
            observations.append({"id": observation["id"], "role": "open_premise"})
    observations.sort(key=lambda edge: (edge["id"], edge["role"]))

    explicit_claims: list[str] = []
    for claim in sources["claims"]:
        tokens = set(PREMISE_TOKEN_RE.findall(json.dumps(claim, sort_keys=True)))
        if premise_id in tokens:
            explicit_claims.append(claim["claim_id"])
    explicit_claims.sort()

    architecture_decisions = sorted(
        (
            {"architecture_version": version["id"], "decision_id": decision["id"]}
            for version in sources["architecture"]["versions"]
            for decision in version["protocol_decisions"]
            if premise_id in decision["premises"]
        ),
        key=lambda edge: (edge["architecture_version"], edge["decision_id"]),
    )

    non_observation = sorted(
        (
            {"lane_issue": lane, "surface": surface}
            for (lane, mapped_premise), surface in observation_ledger.NON_OBSERVATION_SURFACE_CONSUMERS.items()
            if mapped_premise == premise_id
        ),
        key=lambda edge: (edge["lane_issue"], edge["surface"]),
    )

    evidence_references = []
    unreferenced_evidence_paths = []
    for path in row["evidence"]:
        observation_ids = sorted(
            observation["id"]
            for observation in sources["observation_rows"]
            if path in observation["evidence"]
        )
        claim_ids = sorted(
            claim["claim_id"] for claim in sources["claims"] if path in claim["evidence"]
        )
        if observation_ids or claim_ids:
            evidence_references.append(
                {
                    "path": path,
                    "observation_row_ids": observation_ids,
                    "claim_ids": claim_ids,
                    "role": "shared_evidence_reference_not_semantic_premise_inference",
                }
            )
        else:
            unreferenced_evidence_paths.append(path)

    return {
        "lane_issues": list(row["consuming_lanes"]),
        "observation_rows": observations,
        "claim_ids": explicit_claims,
        "architecture_decisions": architecture_decisions,
        "non_observation_surfaces": non_observation,
        "evidence_references": evidence_references,
        "unreferenced_evidence_paths": unreferenced_evidence_paths,
    }


def _reverse_edge_coverage(items: list[dict], sources: dict) -> dict:
    unannotated_claims = sorted(
        claim["claim_id"]
        for claim in sources["claims"]
        if not PREMISE_TOKEN_RE.findall(json.dumps(claim, sort_keys=True))
    )
    unreferenced_paths = [
        {
            "premise_id": item["premise_id"],
            "paths": item["consumer_edges"]["unreferenced_evidence_paths"],
        }
        for item in items
        if item["consumer_edges"]["unreferenced_evidence_paths"]
    ]
    return {
        "status": "exact_registered_edges_only__claim_semantic_mapping_not_exhaustive",
        "observation_edge_method": "exact membership in observation-ledger premises or open_premises",
        "claim_edge_method": "exact explicit PR-xx token in one canonical claim record; shared evidence is reported separately and never promoted to a semantic edge",
        "surface_edge_method": "exact premise/lane pair in NON_OBSERVATION_SURFACE_CONSUMERS",
        "claim_records_without_explicit_premise_ids": unannotated_claims,
        "evidence_paths_without_registered_reverse_reference": unreferenced_paths,
        "closure_guard": "No queue item may claim claim-exhaustive replay from this queue alone until every potentially affected claim has an explicit premise edge or an audited non-consumer certificate. An evidence path without a canonical reverse reference is reported, never assigned an invented consumer.",
    }


def build_queue(sources: dict) -> dict:
    if not {AUDIT_ISSUE, ISSUE} <= sources["open_issues"]:
        fail("queue and standing-audit owners must both be live")

    queued = _queue_rows(sources["rows"])
    items: list[dict] = []
    for rank, row in enumerate(queued, start=1):
        paths = ["derive"]
        if row["disposition"] == "axiomatize":
            paths.append("axiomatize")
        items.append(
            {
                "triage_rank": rank,
                "queue_id": f"D1-{row['id']}",
                "premise_id": row["id"],
                "premise_name": row["name"],
                "premise_type": row["type"],
                "disposition": row["disposition"],
                "state": STATE,
                "origin_architecture_version": ORIGIN_ARCHITECTURE_VERSION,
                "evaluation_architecture_version": sources["current_version"],
                "origin_premise_record_sha256": _canonical_sha256(row),
                "fanout_band": _fanout_band(len(row["consuming_lanes"])),
                "consumer_edges": _consumer_edges(row, sources),
                "premise_statement": row["statement"],
                "next_action_or_decision": row["notes"],
                "origin_evidence": [
                    copy.deepcopy(_evidence_pin(sources["origin_revision"], path))
                    for path in row["evidence"]
                ],
                "available_exit_paths": paths,
                "decision_custody": DECISION_CUSTODY,
            }
        )

    imports = [
        {
            "premise_id": row["id"],
            "premise_name": row["name"],
            "premise_type": row["type"],
            "disposition": row["disposition"],
            "reason": "Flagged import: outside #739 unless a recorded disposition change creates a new queue item.",
        }
        for row in sources["rows"]
        if row["disposition"] == "import"
    ]

    remove_count = sum(row["disposition"] == "remove" for row in queued)
    axiomatize_count = sum(row["disposition"] == "axiomatize" for row in queued)
    return {
        "schema": SCHEMA,
        "issue": ISSUE,
        "audit_custody_issue": AUDIT_ISSUE,
        "decision_custody": DECISION_CUSTODY,
        "policy": POLICY,
        "input_custody": sources["input_custody"],
        "reverse_edge_coverage": _reverse_edge_coverage(items, sources),
        "summary": {
            "registered_premises": len(sources["rows"]),
            "queued_items": len(queued),
            "remove_items": remove_count,
            "axiomatize_items": axiomatize_count,
            "excluded_imports": len(imports),
            "state": STATE,
            "triage_rule": "descending lane-consumer count, then ascending premise id; routing aid only, not an effort or truth score",
        },
        "path_contracts": PATH_CONTRACTS,
        "inherited_v2_residuals": INHERITED_V2_RESIDUALS,
        "excluded_imports": imports,
        "items": items,
    }


def _validate_content_pin(where: str, pin: object) -> None:
    if not isinstance(pin, dict) or not CONTENT_PIN_KEYS <= set(pin):
        fail(f"{where}: malformed content pin")
    if not isinstance(pin["byte_count"], int) or pin["byte_count"] < 0:
        fail(f"{where}: malformed byte count")
    if not DIGEST_RE.fullmatch(pin["sha256"]):
        fail(f"{where}: malformed SHA-256")
    if not OBJECT_RE.fullmatch(pin["git_blob_sha1"]):
        fail(f"{where}: malformed Git blob id")


def _validate_evidence_pin(where: str, pin: object) -> None:
    if not isinstance(pin, dict) or set(pin) != EVIDENCE_PIN_KEYS:
        fail(f"{where}: malformed origin evidence pin")
    if not REVISION_RE.fullmatch(pin["origin_revision"]):
        fail(f"{where}: malformed origin revision")
    if pin["git_object_type"] not in {"blob", "tree"}:
        fail(f"{where}: malformed Git object type")
    if not OBJECT_RE.fullmatch(pin["git_object_sha1"]):
        fail(f"{where}: malformed Git object id")
    if not DIGEST_RE.fullmatch(pin["content_sha256"]):
        fail(f"{where}: malformed content SHA-256")
    if not isinstance(pin["byte_count"], int) or pin["byte_count"] < 0:
        fail(f"{where}: malformed byte count")
    if not isinstance(pin["descendant_count"], int) or pin["descendant_count"] < 1:
        fail(f"{where}: malformed descendant count")


def validate_queue(queue: dict, sources: dict) -> list[dict]:
    if not isinstance(queue, dict) or set(queue) != TOP_KEYS:
        fail(f"top-level keys must equal {sorted(TOP_KEYS)}")
    expected = build_queue(sources)
    for key in TOP_KEYS - {"items"}:
        if queue[key] != expected[key]:
            fail(f"{key} has drifted from the canonical pinned inputs")

    custody = queue["input_custody"]
    for key in ("premise_register", "architecture_register", "open_issue_snapshot"):
        _validate_content_pin(f"input_custody.{key}", custody[key])
    for key, pin in custody["reverse_edge_inputs"].items():
        _validate_content_pin(f"input_custody.reverse_edge_inputs.{key}", pin)

    items = queue["items"]
    expected_items = expected["items"]
    if not isinstance(items, list) or len(items) != len(expected_items):
        fail(f"items must contain exactly {len(expected_items)} queue rows")
    for index, (item, expected_item) in enumerate(zip(items, expected_items), start=1):
        where = f"item {index}"
        if not isinstance(item, dict) or set(item) != ITEM_KEYS:
            fail(f"{where}: keys must equal {sorted(ITEM_KEYS)}")
        if not isinstance(item["consumer_edges"], dict) or set(item["consumer_edges"]) != CONSUMER_KEYS:
            fail(f"{where}: consumer_edges keys must equal {sorted(CONSUMER_KEYS)}")
        for evidence_index, pin in enumerate(item["origin_evidence"], start=1):
            _validate_evidence_pin(f"{where} evidence {evidence_index}", pin)
        for key in ITEM_KEYS:
            if item[key] != expected_item[key]:
                fail(
                    f"{where} ({expected_item['queue_id']}): {key} has drifted "
                    "from the canonical pinned inputs"
                )
    return items


def _issue_links(issues: list[int]) -> str:
    return ", ".join(f"[#{issue}]({ISSUE_URL}/{issue})" for issue in issues) or "none"


def _paths(paths: list[str]) -> str:
    return ", ".join(f"`{path}`" for path in paths)


def _table_text(value: str) -> str:
    return " ".join(value.split()).replace("|", "\\|")


def _edge_summary(edges: dict) -> str:
    observations = ", ".join(edge["id"] for edge in edges["observation_rows"]) or "none"
    claims = ", ".join(edges["claim_ids"]) or "none"
    surfaces = ", ".join(
        f"#{edge['lane_issue']} {edge['surface']}" for edge in edges["non_observation_surfaces"]
    ) or "none"
    decisions = ", ".join(
        f"{edge['architecture_version']}/{edge['decision_id']}"
        for edge in edges["architecture_decisions"]
    ) or "none"
    return (
        f"lanes {_issue_links(edges['lane_issues'])}; observations {observations}; "
        f"claims {claims}; architecture decisions {decisions}; other surfaces {surfaces}"
    )


def render(queue: dict) -> str:
    summary = queue["summary"]
    custody = queue["input_custody"]
    origin = custody["origin_architecture_snapshot"]
    premise_pin = custody["premise_register"]
    lines = [
        "# OPH V3 Premise Discharge Queue",
        "",
        "Generated by `tools/build_premise_discharge_queue.py`. The AV-0 "
        "premise rows and evidence are recovered from their immutable origin "
        "commit; the observation ledger, claim registry, architecture register, "
        "surface map, and issue snapshot are content-pinned reverse-edge inputs. "
        f"[Issue #{ISSUE}]({ISSUE_URL}/{ISSUE}) owns derivation work. The "
        f"register bootstraps established by [#{ARCHITECTURE_ISSUE}]({ISSUE_URL}/{ARCHITECTURE_ISSUE}) "
        f"and [#{AUDIT_ISSUE}]({ISSUE_URL}/{AUDIT_ISSUE}) retain decisions and audits even after a bootstrap issue closes.",
        "",
        queue["policy"],
        "",
        f"Origin architecture: `{origin['id']}` at `{origin['origin_revision']}`. "
        f"Current evaluation architecture: `{queue['items'][0]['evaluation_architecture_version']}`. "
        f"Premise-register SHA-256: `{premise_pin['sha256']}`; Git blob: "
        f"`{premise_pin['git_blob_sha1']}`; bytes: `{premise_pin['byte_count']}`.",
        "",
        "For origin evidence, the SHA-256 and byte count cover the file bytes "
        "for a Git blob and the canonical recursive Git manifest for a tree; "
        "the full origin commit and Git object ID are recorded in both cases.",
        "",
        "An origin version never changes when a later architecture is evaluated. "
        "This schema fails closed if an AV-0 premise row changes; a future "
        "discharge must add an explicit decision/event schema before changing the origin inventory.",
        "",
        "## Reverse-edge coverage boundary",
        "",
        queue["reverse_edge_coverage"]["closure_guard"],
        "",
        f"- Claim records without an explicit premise ID: {len(queue['reverse_edge_coverage']['claim_records_without_explicit_premise_ids'])}.",
        f"- Premise rows with one or more evidence paths lacking a canonical reverse reference: {len(queue['reverse_edge_coverage']['evidence_paths_without_registered_reverse_reference'])}.",
        "- These are reported custody gaps. The generator does not invent consumers from thematic similarity or a shared evidence path.",
        "",
        "## Queue census",
        "",
        f"- Registered premises: {summary['registered_premises']}",
        f"- Queued items: {summary['queued_items']} ({summary['remove_items']} `remove`, {summary['axiomatize_items']} `axiomatize`)",
        f"- Flagged imports outside the queue: {summary['excluded_imports']}",
        f"- Current queue state: `{summary['state']}`",
        f"- Triage order: {summary['triage_rule']}.",
        "",
        "## Actionable queue",
        "",
        "| Rank | Item | Premise | Disposition | Fanout | Lane consumers | Observation rows | Explicit claim IDs | Allowed exits | Next action or decision |",
        "| ---: | --- | --- | --- | --- | --- | --- | --- | --- | --- |",
    ]
    for item in queue["items"]:
        edges = item["consumer_edges"]
        observations = ", ".join(edge["id"] for edge in edges["observation_rows"]) or "none"
        claims = ", ".join(edges["claim_ids"]) or "none"
        lines.append(
            f"| {item['triage_rank']} | `{item['queue_id']}` | {item['premise_id']} {item['premise_name']} | "
            f"`{item['disposition']}` | `{item['fanout_band']}` | {_issue_links(edges['lane_issues'])} | "
            f"{observations} | {claims} | {_paths(item['available_exit_paths'])} | "
            f"{_table_text(item['next_action_or_decision'])} |"
        )

    lines += [
        "",
        "Claim IDs above are exact explicit `PR-xx` references in the canonical "
        "claim records; shared evidence paths are reported separately and never "
        "invent a semantic premise edge. The explicit non-observation surface "
        "map covers lane consumers not attached to an observation row.",
        "",
        "The rank is a deterministic lane-fanout routing aid. It is not a scientific importance, feasibility, or probability score.",
        "",
        "## Exit contracts",
        "",
    ]
    for path, contract in queue["path_contracts"].items():
        if "owner_issue" in contract:
            owner = f"[#{contract['owner_issue']}]({ISSUE_URL}/{contract['owner_issue']})"
        else:
            owner = (
                f"`{contract['owner_registry']}` with audits in "
                f"`{contract['audit_registry']}` (bootstrap "
                f"[#{contract['bootstrap_issue']}]({ISSUE_URL}/{contract['bootstrap_issue']}))"
            )
        lines += [
            f"### {path}",
            "",
            f"Standing owner: {owner}.",
            "",
            "Closure requires:",
            "",
        ]
        lines.extend(f"- {requirement}" for requirement in contract["closure_requirements"])
        lines += ["", "Does not close the item:", ""]
        lines.extend(f"- {boundary}" for boundary in contract["non_closing_outputs"])
        lines.append("")

    lines += ["## Per-premise custody", ""]
    for item in queue["items"]:
        lines += [
            f"### {item['queue_id']} — {item['premise_name']}",
            "",
            f"- State: `{item['state']}`; origin `{item['origin_architecture_version']}`; current evaluation `{item['evaluation_architecture_version']}`.",
            f"- Origin premise-record SHA-256: `{item['origin_premise_record_sha256']}`.",
            f"- Disposition: `{item['disposition']}`; allowed exits: {_paths(item['available_exit_paths'])}.",
            f"- Exact reverse consumers: {_edge_summary(item['consumer_edges'])}.",
            f"- Standing decision custody: `{item['decision_custody']['architecture_register']}` with independent records in `{item['decision_custody']['audit_register']}`; [#{item['decision_custody']['bootstrap_issue']}]({ISSUE_URL}/{item['decision_custody']['bootstrap_issue']}) is historical bootstrap custody only.",
            f"- Registered statement: {item['premise_statement']}",
            f"- Next action or decision: {item['next_action_or_decision']}",
            "- Immutable origin evidence:",
            "",
        ]
        for pin in item["origin_evidence"]:
            lines.append(
                f"  - `{pin['path']}` — `{pin['git_object_type']}` "
                f"`{pin['git_object_sha1']}`, SHA-256 `{pin['content_sha256']}`, "
                f"{pin['byte_count']} pin-payload bytes at `{pin['origin_revision']}`."
            )
        references = item["consumer_edges"]["evidence_references"]
        if references:
            lines += ["", "  Shared-evidence references (not inferred premise edges):"]
            for reference in references:
                observations = ", ".join(reference["observation_row_ids"]) or "none"
                claims = ", ".join(reference["claim_ids"]) or "none"
                lines.append(
                    f"  - `{reference['path']}` — observation rows {observations}; claims {claims}."
                )
        unreferenced = item["consumer_edges"]["unreferenced_evidence_paths"]
        if unreferenced:
            lines += ["", "  Evidence paths with no canonical reverse reference (no consumer inferred):"]
            lines.extend(f"  - `{path}`" for path in unreferenced)
        lines.append("")

    lines += [
        "## Flagged imports outside #739",
        "",
        "Imports remain explicit inputs and are not silently transformed into derivation obligations. A recorded disposition change would place the affected row into a versioned successor queue.",
        "",
        "| Premise | Type | Reason |",
        "| --- | --- | --- |",
    ]
    for row in queue["excluded_imports"]:
        lines.append(
            f"| {row['premise_id']} {row['premise_name']} | `{row['premise_type']}` | {row['reason']} |"
        )

    lines += [
        "",
        "## Inherited V2 custody",
        "",
        "The live #739 contract preserves these superseded issue records. They remain historical evidence, not automatically one-to-one aliases for the newer atomic premise rows:",
        "",
    ]
    for row in queue["inherited_v2_residuals"]:
        lines.append(
            f"- {row['code']} [#{row['issue']}]({ISSUE_URL}/{row['issue']}): {row['residual']}."
        )
    lines += [
        "",
        "A bounded negative may be retained as evidence, but the item stays open unless a positive derivation lands or an eligible axiomatization is recorded. Deferral is never discharge.",
    ]
    return "\n".join(lines) + "\n"


def _json_bytes(queue: dict) -> bytes:
    return (json.dumps(queue, indent=2, ensure_ascii=False) + "\n").encode("utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail when either committed generated surface has drifted",
    )
    args = parser.parse_args()

    sources = load_sources()
    expected = build_queue(sources)
    queue_bytes = _json_bytes(expected)
    surface_bytes = render(expected).encode("utf-8")

    if args.check:
        current_queue = QUEUE_PATH.read_bytes() if QUEUE_PATH.is_file() else b""
        current_surface = SURFACE_PATH.read_bytes() if SURFACE_PATH.is_file() else b""
        if current_queue != queue_bytes or current_surface != surface_bytes:
            print(
                "premise discharge queue: generated surfaces are stale; run "
                "python tools/build_premise_discharge_queue.py",
                file=sys.stderr,
            )
            return 1
        validate_queue(load_json(QUEUE_PATH), sources)
        print(
            "premise discharge queue: surfaces are current "
            f"({len(expected['items'])} actionable items)"
        )
        return 0

    QUEUE_PATH.write_bytes(queue_bytes)
    SURFACE_PATH.write_bytes(surface_bytes)
    print(
        "premise discharge queue: wrote "
        f"{QUEUE_PATH.relative_to(ROOT)} and {SURFACE_PATH.relative_to(ROOT)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
