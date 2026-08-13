"""Build and validate the V3 observation ledger surface (issue #726).

The machine-readable ledger is ``tracking/observation_ledger.json``: one row
per physical observation the architecture must reproduce, with an adequacy
rung, a conservative status, the owning composition lane, the premise-register
rows the current status consumes, and evidence links. This tool validates the
ledger fail-closed (exact key set, rung and status enums, lane range, premise
ids drawn from the fixed register, evidence paths resolving to committed
files) and renders ``docs/OBSERVATION_LEDGER_V3.md`` from it. The rendered
page is a generated surface: ``--check`` fails when the committed page differs
byte for byte from the render.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path

import build_audit_custody
import build_architecture_versions
import prediction_lineage_custody
import strict_json

ROOT = Path(__file__).resolve().parents[1]
LEDGER_PATH = ROOT / "tracking" / "observation_ledger.json"
SURFACE_PATH = ROOT / "docs" / "OBSERVATION_LEDGER_V3.md"
PREMISE_REGISTER_PATH = ROOT / "tracking" / "premise_register.json"
ARCHITECTURE_REGISTER_PATH = ROOT / "tracking" / "architecture_versions.json"
AUDIT_CUSTODY_PATH = ROOT / "tracking" / "audit_custody.json"
PREDICTION_LINEAGE_PATH = (
    ROOT / "claims" / "frozen_prediction_architecture_lineages.json"
)
FROZEN_PREDICTION_PATH = ROOT / "claims" / "frozen_prediction_register.json"

SCHEMA = "oph.observation_ledger.v3"
ISSUE = 726
REPO_URL = "https://github.com/FloatingPragma/observer-patch-holography"

RUNGS = ("formal_precursor", "structural", "emergent", "predictive")
STATUSES = ("attained", "partial", "owed")
LANE_ISSUES = frozenset((*range(728, 739), 740, 742, 743, 744, 745))

ID_PATTERN = re.compile(r"^OL-[A-N][1-9]$")
PREMISE_PATTERN = re.compile(r"^PR-\d{2}$")
BANNED_CHARACTERS = ("—", "–")

EXPECTED_ROW_IDS = (
    "OL-A1",
    "OL-A2",
    "OL-A3",
    "OL-A4",
    "OL-B1",
    "OL-B2",
    "OL-B3",
    "OL-B4",
    "OL-B5",
    "OL-C1",
    "OL-C2",
    "OL-C3",
    "OL-C4",
    "OL-C5",
    "OL-C6",
    "OL-D1",
    "OL-D2",
    "OL-D3",
    "OL-E1",
    "OL-E2",
    "OL-E3",
    "OL-E4",
    "OL-F1",
    "OL-F2",
    "OL-F3",
    "OL-F4",
    "OL-G1",
    "OL-G2",
    "OL-G3",
    "OL-G4",
    "OL-G5",
    "OL-G8",
    "OL-G6",
    "OL-G7",
    "OL-G9",
    "OL-H1",
    "OL-H2",
    "OL-H3",
    "OL-H4",
    "OL-H5",
    "OL-H6",
    "OL-H7",
    "OL-H8",
    "OL-I1",
    "OL-I2",
    "OL-I3",
    "OL-J1",
    "OL-J2",
    "OL-J3",
    "OL-K1",
    "OL-K2",
    "OL-K3",
    "OL-L1",
    "OL-L2",
    "OL-L3",
    "OL-M1",
    "OL-M2",
    "OL-M3",
    "OL-N1",
)

# These premise packets were the subject of a concrete independent
# hidden-premise audit.  Keep their row-level ancestry explicit: lane-level
# reverse-map coverage alone would not detect moving a premise to another row
# in the same lane.
AUDITED_ROW_PREMISE_CONTRACTS = {
    "OL-E1": {
        "premises": ("PR-07", "PR-15"),
        "open_premises": ("PR-08",),
    },
}

# Premise-register lane ownership is broader than this observation table: a
# lane can consume a premise in a correspondence table, theorem surface,
# constants/frozen register, or architecture protocol without attaching it to
# a particular observation row.  Every such pair is explicit here so a newly
# uncovered reverse-map mismatch fails closed.
NON_OBSERVATION_SURFACE_CONSUMERS = {
    (728, "PR-01"): "architecture and spacetime composition surface",
    (729, "PR-01"): "Einstein composition surface",
    (730, "PR-01"): "quantum composition surface",
    (732, "PR-01"): "thermodynamic composition surface",
    (735, "PR-13"): "SM correspondence and mass-sector boundary",
    (729, "PR-14"): "gravity/constants comparison boundary",
    (735, "PR-15"): "SM correspondence calibration boundary",
    (729, "PR-16"): "gravity ladder attachment boundary",
    (729, "PR-29"): "layered gravity source surface",
    (729, "PR-30"): "layered gravity source surface",
    (729, "PR-31"): "layered gravity source surface",
    (729, "PR-38"): "P-closure gravity/constants crosswalk",
    (729, "PR-39"): "P-closure gravity/constants crosswalk",
    (736, "PR-53"): "constants physical-comparison boundary",
    (738, "PR-53"): "frozen-instrument physical-comparison custody",
    (744, "PR-55"): "baryon-label/proton surface",
    (735, "PR-56"): "SM correspondence operator census",
    (744, "PR-56"): "baryon-operator/proton surface",
    (735, "PR-57"): "SM correspondence proton boundary",
    (744, "PR-57"): "physical proton/effective-action surface",
}

ROW_KEYS = {
    "id",
    "target",
    "rung",
    "status",
    "lane_issue",
    "architecture_version",
    "premises",
    "open_premises",
    "evidence",
    "notes",
}
PREDICTIVE_ROW_KEYS = ROW_KEYS | {"prediction_event"}

GROUPS = (
    ("Spacetime", (728,)),
    ("Gravitation", (729,)),
    ("Quantum", (730,)),
    ("Mechanics", (731,)),
    ("Thermodynamics", (732,)),
    ("Electromagnetism", (733,)),
    ("Standard Model structure", (734,)),
    ("Standard Model Lagrangian", (735,)),
    ("Masses and constants", (736,)),
    ("Cosmology and instruments", (737, 738)),
    ("Common-world integration", (740,)),
    ("Cosmology and astrophysics", (742,)),
    ("Interacting quantum field theory", (743,)),
    ("QCD, hadrons, and nuclei", (744,)),
    ("Electroweak and weak phenomenology", (745,)),
)


def fail(message: str) -> None:
    raise SystemExit(f"observation ledger: {message}")


def load_ledger(path: Path) -> dict:
    try:
        return strict_json.load(path)
    except FileNotFoundError:
        fail(f"missing input {path}")
    except (json.JSONDecodeError, strict_json.DuplicateKeyError) as error:
        fail(f"invalid JSON in {path}: {error}")
    raise AssertionError("unreachable")


def load_premise_register() -> tuple[list[dict], dict[str, dict]]:
    register = load_ledger(PREMISE_REGISTER_PATH)
    rows = register.get("rows")
    if not isinstance(rows, list) or not rows:
        fail("premise register rows must be a nonempty list")
    by_id: dict[str, dict] = {}
    for row in rows:
        if not isinstance(row, dict) or not isinstance(row.get("id"), str):
            fail("premise register contains a malformed row")
        if row["id"] in by_id:
            fail(f"premise register repeats {row['id']}")
        by_id[row["id"]] = row
    return rows, by_id


def _clean_prose(where: str, field: str, value: object) -> str:
    if not isinstance(value, str) or not value.strip():
        fail(f"{where}: {field} must be a nonempty string")
    for character in BANNED_CHARACTERS:
        if character in value:
            fail(f"{where}: {field} carries a banned dash character")
    return value


def _current_evidence_sha256(path: str) -> str:
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def validate(ledger: dict) -> list[dict]:
    if not isinstance(ledger, dict):
        fail("ledger must be an object")
    if set(ledger) != {"schema", "issue", "audit_pointers", "rows"}:
        fail("top-level keys must be exactly schema, issue, audit_pointers, rows")
    if ledger["schema"] != SCHEMA:
        fail(f"schema must equal {SCHEMA}")
    if ledger["issue"] != ISSUE:
        fail(f"issue must equal {ISSUE}")
    rows = ledger["rows"]
    if not isinstance(rows, list) or not rows:
        fail("rows must be a nonempty list")
    row_ids = [row.get("id") if isinstance(row, dict) else None for row in rows]
    if tuple(row_ids) != EXPECTED_ROW_IDS:
        fail("row ids must equal the fixed ordered observation inventory")

    audit_data = build_audit_custody.load_json(AUDIT_CUSTODY_PATH)
    audit_records = build_audit_custody.validate(audit_data)
    audit_by_id = {record["id"]: record for record in audit_records}
    audit_row_indexes = {
        record["id"]: build_audit_custody.reviewed_row_index(record)
        for record in audit_records
    }
    audit_pointers = ledger["audit_pointers"]
    if not isinstance(audit_pointers, dict):
        fail("audit_pointers must be an object")

    _, premise_by_id = load_premise_register()
    architecture_register = build_architecture_versions.load_json(
        ARCHITECTURE_REGISTER_PATH
    )
    build_architecture_versions.validate(architecture_register)
    architecture_ids = {
        version.get("id")
        for version in architecture_register.get("versions", [])
        if isinstance(version, dict)
    }
    if not architecture_ids:
        fail("architecture version register must contain at least one version")
    anchored_architecture_ids = {
        anchor.get("id")
        for anchor in architecture_register.get("version_anchors", [])
        if isinstance(anchor, dict)
    }
    seen_ids: set[str] = set()
    for index, row in enumerate(rows):
        where = f"rows[{index}]"
        if not isinstance(row, dict):
            fail(f"{where}: row must be an object")
        expected_keys = PREDICTIVE_ROW_KEYS if row.get("rung") == "predictive" else ROW_KEYS
        if set(row) != expected_keys:
            missing = expected_keys - set(row)
            extra = set(row) - ROW_KEYS
            if row.get("rung") == "predictive":
                extra = set(row) - PREDICTIVE_ROW_KEYS
            fail(
                f"{where}: keys mismatch "
                f"(missing {sorted(missing)}, extra {sorted(extra)})"
            )
        row_id = row["id"]
        if not isinstance(row_id, str) or not ID_PATTERN.match(row_id):
            fail(f"{where}: id must match OL-<A..N><digit>")
        if row_id in seen_ids:
            fail(f"duplicate row id {row_id}")
        seen_ids.add(row_id)
        where = row_id

        _clean_prose(where, "target", row["target"])
        _clean_prose(where, "notes", row["notes"])
        if row["rung"] not in RUNGS:
            fail(f"{where}: rung must be one of {RUNGS}")
        if row["status"] not in STATUSES:
            fail(f"{where}: status must be one of {STATUSES}")
        if row["rung"] == "predictive":
            pointer = row["prediction_event"]
            if pointer is not None and (
                not isinstance(pointer, dict)
                or set(pointer)
                != prediction_lineage_custody.PREDICTION_EVENT_POINTER_KEYS
            ):
                fail(f"{where}: prediction_event must be null or an exact event pointer")

        lane = row["lane_issue"]
        if not isinstance(lane, int) or isinstance(lane, bool):
            fail(f"{where}: lane_issue must be an integer")
        if lane not in LANE_ISSUES:
            fail(f"{where}: lane_issue {lane} is not a registered V3 lane")

        architecture_version = row["architecture_version"]
        if architecture_version not in architecture_ids:
            fail(
                f"{where}: architecture_version {architecture_version!r} is not "
                "on the architecture version register"
            )
        if row["status"] == "attained" and architecture_version not in anchored_architecture_ids:
            fail(f"{where}: an attained row cannot use an unanchored architecture tip")

        for field in ("premises", "open_premises"):
            premises = row[field]
            if not isinstance(premises, list):
                fail(f"{where}: {field} must be a list")
            if len(premises) != len(set(premises)):
                fail(f"{where}: {field} must be duplicate-free")
            for premise in premises:
                if not isinstance(premise, str) or not PREMISE_PATTERN.match(premise):
                    fail(f"{where}: premise ids must match PR-<two digits>")
                if premise not in premise_by_id:
                    fail(f"{where}: premise {premise} is not on the canonical register")
                consuming = premise_by_id[premise].get("consuming_lanes")
                if not isinstance(consuming, list) or lane not in consuming:
                    fail(
                        f"{where}: premise {premise} does not declare consuming "
                        f"lane #{lane}"
                    )
        if set(row["premises"]) & set(row["open_premises"]):
            fail(f"{where}: consumed and open premise lists must be disjoint")
        if row["status"] == "owed" and row["premises"]:
            fail(f"{where}: an owed row carries no consumed premises")
        contract = AUDITED_ROW_PREMISE_CONTRACTS.get(row["id"])
        if contract is not None:
            for field, expected in contract.items():
                if tuple(row[field]) != expected:
                    fail(
                        f"{where}: {field} must equal the independently audited "
                        f"row-level contract {list(expected)}"
                    )

        evidence = row["evidence"]
        if not isinstance(evidence, list):
            fail(f"{where}: evidence must be a list")
        if len(evidence) != len(set(evidence)):
            fail(f"{where}: evidence must be duplicate-free")
        for path in evidence:
            if not isinstance(path, str) or not path:
                fail(f"{where}: evidence entries must be nonempty strings")
            if path.startswith("/") or ".." in path.split("/"):
                fail(f"{where}: evidence path {path} must be repo-relative")
            if not (ROOT / path).is_file():
                fail(f"{where}: evidence path missing: {path}")

    attained_ids = {row["id"] for row in rows if row["status"] == "attained"}
    if set(audit_pointers) != attained_ids:
        fail("audit_pointers keys must exactly equal the attained observation rows")
    row_by_id = {row["id"]: row for row in rows}
    for row_id, audit_ids in audit_pointers.items():
        if not isinstance(audit_ids, list) or not audit_ids:
            fail(f"{row_id}: audit pointer list must be nonempty")
        if len(audit_ids) != len(set(audit_ids)):
            fail(f"{row_id}: audit pointer list must be duplicate-free")
        current = row_by_id[row_id]
        qualified = False
        for audit_id in audit_ids:
            if not isinstance(audit_id, str) or audit_id not in audit_by_id:
                fail(f"{row_id}: unknown audit record {audit_id!r}")
            record = audit_by_id[audit_id]
            if row_id not in record["reviewed_rows"]:
                fail(f"{row_id}: {audit_id} did not review this row")
            historical = audit_row_indexes[audit_id].get(row_id)
            if historical is None:
                fail(f"{row_id}: {audit_id} has no historical row payload")
            current_projection = dict(current)
            historical_projection = dict(historical)
            historical_projection.pop("audit_records", None)
            if (
                row_id in record["promoted_rows"]
                and historical_projection == current_projection
            ):
                pins = {
                    (pin["revision"], pin["path"]): pin
                    for pin in record["artifact_pins"]
                }
                drifted = [
                    path
                    for path in current["evidence"]
                    if _current_evidence_sha256(path)
                    != pins[(record["repair_commit"], path)]["sha256"]
                ]
                if drifted:
                    fail(
                        f"{row_id}: current evidence bytes drifted from "
                        f"{audit_id} at its repair commit: {drifted}"
                    )
                qualified = True
        if not qualified:
            fail(
                f"{row_id}: no audit pointer qualifies its exact historical row payload"
            )

    grouped_lanes = [lane for _, lanes in GROUPS for lane in lanes]
    if set(grouped_lanes) != set(LANE_ISSUES) or len(grouped_lanes) != len(LANE_ISSUES):
        fail("group table must cover every lane exactly once")
    ledger_pairs = {
        (row["lane_issue"], premise)
        for row in rows
        for premise in (*row["premises"], *row["open_premises"])
    }
    declared_pairs = {
        (lane, premise_id)
        for premise_id, premise in premise_by_id.items()
        for lane in premise.get("consuming_lanes", [])
    }
    uncovered = declared_pairs - ledger_pairs
    exceptions = set(NON_OBSERVATION_SURFACE_CONSUMERS)
    if uncovered != exceptions:
        fail(
            "premise reverse-map exceptions drifted: missing exceptions "
            f"{sorted(uncovered - exceptions)}, stale exceptions "
            f"{sorted(exceptions - uncovered)}"
        )

    lineage_data = prediction_lineage_custody.load_json(PREDICTION_LINEAGE_PATH)
    frozen_data = prediction_lineage_custody.load_json(FROZEN_PREDICTION_PATH)
    lineage_state = prediction_lineage_custody.validate(
        lineage_data,
        frozen_data,
        architecture_register,
        rows,
        root=ROOT,
    )
    prediction_lineage_custody.require_predictive_promotions(rows, lineage_state)
    return rows


def _lane_link(lane: int) -> str:
    return f"[#{lane}]({REPO_URL}/issues/{lane})"


def _issue_link(number: int) -> str:
    return f"[issue #{number}]({REPO_URL}/issues/{number})"


def render(rows: list[dict], audit_pointers: dict[str, list[str]] | None = None) -> str:
    if audit_pointers is None:
        canonical = load_ledger(LEDGER_PATH)
        audit_pointers = canonical.get("audit_pointers", {})
    lines: list[str] = []
    lines.append("# V3 Observation Ledger")
    lines.append("")
    lines.append(
        "Generated by `tools/build_observation_ledger.py` from"
        " `tracking/observation_ledger.json`; edit the JSON, then regenerate."
    )
    lines.append("")
    lines.append(
        f"One row per observation the architecture must reproduce, tracked for"
        f" {_issue_link(726)}. Premise ids PR-01 through"
        f" PR-{len(load_premise_register()[0]):02d} name rows of the"
        f" premise register ({_issue_link(727)}), the anti-cheating"
        f" surface; each row lists the register rows its current status"
        f" consumes. The row set, rungs, and lane assignments follow completion"
        f" plan V3 (`plan/COMPLETION_PLAN_V3.md` in the oph-meta planning"
        f" workspace). The ledger records adequacy status. It is not itself"
        f" evidence: a promotion requires the owning lane's receipts."
    )
    lines.append("")
    lines.append("## Adequacy rungs")
    lines.append("")
    lines.append(
        "- **Formal precursor**: an exact mathematical representation or"
        " helper exists, but no theorem attaches it to the observer"
        " architecture or a physical readout. It cannot by itself satisfy a"
        " physics target."
    )
    lines.append(
        "- **Structural**: the architecture provably carries the law; an exact"
        " theorem from the three axioms plus named register rows,"
        " machine-checked where the statement permits."
    )
    lines.append(
        "- **Emergent**: simulated observers' records exhibit the law; a"
        " preregistered, target-clean instrument run at laptop scale with"
        " pinned receipts."
    )
    lines.append(
        "- **Predictive**: a frozen instrument binds the row to future data"
        " with registered kill bands, under the custody rules of the standing"
        " lane."
    )
    lines.append("")
    lines.append(
        "A structural claim whose premises are unregistered is invalid, an"
        " emergent claim without a preregistered instrument is invalid, and a"
        " predictive claim without custody is invalid."
    )
    lines.append("")
    lines.append("## Status labels")
    lines.append("")
    lines.append(
        "- **attained**: the rung's contract is met inside the stated"
        " boundary, with committed receipts."
    )
    lines.append(
        "- **partial**: committed receipts cover part of the row's contract;"
        " the boundary column states what is owed."
    )
    lines.append(
        "- **owed**: no qualifying receipt is committed; the owning lane"
        " carries the obligation."
    )
    lines.append("")
    lines.append(
        "Statuses are conservative: a conditional result is conditional, a"
        " declared premise is declared, and postdictions are never"
        " predictions."
    )

    lineage_data = prediction_lineage_custody.load_json(PREDICTION_LINEAGE_PATH)
    lines.extend(
        [
            "",
            "## Predictive lineage custody",
            "",
            "Every predictive row has an explicit mapping in "
            "`claims/frozen_prediction_architecture_lineages.json`. Historical "
            "or pending baseline rows do not qualify a promotion. An attained "
            "predictive row requires an append-only freeze event anchored to "
            "its first-appearance commit and bound to the current anchored AV-n.",
            "",
            "| Row | Pending candidates | Historical only | No-candidate boundary |",
            "| --- | --- | --- | --- |",
        ]
    )
    for binding in lineage_data["predictive_observation_bindings"]:
        candidates = ", ".join(binding["candidate_baseline_lineages"]) or "none"
        historical = ", ".join(binding["historical_only_lineages"]) or "none"
        reason = binding["no_candidate_reason"] or "none"
        lines.append(
            f"| {binding['observation_row_id']} | {candidates} | {historical} | {reason} |"
        )

    for title, lanes in GROUPS:
        group_rows = [row for row in rows if row["lane_issue"] in lanes]
        if not group_rows:
            continue
        lane_word = "lane" if len(lanes) == 1 else "lanes"
        lane_list = ", ".join(_lane_link(lane) for lane in lanes)
        lines.append("")
        lines.append(f"## {title} ({lane_word} {lane_list})")
        lines.append("")
        lines.append(
            "| Row | Observation | Rung | Status | Architecture | Audit | Lane | Premises |"
            " Open premises | Evidence | Boundary |"
        )
        lines.append("| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |")
        for row in group_rows:
            premises = ", ".join(row["premises"]) if row["premises"] else "none"
            open_premises = (
                ", ".join(row["open_premises"]) if row["open_premises"] else "none"
            )
            evidence = (
                ", ".join(f"`{path}`" for path in row["evidence"])
                if row["evidence"]
                else "none"
            )
            audits = ", ".join(audit_pointers.get(row["id"], [])) or "none"
            lines.append(
                f"| {row['id']} | {row['target']} | {row['rung']} |"
                f" {row['status']} | {row['architecture_version']} |"
                f" {audits} |"
                f" {_lane_link(row['lane_issue'])} | {premises} |"
                f" {open_premises} | {evidence} |"
                f" {row['notes']} |"
            )

    status_counts = {status: 0 for status in STATUSES}
    rung_counts = {rung: 0 for rung in RUNGS}
    for row in rows:
        status_counts[row["status"]] += 1
        rung_counts[row["rung"]] += 1
    lines.append("")
    lines.append(
        f"Totals: {len(rows)} rows. Status: "
        + ", ".join(f"{status_counts[status]} {status}" for status in STATUSES)
        + ". Rung: "
        + ", ".join(f"{rung_counts[rung]} {rung}" for rung in RUNGS)
        + "."
    )
    lines.append("")
    lines.append("## Premises consumed or still open")
    lines.append("")
    consumed = sorted(
        {
            premise
            for row in rows
            for premise in (*row["premises"], *row["open_premises"])
        }
    )
    premise_by_id = load_premise_register()[1]
    for premise in consumed:
        entry = premise_by_id[premise]
        name = entry["name"]
        kind = entry["type"]
        disposition = entry["disposition"]
        lines.append(f"- **{premise}** {name} ({kind}, {disposition})")
    lines.append("")
    lines.append(
        f"Register rows that no current or open target consumes are omitted"
        f" from this list. `Premises` records hypotheses used by the attained"
        f" portion; `Open premises` records registered hypotheses still needed"
        f" for the row's full target. The full register lives"
        f" under {_issue_link(727)}."
    )
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail when the committed surface differs from the render",
    )
    args = parser.parse_args(argv)

    ledger = load_ledger(LEDGER_PATH)
    rows = validate(ledger)
    surface = render(rows, ledger["audit_pointers"]).encode("utf-8")
    if args.check:
        committed = SURFACE_PATH.read_bytes() if SURFACE_PATH.is_file() else b""
        if committed != surface:
            print(
                "observation ledger: docs/OBSERVATION_LEDGER_V3.md is stale;"
                " run python tools/build_observation_ledger.py",
                file=sys.stderr,
            )
            return 1
        print("observation ledger: surface is current")
        return 0
    SURFACE_PATH.write_bytes(surface)
    print(f"observation ledger: wrote {SURFACE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
