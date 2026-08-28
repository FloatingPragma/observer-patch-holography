"""Build and validate the V3 observation ledger surface (bootstrap issue #726).

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
import json
import re
import sys
from pathlib import Path

import strict_json

ROOT = Path(__file__).resolve().parents[1]
LEDGER_PATH = ROOT / "tracking" / "observation_ledger.json"
SURFACE_PATH = ROOT / "docs" / "OBSERVATION_LEDGER_V3.md"
PREMISE_REGISTER_PATH = ROOT / "tracking" / "premise_register.json"
FROZEN_REGISTER_PATH = ROOT / "claims" / "frozen_prediction_register.json"

SCHEMA = "oph.observation_ledger.v4"
ISSUE = 726
REPO_URL = "https://github.com/FloatingPragma/observer-patch-holography"

RUNGS = ("formal_precursor", "structural", "emergent", "predictive")
STATUSES = ("attained", "partial", "owed")
LANE_ISSUES = frozenset((*range(728, 738), 740, 742, 743, 744, 745))

ID_PATTERN = re.compile(r"^OL-[A-N][1-9]$")
PREMISE_PATTERN = re.compile(r"^PR-\d{2}$")
FROZEN_TARGET_PATTERN = re.compile(r"^FZ-\d{2}$")
BANNED_CHARACTERS = ("—", "–")
FROZEN_TARGET_STATUSES = {
    "frozen_attested",
    "frozen_stamped_upgrade_pending",
    "standing_frozen",
}

# Exact scientific relationship between predictive observation rows and the
# frozen targets that actually address them.  Existence alone is insufficient:
# an unrelated frozen target must never qualify a row.  `all` means every
# listed target is part of the row contract; `any` is reserved for a future
# row whose alternatives are explicitly interchangeable.
PREDICTIVE_TARGET_CONTRACTS = {
    "OL-B4": {"targets": ("FZ-06", "FZ-14"), "attainment": "all"},
    "OL-F4": {"targets": ("FZ-11", "FZ-12"), "attainment": "all"},
    "OL-H7": {"targets": ("FZ-03",), "attainment": "all"},
    "OL-I2": {"targets": ("FZ-09",), "attainment": "all"},
    "OL-I3": {"targets": ("FZ-07", "FZ-08"), "attainment": "all"},
}

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
    "OL-I2",
    "OL-I3",
    "OL-J1",
    "OL-J2",
    "OL-J3",
    "OL-K1",
    "OL-K2",
    "OL-K3",
    "OL-K4",
    "OL-K5",
    "OL-K6",
    "OL-L1",
    "OL-L2",
    "OL-L3",
    "OL-M1",
    "OL-M2",
    "OL-M3",
    "OL-N1",
)

# Keep this load-bearing row-level ancestry explicit: lane-level reverse-map
# coverage alone would not detect moving a premise to another row in the same
# lane.
ROW_PREMISE_CONTRACTS = {
    "OL-D1": {
        # The composed theorem consumes these antecedents, but the row stays
        # partial: the committed witness is one declared one-parameter
        # enrichment, not a source-selected or generally unique least action.
        "premises": ("PR-05", "PR-06", "PR-45"),
        "open_premises": (),
    },
    "OL-E1": {
        # 2026-08-14 composed promotion: the refinement-uniform third law is
        # threaded into fourLaws_composed, so PR-08 moved from open to
        # consumed on this row.
        "premises": ("PR-07", "PR-15", "PR-08"),
        "open_premises": (),
    },
    "OL-G2": {
        "premises": ("PR-35", "PR-59"),
        "open_premises": ("PR-46",),
    },
    "OL-G3": {
        "premises": ("PR-59",),
        "open_premises": ("PR-47", "PR-54"),
    },
    "OL-G4": {
        "premises": ("PR-36", "PR-59"),
        "open_premises": ("PR-47", "PR-54"),
    },
    "OL-G5": {
        "premises": ("PR-59",),
        "open_premises": ("PR-47",),
    },
    "OL-G8": {
        "premises": ("PR-55", "PR-56", "PR-59"),
        "open_premises": ("PR-47", "PR-54", "PR-57"),
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
    (734, "PR-12"): "separately frozen B15 matter-search negative control, not a PR-59 antecedent",
    (735, "PR-12"): "separately frozen B15 matter-search negative control, not a PR-59 antecedent",
    (729, "PR-14"): "gravity/constants comparison boundary",
    (742, "PR-70"): "frozen dispersion comparison contract surface of the prediction ladder, not a cosmology observation row",
    (742, "PR-71"): "frozen dispersion comparison contract surface of the prediction ladder, not a cosmology observation row",
    (742, "PR-72"): "frozen dispersion comparison contract surface of the prediction ladder, not a cosmology observation row",
    (742, "PR-73"): "frozen dispersion comparison contract surface of the prediction ladder, not a cosmology observation row",
    (742, "PR-74"): "frozen dispersion comparison contract surface of the prediction ladder, not a cosmology observation row",
    (742, "PR-75"): "frozen dispersion comparison contract surface of the prediction ladder, not a cosmology observation row",
    (742, "PR-76"): "frozen dispersion comparison contract surface of the prediction ladder, not a cosmology observation row",
    (750, "PR-77"): "adaptive consensus theorem premise, not a physical observation row",
    (750, "PR-78"): "adaptive consensus rank-bound premise, not a physical observation row",
    (735, "PR-15"): "SM correspondence calibration boundary",
    (729, "PR-16"): "gravity ladder attachment boundary",
    (729, "PR-29"): "layered gravity source surface",
    (729, "PR-30"): "layered gravity source surface",
    (729, "PR-31"): "layered gravity source surface",
    (729, "PR-38"): "P-closure gravity/constants crosswalk",
    (729, "PR-39"): "P-closure gravity/constants crosswalk",
    (736, "PR-53"): "constants physical-comparison boundary",
    (744, "PR-55"): "baryon-label/proton surface",
    (735, "PR-56"): "SM correspondence operator census",
    (744, "PR-56"): "baryon-operator/proton surface",
    (735, "PR-57"): "SM correspondence proton boundary",
    (744, "PR-57"): "physical proton/effective-action surface",
    (736, "PR-59"): "mass and constants exterior-table boundary",
    (744, "PR-59"): "QCD and baryon-census exterior-table boundary",
    (745, "PR-59"): "electroweak exterior-table boundary",
    (740, "PR-18"): "common-world island-bridge shell-readout surface; OL-N1 stays owed and carries no consumed premises",
    (740, "PR-19"): "common-world island-bridge shell-readout surface; OL-N1 stays owed and carries no consumed premises",
    (740, "PR-67"): "common-world finite electroweak grammar surface; OL-N1 stays owed and carries no consumed premises",
    (740, "PR-68"): "common-world finite breaking-packet surface; OL-N1 stays owed and carries no consumed premises",
    (740, "PR-69"): "common-world finite Yukawa-line surface; OL-N1 stays owed and carries no consumed premises",
    (745, "PR-67"): "tree-level neutral and charged current dictionaries; physical observation rows stay owed",
    (745, "PR-68"): "tree-level neutral and charged current dictionaries; physical observation rows stay owed",
    (750, "PR-60"): "declared scheduler-class surface; the consuming lane has no observation row",
    (750, "PR-61"): "declared scheduler-class surface; the consuming lane has no observation row",
    (750, "PR-62"): "declared scheduler-class surface; the consuming lane has no observation row",
    (750, "PR-63"): "declared scheduler-class capacity-split surface; the consuming lane has no observation row",
}

BASE_ROW_KEYS = {
    "id",
    "target",
    "rung",
    "status",
    "lane_issue",
    "premises",
    "open_premises",
    "evidence",
    "notes",
}
PREDICTIVE_ROW_KEYS = BASE_ROW_KEYS | {"frozen_targets"}

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
    ("Cosmology and instruments", (737,)),
    ("Common-world integration", (740,)),
    ("Cosmology and astrophysics", (742,)),
    ("Interacting quantum field theory", (743,)),
    ("QCD, hadrons, and nuclei", (744,)),
    ("Electroweak and weak phenomenology", (745,)),
)



def _cell(text: object) -> str:
    """Escape pipes so free text cannot break a Markdown table row."""
    return str(text).replace("|", "\\|")

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


def load_frozen_target_rows() -> dict[str, dict]:
    register = load_ledger(FROZEN_REGISTER_PATH)
    rows = register.get("rows") if isinstance(register, dict) else None
    if not isinstance(rows, list) or not rows:
        fail("frozen-prediction register rows must be a nonempty list")
    by_id: dict[str, dict] = {}
    for row in rows:
        if not isinstance(row, dict):
            fail("frozen-prediction register contains a malformed row")
        row_id = row.get("id")
        status = row.get("status")
        if not isinstance(row_id, str) or not FROZEN_TARGET_PATTERN.fullmatch(row_id):
            fail("frozen-prediction register contains a malformed target id")
        if not isinstance(status, str) or not status:
            fail(f"frozen-prediction register row {row_id} has no status")
        if row_id in by_id:
            fail(f"frozen-prediction register repeats {row_id}")
        by_id[row_id] = row
    return by_id


def _clean_prose(where: str, field: str, value: object) -> str:
    if not isinstance(value, str) or not value.strip():
        fail(f"{where}: {field} must be a nonempty string")
    for character in BANNED_CHARACTERS:
        if character in value:
            fail(f"{where}: {field} carries a banned dash character")
    return value


def validate(ledger: dict) -> list[dict]:
    if not isinstance(ledger, dict):
        fail("ledger must be an object")
    if set(ledger) != {"schema", "issue", "rows"}:
        fail("top-level keys must be exactly schema, issue, rows")
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

    _, premise_by_id = load_premise_register()
    frozen_target_by_id = load_frozen_target_rows()
    seen_ids: set[str] = set()
    for index, row in enumerate(rows):
        where = f"rows[{index}]"
        if not isinstance(row, dict):
            fail(f"{where}: row must be an object")
        expected_keys = (
            PREDICTIVE_ROW_KEYS
            if row.get("rung") == "predictive"
            else BASE_ROW_KEYS
        )
        if set(row) != expected_keys:
            missing = expected_keys - set(row)
            extra = set(row) - expected_keys
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
        if row["status"] == "attained" and row["open_premises"]:
            fail(f"{where}: an attained row cannot retain open premises")
        if row["rung"] == "predictive":
            targets = row["frozen_targets"]
            if not isinstance(targets, list):
                fail(f"{where}: frozen_targets must be a list")
            for target in targets:
                if (
                    not isinstance(target, str)
                    or not FROZEN_TARGET_PATTERN.fullmatch(target)
                ):
                    fail(f"{where}: frozen target ids must match FZ-<two digits>")
                if target not in frozen_target_by_id:
                    fail(f"{where}: frozen target {target} is not on the register")
            if len(targets) != len(set(targets)):
                fail(f"{where}: frozen_targets must be duplicate-free")
            contract = PREDICTIVE_TARGET_CONTRACTS.get(row_id)
            if contract is None:
                fail(f"{where}: predictive row lacks a fixed target contract")
            if contract["attainment"] != "all":
                fail(
                    f"{where}: fixed predictive target contract must require all"
                    " targets"
                )
            if tuple(targets) != contract["targets"]:
                fail(
                    f"{where}: frozen_targets must equal fixed scientific"
                    f" contract {list(contract['targets'])}"
                )
            if row["status"] == "attained":
                locked = [
                    frozen_target_by_id[target]["status"] in FROZEN_TARGET_STATUSES
                    for target in targets
                ]
                qualifies = all(locked)
                if not qualifies:
                    fail(
                        f"{where}: attained predictive row requires"
                        f" {contract['attainment']} targets in its fixed contract"
                        " to be frozen or locked"
                    )
        lane = row["lane_issue"]
        if not isinstance(lane, int) or isinstance(lane, bool):
            fail(f"{where}: lane_issue must be an integer")
        if lane not in LANE_ISSUES:
            fail(f"{where}: lane_issue {lane} is not a registered V3 lane")

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
        contract = ROW_PREMISE_CONTRACTS.get(row["id"])
        if contract is not None:
            for field, expected in contract.items():
                if tuple(row[field]) != expected:
                    fail(
                        f"{where}: {field} must equal the fixed "
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

    grouped_lanes = [lane for _, lanes in GROUPS for lane in lanes]
    if set(grouped_lanes) != set(LANE_ISSUES) or len(grouped_lanes) != len(LANE_ISSUES):
        fail("group table must cover every lane exactly once")
    predictive_ids = {row["id"] for row in rows if row["rung"] == "predictive"}
    if predictive_ids != set(PREDICTIVE_TARGET_CONTRACTS):
        fail("fixed predictive-target contracts must cover every predictive row")
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

    return rows


def _lane_link(lane: int) -> str:
    return f"[#{lane}]({REPO_URL}/issues/{lane})"


def _issue_link(number: int) -> str:
    return f"[issue #{number}]({REPO_URL}/issues/{number})"


def render(rows: list[dict]) -> str:
    lines: list[str] = []
    lines.append("# V3 Observation Ledger")
    lines.append("")
    lines.append(
        "Generated by `tools/build_observation_ledger.py` from"
        " `tracking/observation_ledger.json`; edit the JSON, then regenerate."
    )
    lines.append("")
    lines.append(
        f"One row per observation the architecture must reproduce. The closed"
        f" {_issue_link(726)} is the historical bootstrap; this committed"
        f" ledger is the maintained surface. Premise ids PR-01 through"
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
        "- **Predictive**: a pre-comparison frozen instrument owned by the"
        " relevant physics lane binds the row to future data with registered"
        " kill bands."
    )
    lines.append("")
    lines.append(
        "A structural claim whose premises are unregistered is invalid, an"
        " emergent claim without a preregistered instrument is invalid, and a"
        " predictive claim without its exact pre-comparison frozen target is"
        " invalid."
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
        " predictions. Predictive rows name their related FZ register targets;"
        " an attained predictive row requires the targets specified by its"
        " fixed row-level contract to be frozen or locked."
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
            "| Row | Observation | Rung | Status | Lane | Frozen targets | Premises |"
            " Open premises | Evidence | Boundary |"
        )
        lines.append(
            "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"
        )
        for row in group_rows:
            frozen_targets = (
                ", ".join(row.get("frozen_targets", []))
                if row.get("frozen_targets")
                else "none"
            )
            premises = ", ".join(row["premises"]) if row["premises"] else "none"
            open_premises = (
                ", ".join(row["open_premises"]) if row["open_premises"] else "none"
            )
            evidence = (
                ", ".join(f"`{path}`" for path in row["evidence"])
                if row["evidence"]
                else "none"
            )
            lines.append(
                f"| {row['id']} | {_cell(row['target'])} | {row['rung']} |"
                f" {row['status']} | {_lane_link(row['lane_issue'])} |"
                f" {frozen_targets} | {premises} |"
                f" {open_premises} | {evidence} |"
                f" {_cell(row['notes'])} |"
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
    lines.append("## Premise usage and missing attachments")
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
        f"Register rows that no declared target consumes are omitted"
        f" from this list. `Premises` records hypotheses used by the attained"
        f" portion; `Open premises` records registered hypotheses needed"
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
    surface = render(rows).encode("utf-8")
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
