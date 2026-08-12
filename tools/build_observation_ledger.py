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
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LEDGER_PATH = ROOT / "tracking" / "observation_ledger.json"
SURFACE_PATH = ROOT / "docs" / "OBSERVATION_LEDGER_V3.md"

SCHEMA = "oph.observation_ledger.v1"
ISSUE = 726
REPO_URL = "https://github.com/muellerberndt/reverse-engineering-reality"

RUNGS = ("structural", "emergent", "predictive")
STATUSES = ("attained", "partial", "owed")
LANE_MIN = 728
LANE_MAX = 738

ID_PATTERN = re.compile(r"^OL-[A-I][1-9]$")
PREMISE_PATTERN = re.compile(r"^PR-\d{2}$")
BANNED_CHARACTERS = ("—", "–")

ROW_KEYS = {
    "id",
    "target",
    "rung",
    "status",
    "lane_issue",
    "premises",
    "evidence",
    "notes",
}

# Fixed program-wide premise register seed (issue #727): id -> (name, type,
# disposition). The register artifact under #727 is the authority; this copy
# exists so the ledger fails closed on ids the register does not carry.
PREMISE_REGISTER = {
    "PR-01": (
        "confluent terminating repair contract",
        "structural_rule",
        "axiomatize",
    ),
    "PR-02": (
        "algebra-state representation of public records",
        "representation_choice",
        "axiomatize",
    ),
    "PR-03": ("operational effect additivity", "selection_rule", "remove"),
    "PR-04": ("phase-operation instrument", "structural_rule", "remove"),
    "PR-05": ("counting path reference", "selection_rule", "remove"),
    "PR-06": ("real Legendre enrichment", "representation_choice", "remove"),
    "PR-07": (
        "declared repair law: common faithful reference and repaired-visible fibre",
        "structural_rule",
        "axiomatize",
    ),
    "PR-08": (
        "refinement-uniform low-temperature control",
        "structural_rule",
        "remove",
    ),
    "PR-09": (
        "equal-face-weight and barycentric one-third port-dual rules",
        "selection_rule",
        "axiomatize",
    ),
    "PR-10": (
        "measure-to-metric and nearest-point repair rules",
        "selection_rule",
        "axiomatize",
    ),
    "PR-11": ("compact-Lie classification inputs", "external_mathematics", "import"),
    "PR-12": ("matter candidate grammar", "selection_rule", "remove"),
    "PR-13": ("Koide balance premise", "selection_rule", "remove"),
    "PR-14": ("hadronic transport packet", "empirical_import", "import"),
    "PR-15": ("clock and energy calibration anchors", "empirical_import", "import"),
    "PR-16": (
        "stable causality and open-image conditions",
        "structural_rule",
        "remove",
    ),
    "PR-17": ("numerical inputs P and N", "numerical_input", "import"),
    "PR-18": (
        "spherically symmetric radial readout",
        "structural_rule",
        "remove",
    ),
    "PR-19": ("shell-flux normalization", "structural_rule", "remove"),
}

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
)


def fail(message: str) -> None:
    raise SystemExit(f"observation ledger: {message}")


def load_ledger(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        fail(f"missing input {path}")
    except json.JSONDecodeError as error:
        fail(f"invalid JSON in {path}: {error}")
    raise AssertionError("unreachable")


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

    seen_ids: set[str] = set()
    for index, row in enumerate(rows):
        where = f"rows[{index}]"
        if not isinstance(row, dict):
            fail(f"{where}: row must be an object")
        if set(row) != ROW_KEYS:
            missing = ROW_KEYS - set(row)
            extra = set(row) - ROW_KEYS
            fail(
                f"{where}: keys mismatch "
                f"(missing {sorted(missing)}, extra {sorted(extra)})"
            )
        row_id = row["id"]
        if not isinstance(row_id, str) or not ID_PATTERN.match(row_id):
            fail(f"{where}: id must match OL-<A..I><digit>")
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

        lane = row["lane_issue"]
        if not isinstance(lane, int) or isinstance(lane, bool):
            fail(f"{where}: lane_issue must be an integer")
        if not LANE_MIN <= lane <= LANE_MAX:
            fail(f"{where}: lane_issue must lie in {LANE_MIN}..{LANE_MAX}")

        premises = row["premises"]
        if not isinstance(premises, list):
            fail(f"{where}: premises must be a list")
        if len(premises) != len(set(premises)):
            fail(f"{where}: premises must be duplicate-free")
        for premise in premises:
            if not isinstance(premise, str) or not PREMISE_PATTERN.match(premise):
                fail(f"{where}: premise ids must match PR-<two digits>")
            if premise not in PREMISE_REGISTER:
                fail(f"{where}: premise {premise} is not on the seed register")
        if row["status"] == "owed" and premises:
            fail(f"{where}: an owed row carries no consumed premises")

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
    if sorted(grouped_lanes) != list(range(LANE_MIN, LANE_MAX + 1)):
        fail("group table must cover every lane exactly once")
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
        f"One row per observation the architecture must reproduce, tracked for"
        f" {_issue_link(726)}. Premise ids PR-01 through PR-19 name rows of the"
        f" premise register ({_issue_link(727)}), the audited anti-cheating"
        f" surface; each row lists the register rows its current status"
        f" consumes. The row set, rungs, and lane assignments follow completion"
        f" plan V3 (`plan/COMPLETION_PLAN_V3.md` in the oph-meta planning"
        f" workspace). The ledger records adequacy status. It is not itself"
        f" evidence: a promotion requires the owning lane's receipts and"
        f" survives independent audit before landing here."
    )
    lines.append("")
    lines.append("## Adequacy rungs")
    lines.append("")
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
            "| Row | Observation | Rung | Status | Lane | Premises |"
            " Evidence | Boundary |"
        )
        lines.append("| --- | --- | --- | --- | --- | --- | --- | --- |")
        for row in group_rows:
            premises = ", ".join(row["premises"]) if row["premises"] else "none"
            evidence = (
                ", ".join(f"`{path}`" for path in row["evidence"])
                if row["evidence"]
                else "none"
            )
            lines.append(
                f"| {row['id']} | {row['target']} | {row['rung']} |"
                f" {row['status']} | {_lane_link(row['lane_issue'])} |"
                f" {premises} | {evidence} | {row['notes']} |"
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
    lines.append("## Premises consumed")
    lines.append("")
    consumed = sorted({premise for row in rows for premise in row["premises"]})
    for premise in consumed:
        name, kind, disposition = PREMISE_REGISTER[premise]
        lines.append(f"- **{premise}** {name} ({kind}, {disposition})")
    lines.append("")
    lines.append(
        f"Register rows that no current status consumes are omitted from this"
        f" list; owed compositions name them in their boundary sentences where"
        f" they apply. The full register with types and dispositions lives"
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
