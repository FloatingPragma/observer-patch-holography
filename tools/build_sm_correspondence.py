"""Build and validate the Standard Model Lagrangian correspondence table.

Issue #735's headline deliverable: one row per term of the textbook Standard
Model term-and-sector inventory, each classified by what OPH supplies. The
machine-readable table is ``tracking/sm_lagrangian_correspondence.json``.
This tool validates it fail-closed and renders
``docs/SM_LAGRANGIAN_CORRESPONDENCE.md``; ``--check`` fails when the
committed page differs from the render.

Fail-closed rules: row ids and term labels match the exact ordered nineteen-row
inventory; every
classification is one of the four enum values; a derived_conditional or
registered_premise row names at least one premise id and an absent row names
none; every premise id (in the premises list or cited inline in a boundary)
is a row of the canonical program-wide premise register; every evidence path
resolves to a committed file; a
derived_conditional or partial row cites at least one evidence path; and the
rendered prose carries no banned wording.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

import strict_json

ROOT = Path(__file__).resolve().parents[1]
TABLE_PATH = ROOT / "tracking" / "sm_lagrangian_correspondence.json"
SURFACE_PATH = ROOT / "docs" / "SM_LAGRANGIAN_CORRESPONDENCE.md"
PREMISE_REGISTER_PATH = ROOT / "tracking" / "premise_register.json"

SCHEMA = "oph.sm_lagrangian_correspondence.v2"
ISSUE = 735

CLASSIFICATIONS = (
    "derived_conditional",
    "partial",
    "registered_premise",
    "absent",
)
CLASS_LABELS = {
    "derived_conditional": "derived conditional",
    "partial": "partial",
    "registered_premise": "registered premise",
    "absent": "absent",
}

TOP_KEYS = {"schema", "issue", "policy", "premise_register_source", "rows"}
ROW_KEYS = {
    "id",
    "term",
    "classification",
    "premises",
    "open_premises",
    "evidence",
    "boundary",
}

PREMISE_REGISTER_SOURCE = "tracking/premise_register.json"

EXPECTED_TERMS = (
    "SU(3) gauge kinetic term: -(1/4) G^a_munu G_a^munu",
    "SU(2) gauge kinetic term: -(1/4) W^i_munu W_i^munu",
    "U(1) hypercharge kinetic term: -(1/4) B_munu B^munu",
    "Nonabelian gauge self-interactions: cubic g f A A dA and quartic g^2 f f A A A A structures from the SU(2) and SU(3) field strengths",
    "Quark doublet kinetic/covariant term: Qbar i gamma^mu D_mu Q, Q = (u_L, d_L), Y = 1/6",
    "Up-singlet kinetic/covariant term: ubar_R i gamma^mu D_mu u_R, Y = 2/3",
    "Down-singlet kinetic/covariant term: dbar_R i gamma^mu D_mu d_R, Y = -1/3",
    "Lepton doublet kinetic/covariant term: Lbar i gamma^mu D_mu L, L = (nu_L, e_L), Y = -1/2",
    "Electron-singlet kinetic/covariant term: ebar_R i gamma^mu D_mu e_R, Y = -1",
    "Right-handed neutrino stance: nu_R, no kinetic or mass term in the minimal Standard Model",
    "Higgs kinetic/covariant term: (D_mu H)^dagger (D^mu H), H a weak doublet with Y = 1/2",
    "Higgs potential: mu^2 H^dagger H + lambda (H^dagger H)^2",
    "Up-type Yukawa term: - Qbar Y_u Htilde u_R + h.c.",
    "Down-type Yukawa term: - Qbar Y_d H d_R + h.c.",
    "Lepton Yukawa term: - Lbar Y_e H e_R + h.c.",
    "Theta-QCD term: (theta g3^2 / 32 pi^2) G^a_munu Gtilde_a^munu",
    "Generation triplication: three copies of the chiral fermion content",
    "Gauge couplings: g1, g2, g3 as numerical parameters",
    "CKM/PMNS mixing structure: V_CKM in the charged quark current and U_PMNS in the lepton sector",
)

PR_TOKEN = re.compile(r"PR-\d{2}")
_H_WORD = "hon" + "est"
BANNED_WORDS = re.compile(
    r"\b("
    + "|".join(
        (
            _H_WORD,
            _H_WORD + "ly",
            _H_WORD + "y",
            "now",
            "currently",
            "previously",
            "already",
            "recently",
        )
    )
    + r")\b",
    re.IGNORECASE,
)
BANNED_FRAGMENTS = ("—",)


def fail(message: str) -> None:
    raise SystemExit(f"sm correspondence: {message}")


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


def canonical_premise_register() -> tuple[dict, ...]:
    data = load_json(PREMISE_REGISTER_PATH)
    rows = data.get("rows")
    if not isinstance(rows, list) or not rows:
        fail("canonical premise register must contain rows")
    projection: list[dict] = []
    seen: set[str] = set()
    for row in rows:
        if not isinstance(row, dict):
            fail("canonical premise register contains a malformed row")
        projected = {key: row.get(key) for key in ("id", "name", "type", "disposition")}
        premise_id = projected["id"]
        if not isinstance(premise_id, str) or premise_id in seen:
            fail("canonical premise register ids must be unique strings")
        seen.add(premise_id)
        projection.append(projected)
    return tuple(projection)


def check_prose(where: str, text: str) -> None:
    for fragment in BANNED_FRAGMENTS:
        if fragment in text:
            fail(f"{where}: banned fragment {fragment!r}")
    match = BANNED_WORDS.search(text)
    if match:
        fail(f"{where}: banned word {match.group(0)!r}")


def validate(data: dict) -> list[dict]:
    if not isinstance(data, dict) or set(data) != TOP_KEYS:
        fail(f"top-level keys must equal {sorted(TOP_KEYS)}")
    if data["schema"] != SCHEMA:
        fail(f"schema must equal {SCHEMA}")
    if data["issue"] != ISSUE:
        fail(f"issue must equal {ISSUE}")
    if not isinstance(data["policy"], str) or not data["policy"].strip():
        fail("policy must be a nonempty string")
    check_prose("policy", data["policy"])

    if data["premise_register_source"] != PREMISE_REGISTER_SOURCE:
        fail(f"premise_register_source must equal {PREMISE_REGISTER_SOURCE!r}")
    premise_register = canonical_premise_register()
    register_ids = {entry["id"] for entry in premise_register}

    rows = data["rows"]
    if not isinstance(rows, list) or len(rows) != len(EXPECTED_TERMS):
        fail(f"rows must contain exactly {len(EXPECTED_TERMS)} entries")
    for index, row in enumerate(rows):
        expected_id = f"SML-{index + 1:02d}"
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
        if row["id"] != expected_id:
            fail(f"{where}: id must equal {expected_id}")
        where = row["id"]
        term = row["term"]
        if term != EXPECTED_TERMS[index]:
            fail(f"{where}: term must equal the canonical inventory entry")
        check_prose(f"{where}.term", term)
        classification = row["classification"]
        if classification not in CLASSIFICATIONS:
            fail(
                f"{where}: classification {classification!r} is not one of "
                f"{CLASSIFICATIONS}"
            )
        premises = row["premises"]
        if (
            not isinstance(premises, list)
            or any(not isinstance(p, str) for p in premises)
            or len(premises) != len(set(premises))
        ):
            fail(f"{where}: premises must be a duplicate-free list of strings")
        unknown = set(premises) - register_ids
        if unknown:
            fail(f"{where}: unknown premise ids {sorted(unknown)}")
        if classification in ("derived_conditional", "registered_premise"):
            if not premises:
                fail(
                    f"{where}: a {classification} row must name at least one "
                    "premise register row"
                )
        if classification == "absent" and premises:
            fail(f"{where}: an absent row names no premises")
        open_premises = row["open_premises"]
        if (
            not isinstance(open_premises, list)
            or any(not isinstance(p, str) for p in open_premises)
            or len(open_premises) != len(set(open_premises))
        ):
            fail(f"{where}: open_premises must be a duplicate-free list of strings")
        open_unknown = set(open_premises) - register_ids
        if open_unknown:
            fail(f"{where}: unknown open premise ids {sorted(open_unknown)}")
        if set(premises) & set(open_premises):
            fail(f"{where}: consumed and open premise lists must be disjoint")
        evidence = row["evidence"]
        if (
            not isinstance(evidence, list)
            or any(not isinstance(p, str) for p in evidence)
            or len(evidence) != len(set(evidence))
        ):
            fail(f"{where}: evidence must be a duplicate-free list of paths")
        for path in evidence:
            if not (ROOT / path).is_file():
                fail(f"{where}: evidence path missing: {path}")
        if classification in ("derived_conditional", "partial") and not evidence:
            fail(
                f"{where}: a {classification} row must cite at least one evidence path"
            )
        boundary = row["boundary"]
        if (
            not isinstance(boundary, str)
            or not boundary.strip()
            or not boundary.endswith(".")
        ):
            fail(f"{where}: boundary must be a nonempty sentence ending in '.'")
        check_prose(f"{where}.boundary", boundary)
        inline_unknown = set(PR_TOKEN.findall(boundary)) - register_ids
        if inline_unknown:
            fail(
                f"{where}: boundary cites unknown premise ids {sorted(inline_unknown)}"
            )
    return rows


def referenced_register_rows(rows: list[dict]) -> list[dict]:
    ids: set[str] = set()
    for row in rows:
        ids.update(row["premises"])
        ids.update(row["open_premises"])
        ids.update(PR_TOKEN.findall(row["boundary"]))
    return [entry for entry in canonical_premise_register() if entry["id"] in ids]


def render(data: dict, rows: list[dict]) -> str:
    counts = {cls: 0 for cls in CLASSIFICATIONS}
    for row in rows:
        counts[row["classification"]] += 1

    lines: list[str] = []
    lines.append("# The Standard Model term-and-sector correspondence table")
    lines.append("")
    lines.append(
        "Generated by `tools/build_sm_correspondence.py` from"
        " `tracking/sm_lagrangian_correspondence.json`; edit the JSON, then"
        " regenerate. Issue #735 owns this surface."
    )
    lines.append("")
    lines.append(data["policy"])
    lines.append("")
    lines.append(
        f"Summary: {len(rows)} term-and-sector entries;"
        f" {counts['derived_conditional']} derived conditional,"
        f" {counts['partial']} partial,"
        f" {counts['registered_premise']} registered premise,"
        f" {counts['absent']} absent."
    )
    lines.append("")
    lines.append("| Row | Term or sector | Classification | Premises | Open premises |")
    lines.append("| --- | --- | --- | --- | --- |")
    for row in rows:
        name, expression = row["term"].split(": ", 1)
        premises = ", ".join(row["premises"]) if row["premises"] else "none"
        open_premises = (
            ", ".join(row["open_premises"]) if row["open_premises"] else "none"
        )
        lines.append(
            f"| {row['id']} | {name}: `{expression}` |"
            f" `{row['classification']}` | {premises} | {open_premises} |"
        )
    lines.append("")
    lines.append("## Term boundaries")
    lines.append("")
    lines.append(
        "One paragraph per row: exactly what is supplied, exactly what is"
        " not, and the committed evidence."
    )
    for row in rows:
        name = row["term"].split(": ", 1)[0]
        evidence = (
            ", ".join(f"`{path}`" for path in row["evidence"])
            if row["evidence"]
            else "none"
        )
        lines.append("")
        lines.append(
            f"**{row['id']}, {name}"
            f" ({CLASS_LABELS[row['classification']]}).**"
            f" {row['boundary']} Evidence: {evidence}."
        )
    lines.append("")
    lines.append("## Premise register key")
    lines.append("")
    lines.append(
        "Premise ids name rows of the program-wide premise register"
        " (issue #727). Rows referenced by this table:"
    )
    lines.append("")
    for entry in referenced_register_rows(rows):
        lines.append(
            f"- `{entry['id']}` {entry['name']}"
            f" (type {entry['type']}, disposition {entry['disposition']})"
        )
    lines.append("")
    lines.append(
        "The classification enum is closed: `derived_conditional` is an"
        " exact theorem shape under the named register rows,"
        " `registered_premise` is consumed as a declared register row,"
        " `partial` carries the attained structure stated in its boundary,"
        " and `absent` has no OPH structure. A row upgrade requires"
        " receipts, a premise list, and a regeneration of this surface."
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

    data = load_json(TABLE_PATH)
    rows = validate(data)
    surface = render(data, rows)
    if args.check:
        committed = SURFACE_PATH.read_bytes() if SURFACE_PATH.is_file() else b""
        if committed != surface.encode("utf-8"):
            print(
                "sm correspondence: docs/SM_LAGRANGIAN_CORRESPONDENCE.md is"
                " stale; run python tools/build_sm_correspondence.py",
                file=sys.stderr,
            )
            return 1
        print("sm correspondence: surface is current")
        return 0
    SURFACE_PATH.write_bytes(surface.encode("utf-8"))
    print(f"sm correspondence: wrote {SURFACE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
