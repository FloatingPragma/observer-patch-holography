"""Build and validate the premise register surface.

The machine-readable register is ``tracking/premise_register.json``. This tool
validates it fail-closed and renders ``docs/registers/PREMISE_REGISTER_V3.md``;
``--check`` fails when the committed page differs byte-for-byte from the
render.

The same file carries three further sections. ``sm_coverage`` (the Standard
Model term-and-sector correspondence) and ``gravity_rungs`` (the gravity
premise-elimination ladder) are validated by ``tools/build_sm_correspondence.py``
and ``tools/build_gravity_ladder.py`` and rendered as sections of the page
written here. ``model_supply`` (the model and theorem navigation crosswalk) is
validated and rendered to ``model/PREMISES.md`` by
``tools/build_model_theorem_map.py``; this tool checks only its shape.

Fail-closed rules: the row ids, names, types, and dispositions are fixed
program-wide and must equal the expected inventory exactly, in order; every
type comes from the six-value enum and every disposition from the three-value
enum; every consuming lane is a registered historical lane, duplicate-free and
ascending, rendered beside its live owner from ``tools/closed_lanes.py``;
every evidence path is an existing file in the repository and has exactly one
explicit evidence role; every row carries a discharge record for the current
executable implementation, with a status from the five-value enum, existing
evidence files and a note, and a status that asserts an executed artefact must carry at
least one non-prose evidence artefact; statements and notes are nonempty. The
rendered page is a
generated surface; edit the JSON, then regenerate.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import strict_json
from closed_lanes import is_known, owner_text

import build_gravity_ladder as gravity_tool
import build_model_theorem_map as model_map_tool
import build_sm_correspondence as sm_tool
import check_claim_registry as claim_tool

ROOT = Path(__file__).resolve().parents[1]
REGISTER_PATH = ROOT / "tracking" / "premise_register.json"
SURFACE_PATH = ROOT / "docs" / "registers" / "PREMISE_REGISTER_V3.md"

SCHEMA = "oph.premise_register.v3"
ISSUE = 727
ISSUE_URL = "https://github.com/FloatingPragma/observer-patch-holography/issues"

# Sections folded into the register, in file order, with the schema each keeps.
SECTION_SCHEMAS = {
    "model_supply": model_map_tool.SCHEMA,
    "sm_coverage": sm_tool.SCHEMA,
    "gravity_rungs": gravity_tool.SCHEMA,
}
TOP_KEYS = ("schema", "issue", "rows", *SECTION_SCHEMAS)

TYPES = (
    "structural_rule",
    "representation_choice",
    "selection_rule",
    "external_mathematics",
    "empirical_import",
    "numerical_input",
)
DISPOSITIONS = ("remove", "axiomatize", "import")
LANE_ISSUES = frozenset((*range(728, 738), 740, 742, 743, 744, 745, 750, 763))
EVIDENCE_ROLES = {
    "statement",
    "conditional_consumer",
    "no_go",
    "external_input",
}

ROW_KEYS = {
    "id",
    "name",
    "type",
    "disposition",
    "statement",
    "consuming_lanes",
    "evidence",
    "evidence_roles",
    "notes",
    "discharge",
}

# What the current executable implementation does with each premise.
DISCHARGE_STATUSES = (
    "derived_from_native_records",
    "constructed_on_declared_inputs",
    "obstructed",
    "declared",
    "not_addressed",
)
DISCHARGE_KEYS = {"status", "evidence", "note"}

# Discharge statuses that assert an executed artefact exists, not merely a
# declaration or a stated obligation. Their evidence therefore cannot be prose
# only: a note may describe the run, but it cannot substitute for it.
EXECUTION_ASSERTING_DISCHARGE = (
    "derived_from_native_records",
    "constructed_on_declared_inputs",
)

# Fixed program-wide inventory: (id, name, type, disposition), in order.
EXPECTED_ROWS = (
    ("PR-01", "confluent terminating repair contract", "structural_rule", "axiomatize"),
    (
        "PR-02",
        "algebra-state representation of public records",
        "representation_choice",
        "axiomatize",
    ),
    ("PR-03", "operational effect additivity", "selection_rule", "remove"),
    ("PR-04", "declared phase-sensitive effect", "structural_rule", "axiomatize"),
    ("PR-05", "counting path reference", "selection_rule", "remove"),
    ("PR-06", "real Legendre enrichment", "representation_choice", "remove"),
    (
        "PR-07",
        "declared repair law: common faithful reference and repaired-visible fibre",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-08",
        "refinement-uniform low-temperature control",
        "structural_rule",
        "remove",
    ),
    (
        "PR-09",
        "equal-face-weight and barycentric one-third port-dual rules",
        "selection_rule",
        "axiomatize",
    ),
    (
        "PR-10",
        "measure-to-metric and nearest-point repair rules",
        "selection_rule",
        "axiomatize",
    ),
    ("PR-11", "compact-Lie classification inputs", "external_mathematics", "import"),
    ("PR-12", "matter candidate grammar", "selection_rule", "remove"),
    ("PR-13", "Koide balance premise", "selection_rule", "remove"),
    ("PR-14", "hadronic transport packet", "empirical_import", "import"),
    ("PR-15", "clock and energy calibration anchors", "empirical_import", "import"),
    (
        "PR-16",
        "source-causal continuum and stable-causality certificate",
        "structural_rule",
        "remove",
    ),
    (
        "PR-17",
        "proposed fundamental numerical parameters P and N",
        "numerical_input",
        "import",
    ),
    (
        "PR-18",
        "spherically symmetric radial readout",
        "structural_rule",
        "remove",
    ),
    ("PR-19", "shell-flux normalization", "structural_rule", "remove"),
    (
        "PR-20",
        "equal source-counting measure on the sixty directed source seams",
        "selection_rule",
        "axiomatize",
    ),
    (
        "PR-21",
        "auxiliary oscillator lift of the seam-current generator",
        "representation_choice",
        "axiomatize",
    ),
    (
        "PR-22",
        "physical-frequency identification for the seam-current symbol",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-23",
        "typed modular flow supply",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-24",
        "half-sided inclusion and null translation supply",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-25",
        "null-stress tomography supply",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-26",
        "generalized-entropy stationarity supply",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-27",
        "small-ball geometry supply",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-28",
        "physical tower and calibration supply",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-29",
        "layered shell-cardinality law",
        "structural_rule",
        "remove",
    ),
    (
        "PR-30",
        "layered shell equidistribution",
        "structural_rule",
        "remove",
    ),
    (
        "PR-31",
        "layered steady sourcing",
        "structural_rule",
        "remove",
    ),
    (
        "PR-32",
        "Cartesian joint carrier of an observer pair",
        "representation_choice",
        "axiomatize",
    ),
    (
        "PR-33",
        "two-slot diamond region map",
        "selection_rule",
        "axiomatize",
    ),
    (
        "PR-34",
        "two-observer slot split",
        "representation_choice",
        "axiomatize",
    ),
    (
        "PR-35",
        "same-source loop-to-kernel identity",
        "structural_rule",
        "remove",
    ),
    (
        "PR-36",
        "family band-selection rules: single complete faithful band and operational cost order",
        "selection_rule",
        "remove",
    ),
    (
        "PR-37",
        "fine-structure closure map selection",
        "selection_rule",
        "remove",
    ),
    (
        "PR-38",
        "same-quantity bridge for the P closure",
        "structural_rule",
        "remove",
    ),
    (
        "PR-39",
        "physical Thomson readback contract",
        "structural_rule",
        "remove",
    ),
    (
        "PR-40",
        "payload-coherent anchor-gap premise",
        "selection_rule",
        "remove",
    ),
    (
        "PR-41",
        "electroweak sector selection premises",
        "selection_rule",
        "remove",
    ),
    (
        "PR-42",
        "horizon record central-charge identification",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-43",
        "pointwise-continuous public star-automorphism flow",
        "structural_rule",
        "remove",
    ),
    (
        "PR-44",
        "finite subsystem split and normalized local-operation interface",
        "representation_choice",
        "remove",
    ),
    (
        "PR-45",
        "source-selected real variational dynamics and symmetry data",
        "structural_rule",
        "remove",
    ),
    (
        "PR-46",
        "physical global-form selection attachments",
        "structural_rule",
        "remove",
    ),
    (
        "PR-47",
        "spacetime Spin and fermion-chirality attachment",
        "structural_rule",
        "remove",
    ),
    (
        "PR-48",
        "one-Higgs scalar carrier and kinetic-action attachment",
        "representation_choice",
        "remove",
    ),
    (
        "PR-49",
        "renormalizable Higgs potential and electroweak-breaking law",
        "structural_rule",
        "remove",
    ),
    (
        "PR-50",
        "Yukawa interaction-line and coupling-matrix attachment",
        "selection_rule",
        "remove",
    ),
    (
        "PR-51",
        "common screen-electroweak load and capacity bridge",
        "structural_rule",
        "remove",
    ),
    (
        "PR-52",
        "observer-to-physical-spacetime and causal attachment",
        "structural_rule",
        "remove",
    ),
    (
        "PR-53",
        "physical finite-carrier propagation, frame, and comparison attachment",
        "structural_rule",
        "remove",
    ),
    (
        "PR-54",
        "source gauge-field, current, and action attachment",
        "structural_rule",
        "remove",
    ),
    (
        "PR-55",
        "declared baryon- and lepton-number labels",
        "selection_rule",
        "remove",
    ),
    (
        "PR-56",
        "dimension-six baryon-operator eligibility and contraction grammar",
        "representation_choice",
        "remove",
    ),
    (
        "PR-57",
        "physical proton state and baryon-violating effective-action attachment",
        "structural_rule",
        "remove",
    ),
    (
        "PR-58",
        "time-indexed net evolution and time-slice generation interface",
        "structural_rule",
        "remove",
    ),
    (
        "PR-59",
        "declared five-mode exterior component and mask-selection grammar",
        "selection_rule",
        "remove",
    ),
    (
        "PR-60",
        "admissible scheduler class: mismatch non-increase clause",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-61",
        "admissible scheduler class: single-register locality clause",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-62",
        "admissible scheduler class: unit expected capacity clause",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-63",
        "two-register capacity-split reading",
        "representation_choice",
        "axiomatize",
    ),
    (
        "PR-64",
        "operational phase-instrument channel realization",
        "structural_rule",
        "remove",
    ),
    (
        "PR-65",
        "source-produced phase preparation and readback provenance",
        "structural_rule",
        "remove",
    ),
    (
        "PR-66",
        "declared discrete kinetic term of the temporal Maxwell action",
        "structural_rule",
        "remove",
    ),
    (
        "PR-67",
        "declared finite one-doublet electroweak grammar",
        "representation_choice",
        "remove",
    ),
    (
        "PR-68",
        "declared finite quartic breaking packet",
        "structural_rule",
        "remove",
    ),
    (
        "PR-69",
        "declared finite Yukawa-line coefficient",
        "selection_rule",
        "remove",
    ),
    (
        "PR-70",
        "source homogeneous edge position action",
        "selection_rule",
        "remove",
    ),
    (
        "PR-71",
        "continuous field same-operator sector attachment",
        "selection_rule",
        "remove",
    ),
    (
        "PR-72",
        "cofinal gluing into the compared regime",
        "structural_rule",
        "remove",
    ),
    (
        "PR-73",
        "declared finite carrier scale",
        "numerical_input",
        "axiomatize",
    ),
    (
        "PR-74",
        "coherent frame and boost transport",
        "structural_rule",
        "remove",
    ),
    (
        "PR-75",
        "frozen nuisance and coverage readout rules",
        "representation_choice",
        "axiomatize",
    ),
    (
        "PR-76",
        "source-derived positive exclusivity lower bound",
        "numerical_input",
        "remove",
    ),
    (
        "PR-77",
        "pathwise weak fairness for canonical adaptive repair",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-78",
        "work-conserving canonical adaptive choice",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-79",
        "bounded waste between canonical adaptive changes",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-80",
        "sufficient cumulative attempt budget",
        "numerical_input",
        "axiomatize",
    ),
    (
        "PR-81",
        "certified causal read support",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-82",
        "abstract semantic ancestry acyclicity rank",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-84",
        "pathwise weak fairness for fixed node attempts",
        "structural_rule",
        "axiomatize",
    ),
    (
        "PR-85",
        "bounded waste between fixed node changes",
        "structural_rule",
        "axiomatize",
    ),
)

# Retired ids keep their number: no later row may reuse one.
RETIRED_ROWS = {
    "PR-83": (
        "commit-populated event base",
        "retired when the finite causal carrier was derived from authenticated"
        " source order (PR-81, PR-82)",
    ),
}

DISPOSITION_MEANING = {
    "remove": (
        "A remove row is consumed as a declared premise until the executable "
        "specification supplies it and observers recover that supply from "
        "their records. The row leaves the register only through such a "
        "supply on the same domain; every composition citing it is "
        "conditional on it, and any bounded no-go attached to the row records "
        "exactly what the specification fails to supply."
    ),
    "axiomatize": (
        "An axiomatize row is a candidate for adoption into the executable "
        "specification. Adoption is a recorded register decision with "
        "rationale and a specification-wide dependency and countermodel "
        "review, never a silent edit; the row keeps its source-selection "
        "no-gos as the evidence of what the specification does not supply, "
        "and an adopted row stays on the register with its decision record."
    ),
    "import": (
        "An import row is a flagged external input that stays on the "
        "register: classical external mathematics, pinned empirical data "
        "packets, and the two numerical inputs P and N. An import is never "
        "presented as a derivation; replacing an import with a source-side "
        "result is a recorded disposition change."
    ),
}


def fail(message: str) -> None:
    raise SystemExit(f"premise register: {message}")


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


def validate(register: dict) -> list[dict]:
    if not isinstance(register, dict):
        fail("register must be an object")
    if register.get("schema") != SCHEMA:
        fail(f"schema must equal {SCHEMA}")
    if register.get("issue") != ISSUE:
        fail(f"issue must equal {ISSUE}")
    if set(register) != set(TOP_KEYS):
        fail(f"top-level keys must be exactly {', '.join(TOP_KEYS)}")
    rows = register.get("rows")
    if not isinstance(rows, list):
        fail("rows must be a list")
    if len(rows) != len(EXPECTED_ROWS):
        fail(f"rows must contain exactly {len(EXPECTED_ROWS)} entries")
    reused = sorted(
        row.get("id") for row in rows if isinstance(row, dict) and row.get("id") in RETIRED_ROWS
    )
    if reused:
        fail(f"retired premise ids must not be reused: {reused}")

    for index, row in enumerate(rows):
        expected_id, expected_name, expected_type, expected_disposition = EXPECTED_ROWS[
            index
        ]
        where = f"row {expected_id}"
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
            fail(f"{where}: id must equal {expected_id}, found {row['id']!r}")
        if row["name"] != expected_name:
            fail(f"{where}: name must equal {expected_name!r}")
        if row["type"] not in TYPES:
            fail(f"{where}: type {row['type']!r} is not in the enum {TYPES}")
        if row["type"] != expected_type:
            fail(f"{where}: type must equal {expected_type!r}")
        if row["disposition"] not in DISPOSITIONS:
            fail(
                f"{where}: disposition {row['disposition']!r} is not in the "
                f"enum {DISPOSITIONS}"
            )
        if row["disposition"] != expected_disposition:
            fail(f"{where}: disposition must equal {expected_disposition!r}")

        statement = row["statement"]
        if not isinstance(statement, str) or not statement.strip():
            fail(f"{where}: statement must be a nonempty string")
        notes = row["notes"]
        if not isinstance(notes, str) or not notes.strip():
            fail(f"{where}: notes must be a nonempty string")

        lanes = row["consuming_lanes"]
        if not isinstance(lanes, list) or not lanes:
            fail(f"{where}: consuming_lanes must be a nonempty list")
        for lane in lanes:
            if not isinstance(lane, int) or isinstance(lane, bool):
                fail(f"{where}: consuming lane {lane!r} must be an integer")
            if lane not in LANE_ISSUES or not is_known(lane):
                fail(f"{where}: consuming lane {lane} is not a registered historical lane")
        if lanes != sorted(set(lanes)):
            fail(f"{where}: consuming_lanes must be strictly ascending")

        evidence = row["evidence"]
        if not isinstance(evidence, list) or not evidence:
            fail(f"{where}: evidence must be a nonempty list")
        if len(evidence) != len(set(evidence)):
            fail(f"{where}: evidence paths must be duplicate-free")
        for path in evidence:
            if not isinstance(path, str) or not path.strip():
                fail(f"{where}: evidence entries must be nonempty strings")
            if not (ROOT / path).is_file():
                fail(f"{where}: evidence path missing or not a file: {path}")
        evidence_roles = row["evidence_roles"]
        if not isinstance(evidence_roles, dict):
            fail(f"{where}: evidence_roles must be an object")
        if set(evidence_roles) != set(evidence):
            missing = set(evidence) - set(evidence_roles)
            extra = set(evidence_roles) - set(evidence)
            fail(
                f"{where}: evidence_roles must have exact evidence-path parity "
                f"(missing {sorted(missing)}, extra {sorted(extra)})"
            )
        for path, role in evidence_roles.items():
            if role not in EVIDENCE_ROLES:
                fail(
                    f"{where}: evidence role for {path} must be one of "
                    f"{sorted(EVIDENCE_ROLES)}"
                )
        discharge = row["discharge"]
        if not isinstance(discharge, dict) or set(discharge) != DISCHARGE_KEYS:
            fail(f"{where}: discharge must have exactly the keys {sorted(DISCHARGE_KEYS)}")
        if discharge["status"] not in DISCHARGE_STATUSES:
            fail(f"{where}: discharge status must be one of {DISCHARGE_STATUSES}")
        if not isinstance(discharge["evidence"], list) or not discharge["evidence"]:
            fail(f"{where}: discharge evidence must be a nonempty list")
        for path in discharge["evidence"]:
            if not isinstance(path, str) or not (ROOT / path).is_file():
                fail(f"{where}: discharge evidence path missing or not a file: {path}")
        if discharge["status"] in EXECUTION_ASSERTING_DISCHARGE:
            media = []
            for path in discharge["evidence"]:
                medium = claim_tool.evidence_medium(path)
                if medium is None:
                    fail(
                        f"{where}: discharge evidence medium is unclassified; add "
                        f"the suffix to EVIDENCE_MEDIA_BY_SUFFIX: {path}"
                    )
                media.append(medium)
            if all(medium in claim_tool.PROSE_EVIDENCE_MEDIA for medium in media):
                fail(
                    f"{where}: discharge status {discharge['status']!r} asserts an "
                    f"executed artefact, so its evidence cannot be prose only: "
                    f"{discharge['evidence']}"
                )
        if not isinstance(discharge["note"], str) or not discharge["note"].strip():
            fail(f"{where}: discharge note must be a nonempty string")
    validate_sections(register)
    return rows


def validate_sections(register: dict) -> None:
    """Validate the folded sections; the section tools own their rules."""
    for section, schema in SECTION_SCHEMAS.items():
        data = register[section]
        if not isinstance(data, dict):
            fail(f"section {section} must be an object")
        if data.get("schema") != schema:
            fail(f"section {section}: schema must equal {schema}")
    sm_tool.validate(register["sm_coverage"])
    gravity_tool.validate(register["gravity_rungs"])


def lane_links(lanes: list[int]) -> str:
    """Render the live owner beside the historical consuming lanes."""
    return owner_text(lanes)


DISCHARGE_LABELS = {
    "derived_from_native_records": "derived from native records",
    "constructed_on_declared_inputs": "constructed on declared inputs",
    "obstructed": "obstructed",
    "declared": "declared input",
    "not_addressed": "not addressed",
}


def render(rows: list[dict], register: dict) -> str:
    """Render the page; ``register`` supplies the validated folded sections."""
    counts = {
        disposition: sum(1 for row in rows if row["disposition"] == disposition)
        for disposition in DISPOSITIONS
    }

    lines: list[str] = []
    lines.append("# OPH Premise Register")
    lines.append("")
    lines.append(
        "Generated by `tools/build_premise_register.py` from"
        " `tracking/premise_register.json`; edit the JSON, then regenerate."
        f" Established under the historical issue #{ISSUE} (closed);"
        " maintained as a scientific register. Consuming lanes are the"
        " historical V3 lanes; each is rendered beside the open issue that"
        " owns its remaining work."
    )
    lines.append("")
    lines.append(
        "The same source file carries three further sections. `sm_coverage`"
        " classifies each Standard Model term and sector by what OPH supplies,"
        " and `gravity_rungs` types each gravity interface with a premise"
        " packet, a theorem or countermodel, and a deletion test; both are"
        " rendered at the end of this page as the"
        " [Standard Model correspondence](#sm-coverage) and the"
        " [gravity premise ladder](#gravity-rungs)."
        " `model_supply` holds the model and theorem navigation crosswalk,"
        " which `tools/build_model_theorem_map.py` renders to"
        " [model/PREMISES.md](../../model/PREMISES.md)."
    )
    lines.append("")
    lines.append(
        "This register records the canonical PR-numbered premises, with their"
        " types, intended dispositions, consuming lanes, statements, and evidence"
        " paths. It is not an exhaustive list of the hypotheses of every result."
        " Result-specific assumptions and imported results are declared in the"
        " [claim registry](../../claims/claim_registry.yaml), with named assumptions"
        " defined in the [assumption dictionary](../../claims/assumption_dictionary.md)."
        " Claim-registry validation rejects unknown premise references, missing"
        " per-claim premise classifications, and undefined assumption names. An"
        " explicit non-consumer of these canonical rows can require its"
        " own result-specific assumptions. The full input budget must be read"
        " from each result and its dependencies, not from this register alone."
        " Dispositions such as `axiomatize` and `remove` describe intended"
        " treatment; they do not certify a completed adoption into the"
        " specification or a completed removal."
        " P and N are the only proposed fundamental free"
        " numerical parameters; measured calibrations, fits, targets, and"
        " transport payloads remain separate flagged empirical inputs and must"
        " be included in the input budget of the results that use them."
        " Registering a premise legitimizes its use inside conditional"
        " compositions; it does not make the premise true, derived, or"
        " physical. Disposition changes, including adoption into the"
        " specification, require an explicit scientific rationale, never a"
        " silent edit."
    )
    lines.append("")
    lines.append(
        "Each evidence path has one explicit primary role: `statement` locates"
        " the premise interface, `conditional_consumer` uses that premise in a"
        " conditional result, `no_go` proves a source-selection or derivation"
        " obstruction, and `external_input` identifies imported mathematics or"
        " data. These roles describe evidential function only; they create no"
        " implicit reverse-consumer edge or lane ownership."
    )
    lines.append("")
    lines.append(
        "Retired ids keep their number and are never reused: "
        + "; ".join(
            f"{rid} ({name}), {reason}" for rid, (name, reason) in RETIRED_ROWS.items()
        )
        + "."
    )
    lines.append("")
    lines.append(
        "The current-implementation column says what the executable"
        " construction in `model/` does with each premise: derives it from"
        " native records, constructs an instance on declared inputs, proves an"
        " obstruction, takes it as a declared input, or does not address it."
        " It certifies no physical premise."
    )
    lines.append("")
    lines.append(
        "| Row | Premise | Type | Disposition | Current implementation |"
        " Live owner (historical lanes) |"
    )
    lines.append("| --- | --- | --- | --- | --- | --- |")
    for row in rows:
        lines.append(
            f"| {row['id']} | {row['name']} | `{row['type']}` |"
            f" `{row['disposition']}` |"
            f" {DISCHARGE_LABELS[row['discharge']['status']]} |"
            f" {lane_links(row['consuming_lanes'])} |"
        )
    lines.append("")
    lines.append(
        f"Totals: {len(rows)} premises. {counts['remove']} remove,"
        f" {counts['axiomatize']} axiomatize, {counts['import']} import."
        " Current implementation: "
        + ", ".join(
            f"{sum(1 for row in rows if row['discharge']['status'] == status)}"
            f" {DISCHARGE_LABELS[status]}"
            for status in DISCHARGE_STATUSES
        )
        + "."
    )
    lines.append("")

    lines.append("## Row statements and evidence")
    lines.append("")
    for row in rows:
        lines.append(f"### {row['id']} {row['name']}")
        lines.append("")
        lines.append(row["statement"])
        lines.append("")
        lines.append(
            f"- Type `{row['type']}`; disposition `{row['disposition']}`;"
            f" consumed by {lane_links(row['consuming_lanes'])}."
        )
        evidence = ", ".join(
            f"`{path}` (`{row['evidence_roles'][path]}`)"
            for path in row["evidence"]
        )
        lines.append(f"- Evidence: {evidence}.")
        lines.append(f"- Disposition note: {row['notes']}")
        discharge = row["discharge"]
        discharge_evidence = ", ".join(f"`{path}`" for path in discharge["evidence"])
        lines.append(
            f"- Current implementation: {DISCHARGE_LABELS[discharge['status']]}."
            f" {discharge['note']} Evidence: {discharge_evidence}."
        )
        lines.append("")

    lines.append("## What the dispositions mean")
    lines.append("")
    for disposition in DISPOSITIONS:
        lines.append(f"### {disposition}")
        lines.append("")
        lines.append(DISPOSITION_MEANING[disposition])
        lines.append("")

    lines.append(
        "Register rows keep the source-selection no-gos proved against them"
        " as permanent evidence of what the executable specification does not"
        " supply."
    )
    sm_data = register["sm_coverage"]
    gravity_data = register["gravity_rungs"]
    lines.append("")
    lines.append(sm_tool.render(sm_data, sm_data["rows"]).rstrip("\n"))
    lines.append("")
    lines.append(gravity_tool.render(gravity_data, gravity_data["rungs"]).rstrip("\n"))
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="fail when the committed surface differs from the render",
    )
    args = parser.parse_args(argv)

    register = load_json(REGISTER_PATH)
    rows = validate(register)
    surface = render(rows, register).encode("utf-8")
    if args.check:
        committed = SURFACE_PATH.read_bytes() if SURFACE_PATH.is_file() else b""
        if committed != surface:
            print(
                "premise register: docs/registers/PREMISE_REGISTER_V3.md is stale; run "
                "python tools/build_premise_register.py",
                file=sys.stderr,
            )
            return 1
        print("premise register: surface is current")
        return 0
    SURFACE_PATH.write_bytes(surface)
    print(f"premise register: wrote {SURFACE_PATH.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
