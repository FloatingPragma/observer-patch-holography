#!/usr/bin/env python3
"""Typed closure verdicts for the discriminator producer strata.

Issues #642, #643, #645, and #646 are producer strata inside the issue #647
invariant-mining registry. The registry carries the provisional
source-identity inventory type with physical scoring sealed, and each
stratum's outcome is a typed verdict that consumes the registry and the
pinned source receipts without opening any comparison surface. This
producer emits the four verdicts:

* issue #642: the diagonal Z6 packet converts to an exact conditional kill
  packet. The registry carries the descent congruence, the charge
  corollaries, the flux refinement, and the electric line-class exclusion;
  this verdict adds the exact line-lattice arithmetic for all four untilted
  global forms on the Z6 charge torus, the magnetic-sector statement, and
  the per-form distinguishing observations, and types the physical
  line-operator attachment open.
* issue #643: the stratum verdict is superseded by the angular-sprint
  template receipt, which is pinned by bytes. The registered angular
  candidates stay the retained source template. The operative #643 record
  is the angular receipt's narrowed transfer decision: static base-port
  underdetermination, with the canonical degree-three band exactly
  identifiable from the level-one refinement vertex set per the pinned
  refinement-transfer receipt, and the refinement/repair sky transfer
  typed open rather than nonidentifiable.
* issue #645: the normalized overlap cross-spectrum exists as registered
  source structure, and no map from the declared source interface to an
  interferometer readout exists, so the comparison contract is correctly
  not frozen.
* issue #646: the no-go label is withdrawn and the combination table is
  retained. Every tested pole, width, residue, and asymmetry combination
  stays blocked by a pinned receipt (the no-pole-promotion boundary, plus
  the serialized unit verdict for physical-unit rows). The kinetic lane
  records the exact form dichotomy from the pinned selection receipt: the
  ad-invariant port-response pullback and the rank-fifteen matter trace
  carry Killing-relative su(2):su(3) ratios six and three halves, the
  matter-branch determinant statistic is frozen with exact cofactors, and
  the selection between the two forms is a named open source premise.

Every verdict pins its inputs by exact bytes, fails closed on drift, reads
no public measurement, and permits no comparison.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parent.parent
RUNTIME = HERE / "runtime"
REGISTRY_PATH = (
    REPO_ROOT / "code" / "invariant_mining" / "outputs" / "candidate_registry.json"
)
WZ_MANIFEST_PATH = (
    REPO_ROOT
    / "code"
    / "particles"
    / "calibration"
    / "wz_upstream_completion"
    / "INTEGRATION_MANIFEST.json"
)

VERDICTS = {
    642: "CONDITIONAL_KILL_PACKET__PHYSICAL_LINE_ATTACHMENT_OPEN",
    643: (
        "SUPERSEDED_BY_ANGULAR_SPRINT__"
        "SEE_ANGULAR_TEMPLATE_RECEIPT"
    ),
    645: (
        "NOT_EVALUABLE_NO_REGISTERED_SOURCE_TO_READOUT_MAP__"
        "FINITE_TWIST_SECTOR_SPECTRUM_ONLY"
    ),
    646: (
        "WITHDRAWN_NO_GO_LABEL__COMBINATION_TABLE_RETAINED__"
        "KINETIC_DICHOTOMY_RECORDED__SELECTION_PREMISE_OPEN"
    ),
}
CLOCK_UNIT_VERDICT_PATH = (
    REPO_ROOT / "code" / "a5_closure" / "manifests" / "clock_unit_verdict.json"
)
CLASSICAL_RECEIPT_PATH = (
    REPO_ROOT
    / "code"
    / "a5_closure"
    / "manifests"
    / "classical_realization_receipt.json"
)
FZ_REGISTER_PATH = REPO_ROOT / "claims" / "frozen_prediction_register.json"
ANGULAR_RECEIPT_PATH = (
    REPO_ROOT
    / "code"
    / "angular_sprint"
    / "runtime"
    / "angular_template_receipt.json"
)
REFINEMENT_RECEIPT_PATH = (
    REPO_ROOT
    / "code"
    / "angular_sprint"
    / "runtime"
    / "refinement_transfer_receipt.json"
)
KINETIC_SELECTION_RECEIPT_PATH = (
    REPO_ROOT
    / "code"
    / "angular_sprint"
    / "runtime"
    / "kinetic_form_selection_receipt.json"
)


class VerdictError(ValueError):
    """A stratum verdict refused to build."""


def require(condition: bool, message: str) -> None:
    if not condition:
        raise VerdictError(message)


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


def _pin(path: Path) -> dict[str, Any]:
    payload = path.read_bytes()
    return {
        "path": path.relative_to(REPO_ROOT).as_posix(),
        "bytes": len(payload),
        "sha256": tagged_sha256(payload),
    }


def _load_registry() -> dict[str, Any]:
    registry = json.loads(REGISTRY_PATH.read_text(encoding="utf-8"))
    require(
        registry.get("schema") == "oph.invariant_mining.candidate_registry.v1",
        "candidate registry schema drift",
    )
    require(
        registry.get("status")
        == "PROVISIONAL_SOURCE_IDENTITY_INVENTORY__PHYSICAL_SCORING_SEALED",
        "candidate registry status drift",
    )
    boundary = registry.get("scoring_boundary", {})
    require(
        boundary.get("physical_scoring_permitted") is False
        and boundary.get("comparison_access_permitted") is False,
        "candidate registry scoring boundary drift",
    )
    return registry


def _slot_candidates(registry: dict[str, Any], slot_id: str) -> list[dict[str, Any]]:
    rows = [
        candidate
        for candidate in registry["candidates"]
        if candidate["slot_id"] == slot_id
    ]
    require(bool(rows), f"no registered candidates for slot {slot_id}")
    return rows


def _candidate_digest(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "candidate_id": row["candidate_id"],
            "grammar_class": row["grammar_class"],
            "rank": row["rank"],
            "score": row["score"],
            "kill_rule": row["kill_rule"],
        }
        for row in sorted(rows, key=lambda row: row["rank"])
    ]


# ---------------------------------------------------------------------------
# Issue #642: exact dyonic line-lattice arithmetic on the Z6 charge torus
# ---------------------------------------------------------------------------


def _dirac_pairing(e1: int, m1: int, e2: int, m2: int) -> int:
    return (e1 * m2 - e2 * m1) % 6


def _all_subgroups_z6_squared() -> list[frozenset[tuple[int, int]]]:
    """Every subgroup of Z6 x Z6 by closure of two-element generation."""

    elements = [(e, m) for e in range(6) for m in range(6)]
    subgroups: set[frozenset[tuple[int, int]]] = set()
    for first in elements:
        for second in elements:
            group = {(0, 0)}
            frontier = [first, second]
            while frontier:
                candidate = frontier.pop()
                if candidate in group:
                    continue
                group.add(candidate)
                for member in list(group):
                    total = (
                        (candidate[0] + member[0]) % 6,
                        (candidate[1] + member[1]) % 6,
                    )
                    if total not in group:
                        frontier.append(total)
            subgroups.add(frozenset(group))
    return sorted(subgroups, key=lambda group: (len(group), sorted(group)))


def dyonic_line_lattice_certificate() -> dict[str, Any]:
    """Exact line-class arithmetic across all four untilted global forms.

    Line classes live on the Z6 x Z6 charge torus with the Dirac pairing
    ``<(e1, m1), (e2, m2)> = e1*m2 - e2*m1 mod 6``. For the quotient by a
    central subgroup of order ``g``, the untilted line lattice is the
    product of the electric classes trivial on the subgroup with the
    magnetic classes it supplies. The certificate enumerates every maximal
    isotropic subgroup of the torus by exhaustion, verifies there are
    exactly twelve, locates the four untilted global-form lattices among
    them, and records the distinguishing observation for each form.
    """

    maximal_isotropic = []
    for group in _all_subgroups_z6_squared():
        if len(group) != 6:
            continue
        if any(
            _dirac_pairing(*left, *right) != 0
            for left in group
            for right in group
        ):
            continue
        outside_all_fail = all(
            any(_dirac_pairing(e, m, *member) != 0 for member in group)
            for e in range(6)
            for m in range(6)
            if (e, m) not in group
        )
        if outside_all_fail:
            maximal_isotropic.append(sorted(group))

    untilted_forms = {
        "unquotiented": {
            "quotient_order": 1,
            "lattice": sorted((e, 0) for e in range(6)),
            "distinguishing_observation": (
                "every electric line class occurs and only the trivial "
                "flux class exists"
            ),
        },
        "z2_quotient": {
            "quotient_order": 2,
            "lattice": sorted((e, m) for e in (0, 2, 4) for m in (0, 3)),
            "distinguishing_observation": (
                "electric line classes of even center charge with two flux "
                "classes"
            ),
        },
        "z3_quotient": {
            "quotient_order": 3,
            "lattice": sorted((e, m) for e in (0, 3) for m in (0, 2, 4)),
            "distinguishing_observation": (
                "electric line classes of center charge zero or three with "
                "three flux classes"
            ),
        },
        "z6_quotient": {
            "quotient_order": 6,
            "lattice": sorted((0, m) for m in range(6)),
            "distinguishing_observation": (
                "only descended electric line classes with all six flux "
                "classes; the declared branch"
            ),
        },
    }
    lattices_found = {
        name: sorted(tuple(pair) for pair in row["lattice"]) in [
            [tuple(pair) for pair in group] for group in maximal_isotropic
        ]
        for name, row in untilted_forms.items()
    }
    quotient_lattice = [(0, m) for m in range(6)]
    magnetic_statement = {
        "z6_quotient_pure_magnetic_classes_allowed": sorted(
            m for (e, m) in quotient_lattice
        ),
        "statement": (
            "in the declared quotient branch every pure magnetic flux class "
            "is allowed and every electric line class outside the descended "
            "lattice is forbidden"
        ),
    }
    return {
        "charge_torus": "Z6 x Z6",
        "dirac_pairing": "e1*m2 - e2*m1 mod 6",
        "maximal_isotropic_count": len(maximal_isotropic),
        "untilted_global_forms": {
            name: {
                "quotient_order": row["quotient_order"],
                "lattice": [list(pair) for pair in row["lattice"]],
                "present_among_maximal_isotropic": lattices_found[name],
                "distinguishing_observation": row["distinguishing_observation"],
            }
            for name, row in untilted_forms.items()
        },
        "magnetic_sector": magnetic_statement,
        "declared_branch_lattice": [list(pair) for pair in quotient_lattice],
        "allowed_classes": len(quotient_lattice),
        "forbidden_classes": 36 - len(quotient_lattice),
        "theta_dependence": (
            "the remaining eight maximal isotropic subgroups are the "
            "theta-tilted variants of the quotient lattices; the tilt under "
            "a vacuum angle is a named open direction owned by the "
            "strong-CP issue, and this certificate fixes the untilted "
            "normalization"
        ),
    }


def build_642(registry: dict[str, Any]) -> dict[str, Any]:
    rows = _slot_candidates(registry, "z6_charge_line_congruences")
    lattice = dyonic_line_lattice_certificate()
    require(
        lattice["maximal_isotropic_count"] == 12
        and all(
            row["present_among_maximal_isotropic"]
            for row in lattice["untilted_global_forms"].values()
        )
        and lattice["allowed_classes"] == 6,
        "dyonic lattice certificate drift",
    )
    return {
        "schema": "oph.discriminator_stratum_verdict.v1",
        "issue": 642,
        "status": VERDICTS[642],
        "registered_kill_rules": _candidate_digest(rows),
        "dyonic_line_lattice": lattice,
        "kill_rule_summary": [
            "a physical state violating the descent congruence",
            "a free fractionally charged color singlet",
            "a color triplet outside the descending charge classes",
            "an electric line class outside the descended lattice",
            "a dyonic line class outside the maximal isotropic lattice",
            "a flux and line refinement with index other than six",
        ],
        "physical_attachment": (
            "OPEN: the line-operator and flux observables have no laboratory "
            "attachment map on the declared source interface; every kill "
            "rule is conditional on the named z6 global-form branch"
        ),
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
        },
        "reopen_condition": (
            "a physical line-operator or flux attachment map, at which point "
            "the registered kill rules become scorable through the issue "
            "#639 protocol"
        ),
    }


def build_643(registry: dict[str, Any]) -> dict[str, Any]:
    rows = _slot_candidates(registry, "a5_angular_rules")
    register = json.loads(FZ_REGISTER_PATH.read_text(encoding="utf-8"))
    fz02 = next(
        row
        for row in register.get("candidates", register.get("rows", []))
        if row.get("id") == "FZ-02"
    )
    require(
        "Status correction 2026-07-30" in fz02.get("content", ""),
        "FZ-02 frame-lock status correction is absent",
    )
    return {
        "schema": "oph.discriminator_stratum_verdict.v1",
        "issue": 643,
        "status": VERDICTS[643],
        "registered_source_template": _candidate_digest(rows),
        "frame_lock_disposition": {
            "clause": "FZ02-R03b",
            "disposition": (
                "retired from the scientific target through the append-only "
                "FZ-02 status correction; the clause reopens only with a "
                "source-derived screen-to-sky template map"
            ),
            "custody": "claims/frozen_prediction_register.json",
        },
        "decision_boundary": (
            "the declared source interface carries the A5 action on the "
            "port coefficient module as internal algebra; deciding a "
            "physical rotation-remnant reading requires a screen-to-sky map "
            "that no registered source observable supplies, so the question "
            "is not decidable on the declared interface"
        ),
        "transfer_typing": {
            "static_result": (
                "the twelve base-port values do not select the sky "
                "completion: the canonical degree-three band carries a "
                "four-dimensional evaluation kernel on the base ports, "
                "entirely in the odd sector"
            ),
            "refinement_result": (
                "the same band has a zero evaluation kernel on the "
                "level-one refinement vertex set, so the static ambiguity "
                "is a truncation property of the twelve-port readout, not "
                "of the refinement tower"
            ),
            "open_direction": (
                "the refinement/repair sky transfer: whether the physical "
                "readout exposes refined-port values is an open source "
                "premise owned by the repair and refinement law; any "
                "future nonidentifiability verdict must exhibit two "
                "source-admissible completions intertwining refinement, "
                "coarse-graining, and the repair semigroup"
            ),
            "producer": "code/angular_sprint/refinement_transfer_certificate.py",
        },
        "forced_sky_statistic": (
            "not constructible on the static base-port interface: a forced "
            "nonzero sky statistic requires the open refinement/repair "
            "transfer, so no comparison contract is frozen"
        ),
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
        },
        "reopen_condition": (
            "a source-derived screen-to-sky template map, for which the "
            "level-one refinement identifiability of the canonical band is "
            "the registered frontier; the registered angular candidates "
            "are the retained source template for that event"
        ),
    }


def build_645(registry: dict[str, Any]) -> dict[str, Any]:
    rows = _slot_candidates(registry, "observer_overlap_cross_spectra")
    return {
        "schema": "oph.discriminator_stratum_verdict.v1",
        "issue": 645,
        "status": VERDICTS[645],
        "registered_source_spectrum": _candidate_digest(rows),
        "pinned_spectrum_receipt": (
            "code/a5_closure/manifests/classical_realization_receipt.json"
        ),
        "feasibility_result": (
            "the normalized cross-spectrum exists as registered source "
            "structure: the unique flat twist sector, the conjugate-sector "
            "symmetry, and the ladder stability, all bound to the pinned "
            "classical realization receipt"
        ),
        "readout_boundary": (
            "no map from the declared source interface to an interferometer "
            "readout exists: the seam physicalization and laboratory "
            "attachment are open producer obligations, so the comparison "
            "contract is correctly not frozen"
        ),
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
        },
        "reopen_condition": (
            "a physical seam-readout attachment; the registered spectrum "
            "rows are the retained source side of that future contract"
        ),
    }


def build_646(registry: dict[str, Any]) -> dict[str, Any]:
    rows = _slot_candidates(registry, "wz_scale_free_response")
    manifest = json.loads(WZ_MANIFEST_PATH.read_text(encoding="utf-8"))
    require(
        manifest.get("promotion_allowed") is False,
        "wz integration manifest promotion drift",
    )
    clock_unit = json.loads(CLOCK_UNIT_VERDICT_PATH.read_text(encoding="utf-8"))
    units_status = clock_unit.get("status", clock_unit.get("verdict"))
    require(
        isinstance(units_status, str)
        and "PHYSICAL_UNITS_NOT_EVALUABLE" in units_status,
        "clock-unit verdict drift",
    )
    tested_combinations = [
        {
            "combination": "M_W / M_Z",
            "dimensionless": True,
            "blocked_by": "no_pole_promotion",
            "receipt": (
                "the upstream stack forbids promoting any pole value, so "
                "neither mass enters a source-evaluable ratio"
            ),
        },
        {
            "combination": "Gamma_W / M_W",
            "dimensionless": True,
            "blocked_by": "no_pole_promotion",
            "receipt": "no source-evaluable width or pole value exists",
        },
        {
            "combination": "Gamma_Z / M_Z",
            "dimensionless": True,
            "blocked_by": "no_pole_promotion",
            "receipt": "no source-evaluable width or pole value exists",
        },
        {
            "combination": "pole-residue ratios",
            "dimensionless": True,
            "blocked_by": "no_pole_promotion",
            "receipt": (
                "the boundary diagnostic keeps every residue promotion flag "
                "false"
            ),
        },
        {
            "combination": "asymmetry combinations",
            "dimensionless": True,
            "blocked_by": "no_pole_promotion",
            "receipt": (
                "no source-evaluable coupling asymmetry exists without the "
                "pole and current attachments"
            ),
        },
        {
            "combination": "any physical-unit member of the vector",
            "dimensionless": False,
            "blocked_by": "physical_units_not_evaluable",
            "receipt": units_status,
        },
    ]
    return {
        "schema": "oph.discriminator_stratum_verdict.v1",
        "issue": 646,
        "status": VERDICTS[646],
        "tested_combinations": tested_combinations,
        "non_evaluability": (
            "every tested combination is blocked by a pinned receipt: the "
            "dimensionless rows by the upstream no-pole-promotion boundary, "
            "and the physical-unit rows additionally by the serialized "
            "local-domain unit verdict; the dimensionless rows stay "
            "evaluable in principle and reopen with a pole producer"
        ),
        "nuisance_mapping": (
            "the frozen mining registry carries basis_and_coordinate_choice, "
            "common_dimensionful_scale, physical_sector_choice, "
            "source_admissible_completion, and "
            "detector_transfer_calibration; the clock, scheme, and "
            "threshold directions live in the upstream W/Z stack vocabulary "
            "and stay outside the source cone until a pole producer exists"
        ),
        "adjacent_registered_response_identities": _candidate_digest(rows),
        "adjacency_note": (
            "the registered rows are band-cost, channel, and kinetic-form "
            "identities of the response spectrum, adjacent to and distinct "
            "from the tested pole vector"
        ),
        "kinetic_form_lane": {
            "ad_invariance": (
                "the port-current pairing is ad-invariant; the earlier "
                "su(3) non-invariance wording is superseded by the pinned "
                "selection receipt"
            ),
            "killing_relative_dichotomy": {
                "port_response": ["1", "1/6"],
                "matter_trace": ["1", "2/3"],
                "ratios": ["6", "3/2"],
            },
            "frozen_matter_branch_statistic": {
                "kinetic_column": ["10/3", "2", "2"],
                "beta_column": ["41/6", "-19/6", "-7"],
                "exact_cofactors": ["-23/3", "37", "-218/9"],
                "integer_zero_locus": "69 x1 - 333 x2 + 218 x3 = 0",
            },
            "selection_premise": (
                "which invariant form the repair dynamics selects as the "
                "physical kinetic action is a named open source premise; "
                "the frozen determinant scores only through the issue-639 "
                "custody surface after a certified kinetic-action bridge "
                "selects the matter trace"
            ),
            "producer": (
                "code/angular_sprint/kinetic_form_selection_certificate.py"
            ),
            "lean_check": "Lean/Screen/KineticFormDichotomy.lean",
        },
        "surviving_direction_owners": {
            "physical common-load attachment": 631,
            "field and operator census": 632,
            "matching and schemes": 32,
            "finite-source EFT bridge": 635,
            "family attachment": 569,
            "integrated action": 630,
        },
        "wz_manifest_status": manifest.get("scientific_status"),
        "comparison_boundary": {
            "public_measurement_read": False,
            "comparison_permitted": False,
        },
        "reopen_condition": (
            "a positive unit and pole producer chain, at which point the "
            "minimal surviving source inputs are the registered scale-free "
            "response identities"
        ),
    }


BUILDERS = {642: build_642, 643: build_643, 645: build_645, 646: build_646}


PIN_PATHS = {
    642: (REGISTRY_PATH,),
    643: (
        REGISTRY_PATH,
        FZ_REGISTER_PATH,
        ANGULAR_RECEIPT_PATH,
        REFINEMENT_RECEIPT_PATH,
    ),
    645: (REGISTRY_PATH, CLASSICAL_RECEIPT_PATH),
    646: (
        REGISTRY_PATH,
        WZ_MANIFEST_PATH,
        CLOCK_UNIT_VERDICT_PATH,
        KINETIC_SELECTION_RECEIPT_PATH,
    ),
}


def build_all() -> dict[int, dict[str, Any]]:
    registry = _load_registry()
    verdicts = {}
    for issue, builder in sorted(BUILDERS.items()):
        verdict = builder(registry)
        verdict["parent_pins"] = [_pin(path) for path in PIN_PATHS[issue]]
        verdict["verdict_sha256"] = tagged_sha256(canonical_json_bytes(verdict))
        verdicts[issue] = verdict
    return verdicts


def write_runtime() -> list[Path]:
    RUNTIME.mkdir(exist_ok=True)
    paths = []
    for issue, verdict in build_all().items():
        path = RUNTIME / f"stratum_verdict_{issue}.json"
        path.write_bytes(canonical_json_bytes(verdict))
        paths.append(path)
    return paths


def verify_runtime() -> None:
    for issue, verdict in build_all().items():
        path = RUNTIME / f"stratum_verdict_{issue}.json"
        if path.read_bytes() != canonical_json_bytes(verdict):
            raise SystemExit(f"stratum verdict {issue} is stale")
        committed = json.loads(path.read_text(encoding="utf-8"))
        if committed["status"] != VERDICTS[issue]:
            raise SystemExit(f"stratum verdict {issue} status drift")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    if args.write:
        for path in write_runtime():
            print(path)
    if args.verify:
        verify_runtime()
        print("STRATUM_VERDICTS_VALID")
    if not args.write and not args.verify:
        print(
            canonical_json_bytes(
                {issue: verdict["status"] for issue, verdict in build_all().items()}
            ).decode("ascii"),
            end="",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
