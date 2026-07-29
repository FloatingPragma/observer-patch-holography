#!/usr/bin/env python3
"""Build the bounded issue-506 alpha/HVP accounting verdict.

This packet does not implement an independent HVP calculation.  It checks
one recorded, tabulated KNT19 accounting row against byte-pinned source
artifacts and recomputes its arithmetic from primitive fields.  The raw
dispersive, independent-code, and lattice-HVP classes have no frozen
repository ingests and are therefore not evaluable.

The comparison data were known before this packet was written.  The packet
is retrospective, is not a prospective freeze, and emits no OPH prediction
of the physical fine-structure constant.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from decimal import Decimal, getcontext
from pathlib import Path
from typing import Any

PACKET_DIR = Path(__file__).resolve().parent
REPO_ROOT = PACKET_DIR.parents[2]

PROTOCOL_PATH = PACKET_DIR / "data" / "accounting_protocol_v2.json"
OUT_PATH = PACKET_DIR / "outputs" / "alpha_hvp_class_verdict.json"

SCHEMA = "oph.alpha_hvp_class_verdict.v2"
PROTOCOL_SCHEMA = "oph.alpha_hvp_accounting_protocol.v2"
ISSUE = 506
SHA256_RE = re.compile(r"^sha256:[0-9a-f]{64}$")

ENDPOINT_REL = "code/P_derivation/runtime/empirical_thomson_endpoint_current.json"
BRIDGE_REL = "code/P_derivation/runtime/anchor_scheme_bridge_current.json"
PAYLOAD_REL = "code/particles/runs/hadron/empirical_ee_hadronic_spectral_measure.json"
MEASURE_REL = (
    "code/particles/runs/hadron/empirical_ward_projected_spectral_measure.json"
)

REQUIRED_PIN_PATHS = {
    ENDPOINT_REL,
    BRIDGE_REL,
    PAYLOAD_REL,
    MEASURE_REL,
    "code/particles/hadron/ingest_empirical_ee_hadrons.py",
    "code/particles/hadron/derive_empirical_ward_projected_spectral_measure.py",
    "code/particles/hadron/empirical_ee_hadrons_sources.yaml",
    "code/P_derivation/empirical_thomson_endpoint.py",
    "code/P_derivation/anchor_scheme_bridge.py",
}
UNAVAILABLE_CLASSES = ("raw_dispersive", "independent_code", "lattice_hvp")

getcontext().prec = 90
ONE = Decimal(1)
REFERENCE_TOLERANCE = Decimal("1e-12")


class VerdictInputError(ValueError):
    """A protocol or pinned-source input failed closed."""


def sha256_bytes(data: bytes) -> str:
    return "sha256:" + hashlib.sha256(data).hexdigest()


def canonical_json(value: Any) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n"


def decimal(value: Any) -> Decimal:
    return Decimal(str(value))


def gate(passed: bool, detail: str) -> dict[str, Any]:
    return {
        "status": "PASS" if passed else "FAIL",
        "passed": bool(passed),
        "detail": detail,
    }


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise VerdictInputError(message)


def load_protocol(path: Path = PROTOCOL_PATH) -> dict[str, Any]:
    protocol = json.loads(path.read_text(encoding="utf-8"))
    _require(protocol.get("schema") == PROTOCOL_SCHEMA, "protocol schema mismatch")
    _require(protocol.get("issue") == ISSUE, "protocol issue mismatch")

    scope = protocol.get("scope", {})
    _require(
        scope.get("comparison_timing") == "retrospective",
        "comparison timing must be retrospective",
    )
    _require(
        scope.get("prospective_freeze") is False,
        "issue-506 packet cannot claim a prospective freeze",
    )
    _require(
        scope.get("independent_hvp_implementation") is False,
        "issue-506 packet cannot claim an independent HVP implementation",
    )
    _require(
        scope.get("physical_alpha_prediction") is False,
        "issue-506 packet cannot emit a physical alpha prediction",
    )

    classes = protocol.get("classes")
    _require(
        isinstance(classes, dict)
        and set(classes) == {"tabulated_dispersive", *UNAVAILABLE_CLASSES},
        "protocol class matrix mismatch",
    )
    _require(
        classes["tabulated_dispersive"].get("evaluation_mode")
        == "recorded_accounting_replay",
        "tabulated row must be a recorded accounting replay",
    )
    for name in UNAVAILABLE_CLASSES:
        row = classes[name]
        _require(
            row.get("evaluation_mode") == "not_evaluable_without_frozen_ingest",
            f"{name} must fail closed without a frozen ingest",
        )
        _require(
            isinstance(row.get("ingest_requirement"), str)
            and bool(row["ingest_requirement"].strip()),
            f"{name} needs an exact ingest requirement",
        )

    snapshot = protocol.get("preexisting_source_snapshot", {})
    commit = snapshot.get("commit")
    _require(
        isinstance(commit, str) and re.fullmatch(r"[0-9a-f]{40}", commit) is not None,
        "source snapshot requires a full commit",
    )
    pins = snapshot.get("sha256")
    _require(
        isinstance(pins, dict) and set(pins) == REQUIRED_PIN_PATHS,
        "source snapshot pin set mismatch",
    )
    for relative_path, digest in pins.items():
        _require(
            isinstance(relative_path, str)
            and not Path(relative_path).is_absolute()
            and ".." not in Path(relative_path).parts,
            f"invalid pinned path {relative_path!r}",
        )
        _require(
            isinstance(digest, str) and SHA256_RE.fullmatch(digest) is not None,
            f"invalid SHA-256 pin for {relative_path}",
        )
    return protocol


def verify_exact_pins(
    protocol: dict[str, Any], root: Path = REPO_ROOT
) -> dict[str, bytes]:
    """Load every declared source and reject any byte drift."""

    pinned: dict[str, bytes] = {}
    for relative_path, expected in protocol["preexisting_source_snapshot"][
        "sha256"
    ].items():
        path = root / relative_path
        _require(path.is_file(), f"missing pinned source {relative_path}")
        raw = path.read_bytes()
        actual = sha256_bytes(raw)
        _require(
            actual == expected,
            f"pinned source drift for {relative_path}: {actual} != {expected}",
        )
        pinned[relative_path] = raw
    return pinned


def verify_snapshot_commit(
    protocol: dict[str, Any], root: Path = REPO_ROOT
) -> dict[str, bytes]:
    """Verify that the declared snapshot is ancestral and contains every pin.

    A syntactically valid commit identifier is not evidence that the declared
    files existed with the declared bytes before this audit.  Resolve the
    commit in the repository, require it to be an ancestor of ``HEAD``, and
    hash each blob directly from that commit.
    """

    commit = protocol["preexisting_source_snapshot"]["commit"]

    def run_git(arguments: list[str]) -> subprocess.CompletedProcess[bytes]:
        try:
            return subprocess.run(
                ["git", *arguments],
                cwd=root,
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        except OSError as error:
            raise VerdictInputError(f"cannot execute git: {error}") from error

    resolved = run_git(["cat-file", "-e", f"{commit}^{{commit}}"])
    _require(
        resolved.returncode == 0,
        f"source snapshot commit does not resolve: {commit}",
    )
    ancestral = run_git(["merge-base", "--is-ancestor", commit, "HEAD"])
    _require(
        ancestral.returncode == 0,
        f"source snapshot commit is not an ancestor of HEAD: {commit}",
    )

    blobs: dict[str, bytes] = {}
    for relative_path, expected in protocol["preexisting_source_snapshot"][
        "sha256"
    ].items():
        blob = run_git(["show", f"{commit}:{relative_path}"])
        _require(
            blob.returncode == 0,
            f"source snapshot lacks pinned blob {relative_path} at {commit}",
        )
        actual = sha256_bytes(blob.stdout)
        _require(
            actual == expected,
            (
                f"snapshot blob drift for {relative_path} at {commit}: "
                f"{actual} != {expected}"
            ),
        )
        blobs[relative_path] = blob.stdout
    return blobs


def _load_pinned_json(pinned: dict[str, bytes], relative_path: str) -> dict[str, Any]:
    try:
        value = json.loads(pinned[relative_path].decode("utf-8"))
    except (KeyError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise VerdictInputError(
            f"invalid pinned JSON source {relative_path}: {error}"
        ) from error
    _require(isinstance(value, dict), f"{relative_path} must contain an object")
    return value


def evaluate_tabulated_accounting_replay(
    endpoint: dict[str, Any],
    bridge: dict[str, Any],
    payload: dict[str, Any],
    measure: dict[str, Any],
    protocol: dict[str, Any],
) -> dict[str, Any]:
    """Recompute the recorded accounting row without claiming independence."""

    endpoint_inputs = endpoint["inputs"]
    endpoint_values = endpoint["endpoint"]
    compare = endpoint["compare_only"]
    anchor = decimal(endpoint_inputs["source_anchor_inv_alpha_MZ"])
    lepton = decimal(endpoint_inputs["lepton_transport_delta_inv_alpha"])
    delta_had = decimal(endpoint_inputs["delta_alpha_had_5_MZ"])
    delta_had_unc = decimal(endpoint_inputs["delta_alpha_had_5_MZ_uncertainty"])
    codata = decimal(compare["codata_alpha_inv"])

    a_l = anchor + lepton
    endpoint_central = a_l / (ONE - delta_had)
    endpoint_lower = a_l / (ONE - (delta_had - delta_had_unc))
    endpoint_upper = a_l / (ONE - (delta_had + delta_had_unc))
    stored_endpoint_central = decimal(endpoint_values["alpha_inv_central"])
    stored_endpoint_interval = [
        decimal(value) for value in endpoint_values["alpha_inv_interval"]
    ]

    gap_lower = codata * (ONE - (delta_had + delta_had_unc)) - lepton - anchor
    gap_upper = codata * (ONE - (delta_had - delta_had_unc)) - lepton - anchor
    stored_gap = [
        decimal(value) for value in compare["same_scheme_anchor_gap_interval_inv_alpha"]
    ]
    protocol_gap = [
        decimal(value)
        for value in protocol["accounting_rule"]["same_scheme_gap_interval"]
    ]

    reference = bridge["reference_decomposition_compare_only"]
    reference_alpha_recomputed = decimal(reference["alpha_inv_0"]) * (
        ONE
        - decimal(reference["Delta_lep"])
        - decimal(reference["Delta_had5"])
        - decimal(reference["Delta_top"])
    )
    reference_alpha_stored = decimal(reference["alpha_inv_mz_phys_on_shell"])
    bridge_anchor = decimal(bridge["anchor_provenance"]["a0_oph"])
    deficit_recomputed = reference_alpha_recomputed - bridge_anchor
    deficit_stored = decimal(reference["gap_phys_minus_oph"])
    inside = gap_lower <= deficit_recomputed <= gap_upper
    recorded_flag = bridge["verdict"]["reference_deficit_inside_certified_gap"]

    payload_integral = payload["integral"]
    measure_moment = measure["transport_moments"]["timelike_on_shell_mz"]
    payload_covariance = payload.get("covariance", {})
    source_compilation = payload.get("source_compilation", {})
    payload_kernel = payload.get("kernel", {})
    tabulated_protocol = protocol["classes"]["tabulated_dispersive"]

    checks = {
        "endpoint_recomputed_exactly": gate(
            endpoint_central == stored_endpoint_central
            and [endpoint_lower, endpoint_upper] == stored_endpoint_interval,
            "central value and interval recomputed from anchor, lepton packet, "
            "and tabulated HVP value",
        ),
        "same_scheme_gap_recomputed_exactly": gate(
            [gap_lower, gap_upper] == stored_gap,
            "gap interval recomputed from the recorded CODATA comparison, "
            "anchor, lepton packet, and HVP uncertainty",
        ),
        "protocol_gap_matches_recomputed_gap": gate(
            protocol_gap == [gap_lower, gap_upper],
            "packet accounting interval matches the byte-pinned source row",
        ),
        "reference_decomposition_recomputed": gate(
            abs(reference_alpha_recomputed - reference_alpha_stored)
            <= REFERENCE_TOLERANCE,
            "on-shell reference recomputed from its four primitive decimal "
            "terms; tolerance covers the producer's serialized binary float",
        ),
        "reference_deficit_recomputed": gate(
            abs(deficit_recomputed - deficit_stored) <= REFERENCE_TOLERANCE,
            "reference deficit recomputed rather than trusted",
        ),
        "reference_flag_recomputed": gate(
            isinstance(recorded_flag, bool) and recorded_flag is inside,
            "stored containment flag agrees with recomputation",
        ),
        "anchor_matches_bridge": gate(
            abs(anchor - bridge_anchor) <= REFERENCE_TOLERANCE,
            "endpoint and bridge use the same serialized anchor",
        ),
        "payload_release_matches": gate(
            endpoint_inputs["payload_release"]
            == payload["data_release"]["release_id"]
            == tabulated_protocol["release_id"],
            "endpoint, payload, and protocol identify one KNT19 release",
        ),
        "payload_value_matches": gate(
            delta_had
            == decimal(payload_integral["value"])
            == decimal(tabulated_protocol["delta_alpha_had_5_MZ"]),
            "tabulated HVP central value agrees across all pinned surfaces",
        ),
        "payload_uncertainty_matches": gate(
            delta_had_unc
            == decimal(payload_integral["uncertainty"])
            == decimal(tabulated_protocol["uncertainty"]),
            "tabulated HVP uncertainty agrees across all pinned surfaces",
        ),
        "payload_units_explicit": gate(
            payload_integral.get("unit") == "dimensionless"
            and tabulated_protocol["unit"] == "dimensionless",
            "the HVP moment is dimensionless",
        ),
        "payload_kernel_explicit": gate(
            payload_kernel.get("name") == "subtracted_vacuum_polarization_dispersion"
            and isinstance(payload_kernel.get("formula"), str)
            and bool(payload_kernel["formula"].strip())
            and payload_kernel.get("target") == "Delta_alpha_had_5_MZ",
            "kernel, convention, and target quantity are serialized",
        ),
        "payload_covariance_scope_explicit": gate(
            set(payload_covariance) == {"statistical", "systematic"}
            and isinstance(payload_covariance["statistical"].get("policy"), str)
            and isinstance(payload_covariance["systematic"].get("policy"), str),
            "the tabulated scalar replay records its uncertainty policy; no "
            "raw channel covariance is claimed",
        ),
        "payload_provenance_explicit": gate(
            source_compilation.get("id") == "knt19_pinned_piecewise_v1"
            and source_compilation.get("url") == "https://arxiv.org/abs/1911.00367"
            and "Phys. Rev. D 101, 014029" in source_compilation.get("citation", ""),
            "published compilation, URL, citation, and equation are identified",
        ),
        "exclusions_scope_explicit": gate(
            tabulated_protocol["exclusions_status"]
            == (
                "NOT_APPLICABLE_TO_PUBLISHED_SCALAR_REPLAY__"
                "RAW_CHANNEL_EXCLUSIONS_NOT_INGESTED"
            ),
            "this is a scalar-table replay, not a raw-channel dispersive ingest",
        ),
        "spectral_export_release_matches": gate(
            measure["provenance"]["data_release"]["release_id"]
            == payload["data_release"]["release_id"],
            "spectral export and tabulated payload identify one release",
        ),
        "spectral_export_moment_matches": gate(
            decimal(measure_moment["value"]) == decimal(payload_integral["value"]),
            "spectral export reproduces the pinned tabulated moment",
        ),
        "spectral_export_requadrature_gate": gate(
            measure["consistency"]["within_tolerance"] is True
            and decimal(measure["consistency"]["abs_difference"])
            <= decimal(measure["consistency"]["tolerance"]),
            "the producer's declared requadrature check is attained",
        ),
        "target_and_promotion_guards": gate(
            endpoint["guards"]["promotable_as_oph_source_theorem"] is False
            and endpoint["guards"]["measured_alpha_in_solve_path"] is False
            and bridge["guards"]["public_promotion_allowed"] is False
            and bridge["guards"]["measured_values_in_any_oph_solve_path"] is False
            and protocol["scope"]["physical_alpha_prediction"] is False,
            "serialized producer guards prohibit source or prediction promotion",
        ),
    }
    internally_consistent = all(row["passed"] for row in checks.values())

    if not internally_consistent:
        row_verdict = "INTERNAL_INCONSISTENCY"
    elif inside:
        row_verdict = "COMPATIBLE_RECORDED_ACCOUNTING_REPLAY"
    else:
        row_verdict = "REFUTED_RECORDED_ACCOUNTING_REPLAY"

    return {
        "class": "tabulated_dispersive",
        "evaluation_mode": "recorded_accounting_replay",
        "independent_hvp_implementation": False,
        "prospective_comparison": False,
        "endpoint_alpha_inv_central": str(endpoint_central),
        "endpoint_interval": [str(endpoint_lower), str(endpoint_upper)],
        "same_scheme_gap_interval": [str(gap_lower), str(gap_upper)],
        "reference_deficit": str(deficit_recomputed),
        "reference_deficit_inside_gap": bool(inside),
        "gates": checks,
        "all_internal_gates_pass": internally_consistent,
        "class_verdict": row_verdict,
        "claim_boundary": (
            "A secondary arithmetic implementation reproduces one recorded "
            "KNT19-based accounting row. It is not an independent HVP "
            "implementation, not a blind comparison, and not a physical "
            "alpha prediction."
        ),
    }


def _overall_verdict(tabulated_verdict: str) -> str:
    if tabulated_verdict == "COMPATIBLE_RECORDED_ACCOUNTING_REPLAY":
        return "MULTI_CLASS_NOT_EVALUABLE__ONE_RECORDED_ACCOUNTING_REPLAY_COMPATIBLE"
    if tabulated_verdict == "REFUTED_RECORDED_ACCOUNTING_REPLAY":
        return "MULTI_CLASS_NOT_EVALUABLE__RECORDED_ACCOUNTING_REPLAY_REFUTED"
    return "PROTOCOL_INCONSISTENT"


def build_verdict(
    protocol_path: Path = PROTOCOL_PATH, root: Path = REPO_ROOT
) -> dict[str, Any]:
    protocol = load_protocol(protocol_path)
    verify_snapshot_commit(protocol, root)
    pinned = verify_exact_pins(protocol, root)
    endpoint = _load_pinned_json(pinned, ENDPOINT_REL)
    bridge = _load_pinned_json(pinned, BRIDGE_REL)
    payload = _load_pinned_json(pinned, PAYLOAD_REL)
    measure = _load_pinned_json(pinned, MEASURE_REL)

    tabulated = evaluate_tabulated_accounting_replay(
        endpoint, bridge, payload, measure, protocol
    )
    class_matrix: dict[str, Any] = {"tabulated_dispersive": tabulated}
    for name in UNAVAILABLE_CLASSES:
        row = protocol["classes"][name]
        class_matrix[name] = {
            "class": name,
            "evaluation_mode": row["evaluation_mode"],
            "independent_hvp_implementation": False,
            "class_verdict": "NOT_EVALUABLE_MISSING_FROZEN_INGEST",
            "ingest_requirement": row["ingest_requirement"],
            "fabrication_excluded": True,
        }

    replay_verdicts = {
        "COMPATIBLE_RECORDED_ACCOUNTING_REPLAY",
        "REFUTED_RECORDED_ACCOUNTING_REPLAY",
    }
    accounting_replay_count = sum(
        row["class_verdict"] in replay_verdicts for row in class_matrix.values()
    )
    independent_evaluated_count = sum(
        row.get("independent_hvp_implementation") is True
        and not row["class_verdict"].startswith("NOT_EVALUABLE")
        for row in class_matrix.values()
    )

    payload_out = {
        "schema": SCHEMA,
        "issue": ISSUE,
        "row_class": "retrospective_empirical_same_scheme_accounting_audit",
        "source_snapshot": {
            "commit": protocol["preexisting_source_snapshot"]["commit"],
            "sha256": protocol["preexisting_source_snapshot"]["sha256"],
            "commit_ancestor_verified": True,
            "commit_blob_pins_verified": True,
            "all_exact_pins_verified": True,
        },
        "scope": {
            "comparison_timing": "retrospective",
            "prospective_freeze": False,
            "independent_hvp_implementation_supplied": False,
            "empirical_input_promoted_to_source_output": False,
            "physical_alpha_prediction_emitted": False,
        },
        "serialized_acceptance_gates": {
            "exact_source_pins": gate(
                True,
                "the snapshot commit is ancestral and all declared commit "
                "blobs and working-tree source bytes match its pins",
            ),
            "recorded_accounting_replay": gate(
                tabulated["all_internal_gates_pass"],
                "one tabulated accounting row is arithmetically replayed",
            ),
            "independent_multi_class_evaluation": {
                "status": "NOT_ATTAINED",
                "passed": False,
                "detail": (
                    "zero independent HVP implementations are evaluable; "
                    "three required frozen ingests are absent"
                ),
            },
            "prediction_or_freeze_promotion": {
                "status": "PROHIBITED",
                "passed": True,
                "detail": (
                    "comparison data predate the packet, so no prospective "
                    "freeze or prediction status is available"
                ),
            },
        },
        "class_matrix": class_matrix,
        "cross_class_agreement": {
            "recorded_accounting_replay_count": accounting_replay_count,
            "independently_evaluated_class_count": independent_evaluated_count,
            "verdict": "NOT_EVALUABLE_NO_INDEPENDENT_CLASS",
            "note": (
                "cross-class agreement requires at least two genuinely "
                "independent evaluated classes under one protocol"
            ),
        },
        "verdict": _overall_verdict(tabulated["class_verdict"]),
        "reactivation": {
            name: protocol["classes"][name]["ingest_requirement"]
            for name in UNAVAILABLE_CLASSES
        }
        | {
            "rule": (
                "A new frozen ingest creates a new versioned evaluation. "
                "This retrospective packet and its source pins do not move."
            )
        },
        "claim_boundary": (
            "The multi-class independent alpha/HVP test is not evaluable. "
            "One byte-pinned KNT19 accounting row is compatible under a "
            "secondary arithmetic replay. Raw-dispersive, independent-code, "
            "and lattice-HVP classes lack frozen repository ingests. The "
            "result is retrospective and supplies neither a prospective "
            "freeze nor a physical OPH alpha prediction."
        ),
    }
    payload_out["verdict_sha256"] = sha256_bytes(
        canonical_json(
            {
                key: value
                for key, value in payload_out.items()
                if key != "verdict_sha256"
            }
        ).encode("utf-8")
    )
    return payload_out


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    parser.add_argument(
        "--source-root",
        type=Path,
        default=REPO_ROOT,
        help="repository checkout whose pinned source bytes are verified",
    )
    args = parser.parse_args(argv)
    try:
        payload = build_verdict(root=args.source_root)
    except (OSError, KeyError, TypeError, VerdictInputError) as error:
        print(f"VERDICT_INPUT_ERROR: {error}", file=sys.stderr)
        return 1
    if args.verify:
        stored = json.loads(OUT_PATH.read_text(encoding="utf-8"))
        if stored != payload:
            print("VERDICT_DRIFT", file=sys.stderr)
            return 1
        print("VERDICT_VERIFIED")
        return 0
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(
        json.dumps(payload, sort_keys=True, indent=1) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"verdict": payload["verdict"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
