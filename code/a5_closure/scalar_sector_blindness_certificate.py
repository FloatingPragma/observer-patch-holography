"""Limited scalar-blindness audit of two simulator artifacts (issue #623).

The #616 receipt proves that scalar existence and multiplicity are not
source-determined by its enumerated grammar-visible checks: the empty,
duplicate-doublet, and inert-doublet rows remain exact countermodels at that
scope.  This certificate audits a narrower object: the selected charged
response and pole-residue simulator artifacts.

WHAT IS EXACT:

* Both selected artifacts have valid schemas and self-hashes.  Their carrier
  and parent pins agree, and their serialized interfaces contain no explicit
  scalar-sector input.
* This proves only that the selected two-artifact screen-response subchain
  exposes no scalar-completion coordinate.  It is not an exhaustive producer
  inventory, a transitive source-dependency proof, or a theorem on a typed
  category of physical scalar completions.
* Scalar-sensitive observables need not factor through Yukawa channels.
  Relevant missing interfaces include scalar poles and residues,
  multiplet-resolved spectral rank, gauge-covariant response, potential and
  vacuum observables, and Yukawa channels.  The selected pair supplies no
  source-bound physical scalar attachment; other corpus producers are outside
  this audit.

The physical scalar attachment remains a named conditional interface.  The
declared one-doublet completion stays declared; no pole, potential, vacuum
value, or multiplicity discrimination is claimed, and no Higgs-target value
enters any gate.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Mapping, Sequence

MODULE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(MODULE_DIR))

import echosahedral_selector_certificate as e565  # noqa: E402

CertificateError = e565.CertificateError
require = e565.require
sha256_json = e565.sha256_json
load_json = e565.load_json
write_json = e565.write_json

SCHEMA = "oph.scalar_sector_blindness_certificate.v2"
MANIFEST_PATH = MODULE_DIR / "manifests" / "scalar_sector_blindness_reference.json"
CARRIER_MANIFEST_NAME = "echosahedral_federation_reference.json"
WINDOW_MANIFEST_NAME = "multiplicity_window_reference.json"
MATTER_MANIFEST_NAME = "super_tannakian_matter_reference.json"
RESPONSE_ARTIFACT_NAME = "charged_response_semantic_artifact.json"
POLE_RESIDUE_ARTIFACT_NAME = "charged_response_pole_residue_artifact.json"

ISSUE = 623

SCALAR_COMPLETIONS = ("empty", "one_doublet_declared", "duplicate_doublet", "inert_doublet")

SCALAR_TOKENS = ("scalar", "higgs", "yukawa", "vacuum", "doublet", "potential")


def artifact_self_hash(artifact: Mapping[str, Any]) -> str:
    body = {
        key: value for key, value in artifact.items() if key != "artifact_sha256"
    }
    return "sha256:" + sha256_json(body)


def pin(name: str) -> tuple[dict[str, Any], dict[str, Any]]:
    artifact = load_json(MODULE_DIR / "manifests" / name)
    return artifact, {"path": f"manifests/{name}", "sha256": sha256_json(artifact)}


def validate_selected_subchain(
    carrier_manifest: Mapping[str, Any],
    response: Mapping[str, Any],
    pole: Mapping[str, Any],
) -> None:
    """Validate the two selected artifacts without promoting their scope."""

    e565.validate_carrier(carrier_manifest)
    require(
        response.get("schema") == "oph.charged_response_semantic_artifact.v3"
        and response.get("issue") == 599,
        "RESPONSE_SCHEMA",
        "the selected response artifact must carry its exact schema and issue",
    )
    require(
        pole.get("schema") == "oph.charged_response_pole_residue.v2"
        and pole.get("issue") == 569,
        "POLE_SCHEMA",
        "the selected pole artifact must carry its exact schema and issue",
    )
    require(
        response.get("artifact_sha256") == artifact_self_hash(response),
        "RESPONSE_SELF_HASH",
        "the selected response artifact self-hash does not match",
    )
    require(
        pole.get("artifact_sha256") == artifact_self_hash(pole),
        "POLE_SELF_HASH",
        "the selected pole artifact self-hash does not match",
    )
    carrier_sha = sha256_json(carrier_manifest)
    require(
        response["carrier_binding"]["carrier_manifest_sha256"] == carrier_sha,
        "INPUT_CLOSURE",
        "the response artifact must pin the carrier manifest as its input",
    )
    require(
        pole["carrier_binding"]["carrier_manifest_sha256"] == carrier_sha
        and pole["carrier_binding"]["parent_artifact_sha256"]
        == response["artifact_sha256"],
        "INPUT_CLOSURE",
        "the pole artifact must pin the carrier and its parent artifact",
    )

    serialized_interface = json.dumps(
        {"response": response, "pole": pole},
        sort_keys=True,
        separators=(",", ":"),
    ).lower()
    for token in SCALAR_TOKENS:
        require(
            token not in serialized_interface,
            "EXPLICIT_SCALAR_INPUT",
            f"the selected artifact interfaces contain scalar token {token!r}",
        )


def selected_subchain_audit() -> dict[str, Any]:
    """Audit the selected pair and state the unproved closure boundary."""

    carrier_manifest, carrier_pin = pin(CARRIER_MANIFEST_NAME)
    response, response_pin = pin(RESPONSE_ARTIFACT_NAME)
    pole, pole_pin = pin(POLE_RESIDUE_ARTIFACT_NAME)
    validate_selected_subchain(carrier_manifest, response, pole)

    return {
        "chain": [
            {"artifact": "charged response", "pin": response_pin,
             "declared_inputs": ["carrier manifest", "runtime dynamics report"]},
            {"artifact": "pole residue", "pin": pole_pin,
             "declared_inputs": [
                 "carrier manifest",
                 "parent response artifact",
                 "dynamics module byte digest",
             ]},
        ],
        "carrier": carrier_pin,
        "statement": (
            "the selected two-artifact finite simulator screen-response "
            "subchain has valid self-hashes and carrier/parent pins and "
            "exposes no explicit scalar-completion input"
        ),
        "scope_limits": [
            "not an exhaustive corpus producer inventory",
            "not a transitive source or import closure proof",
            "not a completion-indexed observation functor",
            "not a physical scalar-sector non-identifiability theorem",
        ],
        "dynamics_hash_semantics": {
            "charged_response": "runtime-report payload digest",
            "pole_residue": "dynamics-module byte digest",
            "interchangeable": False,
        },
        "grammar_countermodel_labels_not_physical_completions": list(
            SCALAR_COMPLETIONS
        ),
    }


def countermodels_retained() -> dict[str, Any]:
    """Re-pin the #616 non-determination receipt and the declared scalar."""

    window, window_pin = pin(WINDOW_MANIFEST_NAME)
    scalar_block = window["scalar_response_multiplicity"]
    matter, matter_pin = pin(MATTER_MANIFEST_NAME)
    contract = matter["exterior_matter_contract"]
    require(
        contract["one_scalar"] == "weak_block" and len(contract["yukawa_channels"]) == 3,
        "DECLARED_SCALAR",
        "the declared completion must carry one weak-block scalar with three channels",
    )
    return {
        "grammar_receipt": window_pin,
        "grammar_statement": (
            "scalar existence and multiplicity are not source-determined: "
            "the empty, duplicate-doublet, and inert-doublet configurations "
            "pass every grammar-visible check (pinned #616 receipt)"
        ),
        "scalar_block_keys": sorted(scalar_block) if isinstance(scalar_block, Mapping) else "pinned",
        "declared_completion": {
            "pin": matter_pin,
            "scalar": "one color-singlet weak-doublet scalar on the weak block",
            "yukawa_channels": contract["yukawa_channels"],
        },
    }


def discriminating_interface(matter_channels: Sequence[Sequence[str]]) -> dict[str, Any]:
    return {
        "status": "non_exhaustive_missing_producer_classes",
        "classes": [
            "scalar two-point pole and residue observables",
            "multiplet-resolved scalar spectral rank",
            "gauge-covariant kinetic and current response",
            "scalar potential and vacuum observables",
            "Yukawa-channel observables",
        ],
        "declared_yukawa_channels": [list(channel) for channel in matter_channels],
        "statement": (
            "the selected screen-response pair emits none of these "
            "scalar-sensitive observables. Yukawa channels are one possible "
            "class, not an exhaustive discriminator. The selected pair "
            "supplies no source-bound physical scalar attachment; the status "
            "of other corpus producers is outside this audit"
        ),
        "interface_id": "physical_scalar_attachment",
        "interface_class": "conditional_open_interface",
    }


# ---------------------------------------------------------------------------
# Controls
# ---------------------------------------------------------------------------


def control_scope_promotion() -> dict[str, Any]:
    """The limited interface audit must not be promoted to a category theorem."""

    proved_scope = "selected_two_artifact_subchain"
    requested_scope = "all_scalar_completions_in_typed_source_category"
    try:
        require(
            proved_scope == requested_scope,
            "SCOPE_PROMOTION",
            "the selected artifact audit is not a typed-category non-identifiability theorem",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "SCOPE_PROMOTION",
        }
    return {"expected_failure": True, "failed": False}


def mutate_response_input(key: str, value: Any) -> tuple[dict[str, Any], dict[str, Any]]:
    """Return a self-consistently rehashed response and dependent pole."""

    response = load_json(MODULE_DIR / "manifests" / RESPONSE_ARTIFACT_NAME)
    pole = load_json(MODULE_DIR / "manifests" / POLE_RESIDUE_ARTIFACT_NAME)
    response["explicit_inputs"] = {key: value}
    response["artifact_sha256"] = artifact_self_hash(response)
    pole["carrier_binding"]["parent_artifact_sha256"] = response["artifact_sha256"]
    pole["artifact_sha256"] = artifact_self_hash(pole)
    return response, pole


def control_target_injection() -> dict[str, Any]:
    """A rehashed artifact importing the measured Higgs mass must be refused."""

    carrier = load_json(MODULE_DIR / "manifests" / CARRIER_MANIFEST_NAME)
    response, pole = mutate_response_input("measured_higgs_mass_gev", 125.25)
    try:
        validate_selected_subchain(carrier, response, pole)
    except CertificateError as error:
        return {
            "expected_failure": True,
            "failed": True,
            "code": error.code,
        }
    return {"expected_failure": True, "failed": False}


def control_wrong_multiplicity() -> dict[str, Any]:
    """A claimed scalar-multiplicity selection must be refused: the pinned
    #616 scalar verdict carries the three countermodels and keeps the
    declared completion declared."""

    window = load_json(MODULE_DIR / "manifests" / WINDOW_MANIFEST_NAME)
    verdict = window["scalar_response_multiplicity"]["verdict"]
    countermodels = set(verdict["countermodels"])
    expected = {"n0_no_scalar", "n2_duplicate_identical_charge", "n2_one_inert"}
    try:
        require(
            countermodels != expected,
            "WRONG_MULTIPLICITY",
            "the scalar countermodel battery stands, so a multiplicity selection claim is refused",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "WRONG_MULTIPLICITY",
            "countermodels": sorted(countermodels),
        }
    return {"expected_failure": True, "failed": False}


def control_scalar_key_mutation() -> dict[str, Any]:
    """A self-consistently rehashed explicit scalar input must still fail."""

    carrier = load_json(MODULE_DIR / "manifests" / CARRIER_MANIFEST_NAME)
    response, pole = mutate_response_input(
        "scalar_completion", "inert_doublet"
    )
    try:
        validate_selected_subchain(carrier, response, pole)
    except CertificateError as error:
        return {
            "expected_failure": True,
            "failed": True,
            "code": error.code,
        }
    return {"expected_failure": True, "failed": False}


def control_parent_pin_mutation() -> dict[str, Any]:
    carrier = load_json(MODULE_DIR / "manifests" / CARRIER_MANIFEST_NAME)
    response = load_json(MODULE_DIR / "manifests" / RESPONSE_ARTIFACT_NAME)
    pole = load_json(MODULE_DIR / "manifests" / POLE_RESIDUE_ARTIFACT_NAME)
    doctored = json.loads(json.dumps(pole))
    doctored["carrier_binding"]["parent_artifact_sha256"] = "sha256:" + "0" * 64
    doctored["artifact_sha256"] = artifact_self_hash(doctored)
    try:
        validate_selected_subchain(carrier, response, doctored)
    except CertificateError as error:
        return {"expected_failure": True, "failed": True, "code": error.code}
    return {"expected_failure": True, "failed": False}


def control_self_hash_mutation() -> dict[str, Any]:
    carrier = load_json(MODULE_DIR / "manifests" / CARRIER_MANIFEST_NAME)
    response = load_json(MODULE_DIR / "manifests" / RESPONSE_ARTIFACT_NAME)
    pole = load_json(MODULE_DIR / "manifests" / POLE_RESIDUE_ARTIFACT_NAME)
    doctored = json.loads(json.dumps(response))
    doctored["orientation_convention"] = "doctored"
    try:
        validate_selected_subchain(carrier, doctored, pole)
    except CertificateError as error:
        return {"expected_failure": True, "failed": True, "code": error.code}
    return {"expected_failure": True, "failed": False}


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def build_payload() -> dict[str, Any]:
    closure = selected_subchain_audit()
    retained = countermodels_retained()
    interface = discriminating_interface(
        retained["declared_completion"]["yukawa_channels"]
    )

    controls = {
        "scope_promotion": control_scope_promotion(),
        "target_injection": control_target_injection(),
        "wrong_multiplicity": control_wrong_multiplicity(),
        "scalar_key_mutation": control_scalar_key_mutation(),
        "parent_pin_mutation": control_parent_pin_mutation(),
        "self_hash_mutation": control_self_hash_mutation(),
    }
    for name, verdict in controls.items():
        require(
            verdict["expected_failure"] is True and verdict["failed"] is True,
            "CONTROL_NOT_FAILED",
            f"control {name} did not record its required failure",
        )

    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "issue": ISSUE,
        "claim_boundary": (
            "Limited two-artifact finite simulator screen-response audit: "
            "the selected charged-response and pole-residue artifacts expose "
            "no explicit scalar-completion input. This is not an exhaustive "
            "producer inventory, transitive source closure, or typed-category "
            "non-identifiability theorem. Scalar poles, gauge response, "
            "potential/vacuum observables, and Yukawa channels remain possible "
            "discriminators. The physical scalar attachment remains open with "
            "the #616 grammar-scope countermodels retained."
        ),
        "input_closure": closure,
        "countermodels": retained,
        "discriminating_interface": interface,
        "controls": controls,
        "bounded_exit": "limited_subchain_audit_physical_interface_open",
    }
    return payload


def build_manifest() -> dict[str, Any]:
    payload = build_payload()
    manifest = dict(payload)
    manifest["manifest_sha256"] = "sha256:" + sha256_json(payload)
    return manifest


def verify_stored() -> dict[str, Any]:
    stored = load_json(MANIFEST_PATH)
    body = {key: value for key, value in stored.items() if key != "manifest_sha256"}
    require(
        stored.get("manifest_sha256") == "sha256:" + sha256_json(body),
        "MANIFEST_HASH",
        "stored manifest hash does not match its body",
    )
    require(
        body == build_payload(),
        "MANIFEST_DRIFT",
        "stored manifest does not match a deterministic rebuild",
    )
    return {"status": "PASS", "manifest": str(MANIFEST_PATH)}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Scalar-sector blindness certificate for issue #623")
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--output", type=Path, default=MANIFEST_PATH)
    args = parser.parse_args(argv)
    if args.verify:
        print(json.dumps(verify_stored(), indent=2))
        return 0
    manifest = build_manifest()
    write_json(args.output, manifest)
    print(json.dumps({"status": "WROTE", "manifest": str(args.output), "manifest_sha256": manifest["manifest_sha256"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
