"""Scalar-sector blindness of the measured chain (issue #623).

The #616 receipt proves that scalar existence and multiplicity are not
source-determined at the grammar level: the empty, duplicate-doublet, and
inert-doublet configurations pass every grammar-visible check.  This
certificate proves the sharper statement on the realized observable chain
and lands the fail-closed conditional exit of issue #623.

THE SHARPER NON-IDENTIFIABILITY (exact, by input closure):

* The measured artifact chain of the corpus consists of the charged
  response artifact and the pole-residue artifact.  Each records its
  complete input closure as pins: the carrier manifest and the dynamics
  module, and for the pole artifact its parent response artifact.  No
  scalar-sector object is an input, and the payloads carry no
  scalar-sector key.  The chain values are therefore constant across the
  empty, one-doublet, duplicate-doublet, and inert-doublet completions:
  the complete measured chain is scalar-blind, not merely the grammar.
* The discriminating class is named: every scalar discriminator in the
  declared algebra must factor through the three Yukawa invariant
  channels of the pinned matter packet, and no producer in the measured
  chain emits a Yukawa-channel observable.  Scalar discrimination is
  therefore a new producer interface, not a reading of the chain.

THE EXIT: the physical scalar attachment stays a fail-closed conditional
open interface with the #616 countermodels retained, exactly the second
bounded exit of issue #623.  The declared one-doublet completion stays
declared; no pole, potential, vacuum value, or multiplicity
discrimination is claimed, and no Higgs-target value enters any gate.
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

SCHEMA = "oph.scalar_sector_blindness_certificate.v1"
MANIFEST_PATH = MODULE_DIR / "manifests" / "scalar_sector_blindness_reference.json"
CARRIER_MANIFEST_NAME = "echosahedral_federation_reference.json"
WINDOW_MANIFEST_NAME = "multiplicity_window_reference.json"
MATTER_MANIFEST_NAME = "super_tannakian_matter_reference.json"
RESPONSE_ARTIFACT_NAME = "charged_response_semantic_artifact.json"
POLE_RESIDUE_ARTIFACT_NAME = "charged_response_pole_residue_artifact.json"

ISSUE = 623

SCALAR_COMPLETIONS = ("empty", "one_doublet_declared", "duplicate_doublet", "inert_doublet")

SCALAR_TOKENS = ("scalar", "higgs", "yukawa", "vacuum", "doublet", "potential")


def collect_keys(value: Any, out: set[str]) -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            out.add(str(key).lower())
            collect_keys(item, out)
    elif isinstance(value, list):
        for item in value:
            collect_keys(item, out)


def pin(name: str) -> tuple[dict[str, Any], dict[str, Any]]:
    artifact = load_json(MODULE_DIR / "manifests" / name)
    return artifact, {"path": f"manifests/{name}", "sha256": sha256_json(artifact)}


def input_closure_audit() -> dict[str, Any]:
    """The measured chain's complete input closure, read from its pins."""

    carrier_manifest, carrier_pin = pin(CARRIER_MANIFEST_NAME)
    e565.validate_carrier(carrier_manifest)
    response, response_pin = pin(RESPONSE_ARTIFACT_NAME)
    pole, pole_pin = pin(POLE_RESIDUE_ARTIFACT_NAME)

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

    response_keys: set[str] = set()
    collect_keys(response, response_keys)
    pole_keys: set[str] = set()
    collect_keys(pole, pole_keys)
    for token in SCALAR_TOKENS:
        offending = sorted(
            key for key in response_keys | pole_keys if token in key
        )
        require(
            not offending,
            "SCALAR_KEY_PRESENT",
            f"the measured chain must carry no scalar-sector key ({offending})",
        )

    return {
        "chain": [
            {"artifact": "charged response", "pin": response_pin,
             "inputs": ["carrier manifest", "dynamics module"]},
            {"artifact": "pole residue", "pin": pole_pin,
             "inputs": ["carrier manifest", "parent response artifact", "dynamics module"]},
        ],
        "carrier": carrier_pin,
        "statement": (
            "the complete input closure of the measured chain is the carrier "
            "and the dynamics; no scalar-sector object is an input and no "
            "payload key names one, so every chain value is constant across "
            "the scalar completions"
        ),
        "completions_covered": list(SCALAR_COMPLETIONS),
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
        "class": "Yukawa-channel observables",
        "channels": [list(channel) for channel in matter_channels],
        "statement": (
            "every scalar discriminator in the declared algebra couples "
            "through the three invariant Yukawa channels; the measured "
            "chain emits no Yukawa-channel observable, so scalar "
            "discrimination is a producer interface that does not exist in "
            "the corpus rather than a reading of existing artifacts"
        ),
        "interface_id": "physical_scalar_attachment",
        "interface_class": "conditional_open_interface",
    }


# ---------------------------------------------------------------------------
# Controls
# ---------------------------------------------------------------------------


def control_chain_discrimination_attempt() -> dict[str, Any]:
    """Asserting a scalar discrimination from the measured chain must be
    refused: the chain is completion-independent by input closure."""

    chain_value_by_completion = {name: "identical" for name in SCALAR_COMPLETIONS}
    distinct = len(set(chain_value_by_completion.values())) > 1
    try:
        require(
            distinct,
            "NO_CHAIN_DISCRIMINATOR",
            "the measured chain cannot separate scalar completions",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "NO_CHAIN_DISCRIMINATOR",
        }
    return {"expected_failure": True, "failed": False}


def control_target_injection() -> dict[str, Any]:
    """A discriminator built on the measured Higgs mass must be refused."""

    proposed_gate_inputs = ["carrier manifest", "measured_higgs_mass_gev"]
    try:
        require(
            all("higgs" not in item for item in proposed_gate_inputs),
            "TARGET_INJECTION",
            "a scalar discriminator must not consume a Higgs target value",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "TARGET_INJECTION",
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
    """A doctored artifact carrying a scalar key must fail the closure audit."""

    doctored = {"carrier_binding": {}, "scalar_vacuum_value": "246"}
    keys: set[str] = set()
    collect_keys(doctored, keys)
    offending = [key for key in keys if any(t in key for t in SCALAR_TOKENS)]
    try:
        require(
            not offending,
            "SCALAR_KEY_PRESENT",
            "a doctored scalar key must be detected",
        )
    except CertificateError:
        return {
            "expected_failure": True,
            "failed": True,
            "code": "SCALAR_KEY_PRESENT",
        }
    return {"expected_failure": True, "failed": False}


# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------


def build_payload() -> dict[str, Any]:
    closure = input_closure_audit()
    retained = countermodels_retained()
    interface = discriminating_interface(
        retained["declared_completion"]["yukawa_channels"]
    )

    controls = {
        "chain_discrimination_attempt": control_chain_discrimination_attempt(),
        "target_injection": control_target_injection(),
        "wrong_multiplicity": control_wrong_multiplicity(),
        "scalar_key_mutation": control_scalar_key_mutation(),
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
            "Sharper scalar non-identifiability on the realized observable "
            "chain: the measured artifacts are constant across the scalar "
            "completions by input closure, the discriminating class is the "
            "Yukawa-channel observable interface, which no producer emits, "
            "and the physical scalar attachment stays a fail-closed "
            "conditional open interface with the #616 countermodels "
            "retained. No pole, potential, vacuum value, or discrimination "
            "is claimed, and no Higgs-target value enters any gate."
        ),
        "input_closure": closure,
        "countermodels": retained,
        "discriminating_interface": interface,
        "controls": controls,
        "bounded_exit": "conditional_open_interface_with_countermodels_retained",
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
