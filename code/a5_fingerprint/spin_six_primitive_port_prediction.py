#!/usr/bin/env python3
"""Freeze the primitive twelve-port spin-six propagation prediction.

This producer reads no comparison data.  It combines the exact twelve-port
moment certificate with the invariant-class theorem and freezes the physical
prediction made by one named OPH branch:

* the intrinsic scalar propagation symbol uses the complete primitive
  twelve-port orbit and no additional directional hops at the displayed
  orders;
* proper-carrier covariance acts transitively on those ports, hence every
  invariant coefficient is equal;
* the symbol is normalized to its continuum k^2 term and has one finite,
  nonzero carrier scale a;
* for a photon reading, the symbol acts equally on both transverse
  polarizations and one carrier frame is transported coherently into the
  comparison frame.

Those premises select

    omega^2(k,n) = (1/(2 a^2)) sum_i [1 - cos(a k u_i.n)].

The expansion predicts three linked coefficients.  Once C4 is measured, no
carrier scale or amplitude remains adjustable:

    C4 = -a^2/20,
    B0 = a^4/840,
    B6 = 2 a^4/7875,
    B6/C4^2 = 32/315,
    B0/C4^2 = 10/21,
    B6/B0 = 16/75.

The unknown carrier orientation is fitted in SO(3)/A5.  This is a prospective
conditional prediction of the primitive-port physical branch.  Issue #655
owns the stronger task of deriving the physical sector bridge from the repair
law.  A null remains inconclusive without a source lower bound on a.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import a5_multipole_fixed_point_certificate as base
import spin_six_universality_certificate as universality

HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
RECEIPT_PATH = RUNTIME / "spin_six_primitive_port_prediction_receipt.json"

SCHEMA = "oph.spin_six_primitive_port_prediction.v1"
STATUS = (
    "FROZEN_PROSPECTIVE_PRIMITIVE_TWELVE_PORT_BRANCH_PREDICTION__"
    "PHYSICAL_COMPARISON_UNARMED"
)
FROZEN_UTC = "2026-07-31T18:13:34Z"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise base.FingerprintError(message)


def load_typed_receipt(
    path: Path, *, schema: str, status: str
) -> tuple[bytes, dict[str, Any]]:
    """Fail closed on bytes, schema, status, and the receipt self-digest."""

    raw = path.read_bytes()
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError as error:
        raise base.FingerprintError(f"invalid parent JSON: {path.name}") from error
    require(raw == base.canonical_json_bytes(payload), f"noncanonical {path.name}")
    require(payload.get("schema") == schema, f"schema drift in {path.name}")
    require(payload.get("status") == status, f"status drift in {path.name}")
    body = {key: value for key, value in payload.items() if key != "receipt_sha256"}
    require(
        payload.get("receipt_sha256")
        == base.tagged_sha256(base.canonical_json_bytes(body)),
        f"self-digest drift in {path.name}",
    )
    return raw, payload


def build_receipt() -> dict[str, Any]:
    base_bytes, base_receipt = load_typed_receipt(
        base.RECEIPT_PATH, schema=base.SCHEMA, status=base.STATUS
    )
    universal_bytes, universal_receipt = load_typed_receipt(
        universality.RECEIPT_PATH,
        schema=universality.SCHEMA,
        status=universality.STATUS,
    )

    expansion = base_receipt["kinetic_stencil_conditional"]["expansion"]
    require("- (a^2/20) k^4" in expansion, "C4 parent coefficient drift")
    require("(a^4/840) k^6" in expansion, "B0 parent coefficient drift")
    require("(2 a^4/7875) k^6 I6" in expansion, "B6 parent coefficient drift")
    require(
        universal_receipt["all_order_level_six_coefficient"][
            "strictly_positive"
        ]
        is True,
        "all-order level-six certificate drift",
    )

    c4 = Fraction(-1, 20)
    b0 = Fraction(1, 840)
    b6 = Fraction(2, 7875)
    ratio_b6_c4 = b6 / (c4 * c4)
    ratio_b0_c4 = b0 / (c4 * c4)
    ratio_b6_b0 = b6 / b0
    require(ratio_b6_c4 == Fraction(32, 315), "B6/C4^2 drift")
    require(ratio_b0_c4 == Fraction(10, 21), "B0/C4^2 drift")
    require(ratio_b6_b0 == Fraction(16, 75), "B6/B0 drift")

    receipt: dict[str, Any] = {
        "schema": SCHEMA,
        "status": STATUS,
        "issue": 655,
        "frozen_utc": FROZEN_UTC,
        "parent_pins": [
            {
                "path": (
                    "code/a5_fingerprint/runtime/"
                    "a5_multipole_fixed_point_receipt.json"
                ),
                "bytes": len(base_bytes),
                "sha256": base.tagged_sha256(base_bytes),
            },
            {
                "path": (
                    "code/a5_fingerprint/runtime/"
                    "spin_six_universality_receipt.json"
                ),
                "bytes": len(universal_bytes),
                "sha256": base.tagged_sha256(universal_bytes),
            },
        ],
        "prediction_scope": {
            "name": "primitive twelve-port scalar propagation branch",
            "type": "prospective conditional physical-branch prediction",
            "stronger_derivation_open": (
                "issue #655 must derive or reject the physical sector bridge; "
                "a derivation broadens the scope without changing this target"
            ),
            "spin_language": (
                "spin six denotes spherical-harmonic rank j=6, not particle spin"
            ),
        },
        "branch_premises": {
            "support": (
                "the complete primitive twelve-port orbit is the sole intrinsic "
                "directional support through the displayed derivative order"
            ),
            "equal_weights": (
                "proper-carrier covariance and transitivity force one equal "
                "coefficient on the primitive orbit"
            ),
            "normalization": (
                "the continuum quadratic term is k^2, fixing the common "
                "coefficient to 1/(2 a^2)"
            ),
            "finite_scale": "a is finite and strictly positive",
            "physical_sector": (
                "the tested sector realizes this scalar symbol; a photon test "
                "also requires equal action on both transverse polarizations"
            ),
            "frame": (
                "one carrier rest frame is transported coherently over the "
                "experiment; its constant orientation in SO(3)/A5 is profiled"
            ),
            "exclusivity": (
                "the fitted intrinsic coefficients isolate the carrier term "
                "from source, medium, gravitational, and instrumental effects"
            ),
        },
        "exact_prediction": {
            "operator": (
                "omega^2(k,n) = (1/(2 a^2)) sum_{i=1}^{12} "
                "[1 - cos(a k u_i.n)]"
            ),
            "template": (
                "I6(n) = (25/132) sum_i P6(u_i.n), normalized to one on "
                "the twelve vertex directions"
            ),
            "expansion": (
                "omega^2 = k^2 - (a^2/20) k^4 + (a^4/840) k^6 + "
                "(2 a^4/7875) k^6 I6(n) + O(a^6 k^8)"
            ),
            "coefficients": {
                "C4_over_a2": str(c4),
                "B0_over_a4": str(b0),
                "B6_over_a4": str(b6),
            },
            "scale_free_relations": {
                "B6_over_C4_squared": str(ratio_b6_c4),
                "B0_over_C4_squared": str(ratio_b0_c4),
                "B6_over_B0": str(ratio_b6_b0),
            },
            "signs": {"C4": "negative", "B0": "positive", "B6": "positive"},
            "harmonic_nulls": "all intrinsic anisotropic j=1,2,3,4,5 coefficients vanish",
            "fit_freedom_after_C4": (
                "one orientation in SO(3)/A5; no scale, amplitude, or "
                "additional j=6 shape coefficient remains"
            ),
        },
        "baseline_contrast": {
            "baseline": "minimal locally Lorentz-invariant Standard Model plus General Relativity in local vacuum",
            "baseline_prediction": "C4 = B0 = B6 = 0 for intrinsic vacuum propagation",
            "nonuniqueness": (
                "nonminimal effective operators or another icosahedral medium "
                "can imitate the coefficient pattern; a match distinguishes "
                "the branch from the baseline without identifying OPH uniquely"
            ),
        },
        "prospective_decision_rule": {
            "eligible_data": (
                "a post-freeze data release with a joint likelihood or full "
                "covariance for same-sector C4, isotropic B0, and the j=6 "
                "coefficient vector, plus declared source, medium, gravity, "
                "instrument, boost, and orientation treatment"
            ),
            "trigger": (
                "C4 is measured negative at at least five standard deviations "
                "with enough sensitivity to resolve both linked k^6 terms"
            ),
            "fail": (
                "after profiling the predeclared SO(3)/A5 orientation, the "
                "linked B0 and B6 prediction is excluded at at least five "
                "standard deviations with calibrated joint coverage"
            ),
            "support": (
                "the zero-coefficient baseline is excluded at at least five "
                "standard deviations, the linked OPH manifold agrees within "
                "two standard deviations, named systematics are rejected, and "
                "an independent release replicates the result"
            ),
            "inconclusive": (
                "C4 is not detected, the linked terms are below sensitivity, "
                "the covariance is incomplete, or the carrier contribution "
                "cannot be isolated"
            ),
            "scope_of_failure": (
                "failure rejects the primitive twelve-port physical propagation "
                "branch; it rejects OPH as a whole only if issue #655 proves "
                "that branch is forced and exclusive"
            ),
        },
        "exposure_and_custody_boundary": {
            "new_comparison_data_read": False,
            "comparison_permitted": False,
            "comparison_state": "UNARMED_PENDING_DATASET_SPECIFIC_PREREGISTRATION",
            "prior_related_exposure": (
                "a 2026-07-17 WMAP ILC search tested the A5 l=6 and higher "
                "sky templates and returned a null family-wide p=0.64"
            ),
            "prior_preregistration_sha256": (
                "2b83a001f75d4aa9f5a631b50d0fe8ad51950ae63146118ae369e8cfb80e84b2"
            ),
            "excluded_data_class": (
                "the WMAP search, its CMB likelihood class, and every data "
                "product inspected in that campaign are ineligible for a new "
                "FZ-11 verdict"
            ),
            "uncompared_novel_content": (
                "the linked C4, B0, and B6 coefficient relations have not been "
                "compared with a qualifying physical dataset"
            ),
        },
    }
    receipt["receipt_sha256"] = base.tagged_sha256(base.canonical_json_bytes(receipt))
    return receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    receipt = build_receipt()
    if args.write:
        RUNTIME.mkdir(exist_ok=True)
        RECEIPT_PATH.write_bytes(base.canonical_json_bytes(receipt))
        print(RECEIPT_PATH)
        return 0
    print(json.dumps(receipt["exact_prediction"], indent=1))
    print(receipt["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
