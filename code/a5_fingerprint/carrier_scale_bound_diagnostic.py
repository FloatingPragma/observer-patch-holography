#!/usr/bin/env python3
"""Carrier-scale bound diagnostic from published vacuum-dispersion limits.

The diagnostic translates one published Lorentz-invariance-violation
bound into an upper bound on the carrier scale, on the declared
equal-weight stencil branch with a declared photon-sector assignment.
It is typed exposed-retrospective: the consumed bound is public, the
output constrains an open parameter, and nothing here scores, confirms,
or falsifies the framework.

Chain, exact where the sources are exact:

* the certified stencil dispersion is
  ``omega^2 = k^2 [1 - (a^2/20) k^2 + O(a^4 k^4)]``, subluminal, with
  even powers only, so the class produces no linear (n = 1) vacuum
  dispersion at any order and the strongest published linear bounds
  constrain it not at all;
* the standard quadratic parameterization is
  ``omega^2 = k^2 [1 - (k/E_QG2)^2]``, so the identification is exact:
  ``a = sqrt(20) / E_QG2``;
* the consumed public bound is E_QG2 > 1.3 x 10^11 GeV at 95 percent
  confidence for subluminal quadratic dispersion from GRB 090510
  (Vasileiou et al., Physical Review D 87, 122001 (2013)); the
  conversion constant is hbar c = 1.97327 x 10^-16 GeV m (CODATA,
  declared exposed input).

Reading: the carrier grain on this branch is finer than about
7 x 10^-27 m, roughly 4 x 10^8 Planck lengths above the Planck length,
so a Planck-scale carrier passes with eight orders of headroom and the
diagnostic neither confirms nor endangers the framework. The
anisotropic spin-six residue at the same scale enters at (a k)^4 and is
far below current sensitivity at accessible energies; the isotropic
tower carries the binding constraint, and the spin-six shape stays the
discriminating fingerprint for any future resolved anisotropy.
Off-stencil members of the invariant class change the order-one
isotropic coefficient, so the bound is exact on the declared branch and
order-of-magnitude on the class. Ladder row FZ-11 and the issue #639
comparison budget are untouched.
"""

from __future__ import annotations

import argparse
import json
from fractions import Fraction
from pathlib import Path
from typing import Any

import a5_multipole_fixed_point_certificate as base

HERE = Path(__file__).resolve().parent
RUNTIME = HERE / "runtime"
RECEIPT_PATH = RUNTIME / "carrier_scale_bound_receipt.json"

SCHEMA = "oph.carrier_scale_bound_receipt.v1"
STATUS = "EXPOSED_RETROSPECTIVE_CARRIER_SCALE_BOUND__DIAGNOSTIC_ONLY"

# declared exposed public inputs
E_QG2_GEV = Fraction(13, 10) * 10**11          # Vasileiou et al. 2013, 95% CL
HBARC_GEV_M = Fraction(197327, 10**21)          # 1.97327e-16 GeV m, CODATA
PLANCK_LENGTH_M = Fraction(1616255, 10**41)     # 1.616255e-35 m, CODATA


def require(condition: bool, message: str) -> None:
    if not condition:
        raise base.FingerprintError(message)


def sci(x: Fraction, digits: int = 4) -> str:
    """Deterministic scientific-notation rendering of a positive Fraction."""

    require(x > 0, "scientific rendering needs a positive value")
    exponent = 0
    y = x
    while y >= 10:
        y /= 10
        exponent += 1
    while y < 1:
        y *= 10
        exponent -= 1
    scaled = y * 10 ** (digits - 1)
    mantissa_int = scaled.numerator // scaled.denominator
    mantissa = str(mantissa_int)
    return f"{mantissa[0]}.{mantissa[1:]}e{exponent:+03d}"


def build_receipt() -> dict[str, Any]:
    # pin the parent stencil coefficient from the fingerprint receipt
    parent_bytes = base.RECEIPT_PATH.read_bytes()
    parent = json.loads(parent_bytes)
    expansion = parent["kinetic_stencil_conditional"]["expansion"]
    require(
        "- (a^2/20) k^4" in expansion,
        "stencil isotropic coefficient drift against the pinned receipt",
    )

    # exact identification a = sqrt(20) / E_QG2, reported through a^2
    a_squared_gev = Fraction(20, 1) / (E_QG2_GEV * E_QG2_GEV)
    a_squared_m = a_squared_gev * HBARC_GEV_M * HBARC_GEV_M
    # a = sqrt(a_squared): render via integer square root of a scaled value
    scale = 10**80
    a_m_scaled = _isqrt((a_squared_m * scale * scale).numerator
                        // (a_squared_m * scale * scale).denominator)
    a_m = Fraction(a_m_scaled, scale)
    require(
        a_m * a_m <= a_squared_m < (a_m + Fraction(1, scale)) ** 2,
        "square-root enclosure drift",
    )
    headroom = a_m / PLANCK_LENGTH_M

    return {
        "schema": SCHEMA,
        "status": STATUS,
        "issue": 655,
        "parent_pins": [
            {
                "path": (
                    "code/a5_fingerprint/runtime/"
                    "a5_multipole_fixed_point_receipt.json"
                ),
                "bytes": len(parent_bytes),
                "sha256": base.tagged_sha256(parent_bytes),
            }
        ],
        "declared_premises": {
            "branch": (
                "equal-weight stencil member of the invariant class; "
                "off-stencil members change the order-one isotropic "
                "coefficient, so the class-level bound is "
                "order-of-magnitude"
            ),
            "sector": (
                "photon-sector assignment of the carrier dispersion; the "
                "assignment belongs to the open realization map"
            ),
        },
        "exposed_public_inputs": {
            "E_QG2_GeV": sci(E_QG2_GEV),
            "source": (
                "Vasileiou et al., Physical Review D 87, 122001 (2013): "
                "E_QG2 > 1.3e11 GeV at 95 percent confidence, subluminal "
                "quadratic vacuum dispersion, GRB 090510"
            ),
            "hbar_c_GeV_m": sci(HBARC_GEV_M),
            "planck_length_m": sci(PLANCK_LENGTH_M),
            "exposure_class": "EXPOSED_RETROSPECTIVE",
        },
        "translation": {
            "identification": (
                "omega^2 = k^2[1 - (a^2/20)k^2] matches "
                "omega^2 = k^2[1 - (k/E_QG2)^2] exactly, so "
                "a = sqrt(20)/E_QG2"
            ),
            "no_linear_term": (
                "the hop expansion carries even powers only, so the class "
                "produces no n = 1 vacuum dispersion at any order and the "
                "stronger published linear bounds impose no constraint"
            ),
            "subluminal_branch_matches": True,
        },
        "bound": {
            "carrier_scale_upper_bound_m": sci(a_m),
            "planck_length_headroom": sci(headroom),
            "reading": (
                "the carrier grain on the declared branch is finer than "
                "about 7e-27 m; a Planck-scale carrier sits eight orders "
                "of magnitude below the bound, so current public data "
                "neither confirm nor endanger the framework and every "
                "future tightening lowers the ceiling"
            ),
            "spin_six_note": (
                "at the bound-saturating scale the anisotropic spin-six "
                "residue enters at (a k)^4 and is far below current "
                "sensitivity at accessible energies; the isotropic tower "
                "carries the binding constraint and the rigid spin-six "
                "template stays the discriminating fingerprint for any "
                "future resolved anisotropy"
            ),
        },
        "comparison_boundary": {
            "public_measurement_read": True,
            "reading_type": "published upper-limit consumption only",
            "comparison_permitted": False,
            "scored": False,
            "fz11_untouched": True,
            "issue_639_budget_untouched": True,
            "typing": (
                "diagnostic of an open parameter, in the same class as the "
                "charged-lepton and screen-grain diagnostics; never a "
                "prediction, never evidence"
            ),
        },
    } | {"receipt_sha256": ""}


def _isqrt(n: int) -> int:
    import math

    return math.isqrt(n)


def finalize(receipt: dict[str, Any]) -> dict[str, Any]:
    receipt = dict(receipt)
    receipt.pop("receipt_sha256", None)
    receipt["receipt_sha256"] = base.tagged_sha256(
        base.canonical_json_bytes(receipt)
    )
    return receipt


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args(argv)
    receipt = finalize(build_receipt())
    if args.write:
        RUNTIME.mkdir(exist_ok=True)
        RECEIPT_PATH.write_bytes(base.canonical_json_bytes(receipt))
        print(RECEIPT_PATH)
        return 0
    print(json.dumps(receipt["bound"], indent=1))
    print(receipt["status"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
