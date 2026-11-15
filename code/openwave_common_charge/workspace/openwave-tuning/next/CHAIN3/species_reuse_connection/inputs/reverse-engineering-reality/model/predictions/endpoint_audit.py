"""Retrospective arithmetic audit of the inherited P-response endpoint.

The interval certifies a declared incomplete numerical map. Comparing that
interval with CODATA does not promote it to a physical transport theorem.
All differences below use exact rationals; displayed divisions round outward.
"""
from __future__ import annotations

from copy import deepcopy
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, localcontext
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from model.observer.p_response import experiment as response


ROOT = Path(__file__).resolve().parents[2]
MODE = "thomson_structured_running"
RECEIPT = "model/observer/p_response/results.json"
CONSUMER = "model/observer/p_response/experiment.py"
SOURCE_PINS = {
    **response.PINS,
    CONSUMER: "44049e5702af83ab8a62e40f8f7d3b25b088ed24e47078a0b412c3b2400ed61e",
    RECEIPT: "e050ec01837c59459b0ea1c1f839f8993d5930b97276b972ee82c3366f9ad91f",
}

# Published, previously examined data: never call this a prospective freeze.
_TARGET = {
    "schema": "oph.predictions.empirical_endpoint.v1",
    "observable": "inverse_fine_structure_constant_thomson_limit",
    "unit": "dimensionless",
    "source": "CODATA 2022 recommended values, NIST wall table",
    "citation": "https://physics.nist.gov/cuu/pdf/wall_2022.pdf",
    "value": "137.035999177",
    "standard_uncertainty": "0.000000021",
    "comparison_half_width_in_standard_uncertainties": "5",
    "exposure": "retrospective_known_data",
}


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


TARGET_SHA256 = hashlib.sha256(_canonical(_TARGET)).hexdigest()


def target_fixture():
    """Independent copy of the cited target; edits require a new audit version."""
    return deepcopy(_TARGET)


def _positive_decimal(value, name):
    if not isinstance(value, str):
        raise ValueError(f"{name} must be an exact decimal string")
    try:
        parsed = Decimal(value)
    except ArithmeticError as error:
        raise ValueError(f"Invalid {name}") from error
    if not parsed.is_finite() or parsed <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return Fraction(parsed)


def _target(value):
    if not isinstance(value, dict) or value.keys() != _TARGET.keys():
        raise ValueError("Target must contain the pinned CODATA fixture fields")
    centre = _positive_decimal(value["value"], "target value")
    uncertainty = _positive_decimal(value["standard_uncertainty"], "standard uncertainty")
    multiple = _positive_decimal(value["comparison_half_width_in_standard_uncertainties"], "comparison multiplier")
    try:
        digest = hashlib.sha256(_canonical(value)).hexdigest()
    except (TypeError, ValueError) as error:
        raise ValueError("Invalid target fixture") from error
    if digest != TARGET_SHA256:
        raise ValueError("Target fixture differs from the pinned cited comparison")
    return centre, uncertainty, multiple


def _read_sources(root):
    for path, expected in SOURCE_PINS.items():
        if hashlib.sha256((root / path).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Source pin mismatch: {path}")
    certificate = json.loads((root / response.CERT).read_text())
    receipt = json.loads((root / RECEIPT).read_text())
    if (response.MODE != MODE or receipt["source_branch"]["mode"] != MODE
            or receipt["source_pins"] != response.PINS):
        raise ValueError("Inherited detector and audit do not use the same pinned branch")
    if (certificate["promotion_allowed"] is not False
            or certificate["exact_alpha_promoted"] is not False
            or certificate["consumer_policy"]["may_feed_compare_or_audit_surfaces"] is not True
            or certificate["consumer_policy"]["may_feed_live_particle_predictions"] is not False):
        raise ValueError("Inherited incomplete-map claim boundary changed")
    mode = certificate["modes"][MODE]
    enclosure = mode["certified_enclosure"]
    if receipt["source_branch"]["P_imported_interval"] != enclosure["P"]:
        raise ValueError("Detector P input differs from the inherited closure enclosure")
    lo = _positive_decimal(enclosure["alpha_inv"]["lo"], "source inverse-alpha lower endpoint")
    hi = _positive_decimal(enclosure["alpha_inv"]["hi"], "source inverse-alpha upper endpoint")
    if lo > hi:
        raise ValueError("Reversed source interval")
    return certificate, mode, lo, hi


def _decimal(value, rounding):
    with localcontext() as context:
        context.prec = 80
        context.rounding = rounding
        return format(Decimal(value.numerator) / Decimal(value.denominator), "f")


def _interval(lo, hi):
    if lo > hi:
        raise ValueError("Reversed output interval")
    return {"lo": _decimal(lo, ROUND_FLOOR), "hi": _decimal(hi, ROUND_CEILING)}


def audit(target=None, *, root=None):
    """Compare inherited fixed output with known data without fitting or promotion.

    ``target`` accepts an independently read copy of the exact cited fixture,
    not arbitrary replacement measurements. Different data need a separately
    identified comparison. ``root`` permits source-integrity verification in
    an isolated checkout; no files are written and no source value is changed.
    """
    fixture = target_fixture() if target is None else deepcopy(target)
    centre, uncertainty, multiple = _target(fixture)
    certificate, mode, source_lo, source_hi = _read_sources(ROOT if root is None else Path(root))
    allowance = multiple * uncertainty
    band_lo, band_hi = centre - allowance, centre + allowance
    deficit_lo, deficit_hi = centre - source_hi, centre - source_lo
    correction_lo, correction_hi = band_lo - source_hi, band_hi - source_lo
    disjoint = source_hi < band_lo or band_hi < source_lo
    minimum_correction = max(Fraction(0), correction_lo, -correction_hi)
    return {
        "schema": "oph.predictions.endpoint_audit.v1",
        "comparison_class": "retrospective_known_data",
        "target": fixture,
        "target_sha256": TARGET_SHA256,
        "source": {
            "mode": MODE,
            "pins": dict(SOURCE_PINS),
            "certificate": response.CERT,
            "certificate_replayed_here": False,
            "artifact_hashes_verified": True,
            "inverse_alpha_interval": _interval(source_lo, source_hi),
            "forward_P_interval": deepcopy(mode["certified_enclosure"]["P"]),
            "inherited_claim_status": certificate["claim_status"],
            "inherited_claim_boundary": certificate["claim_boundary"],
            "scope": "Declared structured-running numerical map; its narrow interval is not physical theory uncertainty",
        },
        "arithmetic": {
            "target_comparison_band": _interval(band_lo, band_hi),
            "target_centre_minus_source_interval": _interval(deficit_lo, deficit_hi),
            "relative_inverse_alpha_deficit_ppm": _interval(
                1_000_000 * deficit_lo / centre, 1_000_000 * deficit_hi / centre),
            "correction_needed_to_enter_comparison_band": _interval(correction_lo, correction_hi),
            "minimum_absolute_correction_to_enter_band": _decimal(minimum_correction, ROUND_FLOOR),
            "source_interval_disjoint_from_comparison_band": disjoint,
            "method": "Exact rational arithmetic from decimal input endpoints; displayed quotients rounded outward at 80 digits",
        },
        "verdicts": {
            "inherited_incomplete_transport_contract": "physical_completion_pending",
            "literal_complete_zero_correction_hypothesis": (
                "inconsistent_with_declared_comparison_band" if disjoint else "compatible_with_declared_comparison_band"),
            "OPH_as_a_whole": "not_determined",
        },
        "interpretation": {
            "comparison_band": "Central value plus or minus five reported standard uncertainties; a declared comparison allowance, not a rigorous confidence interval",
            "missing_correction": "Endpoint residual at the fixed inherited forward P; not a derived transport term or a self-consistent solution of a completed closure map",
            "theory_uncertainty": "Incomplete physical transport has no supplied remainder budget; no sigma significance is assigned",
            "counterfactual_scope": "The inconsistent hypothesis additionally declares the incomplete source interval to be the entire physical endpoint with zero missing correction",
            "fit_applied": False,
            "P_recalibrated_from_target": False,
            "prospective_prediction": False,
            "empirical_validation_of_OPH": False,
        },
    }
