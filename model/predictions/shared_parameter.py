"""One fixed forward P in the coupling endpoint and conditional CMB tilt.

This downstream comparison propagates the inherited numerical enclosure;
it supplies neither the incomplete endpoint transport nor a primordial map.
Predictions are sealed before the retained ACT projection is read. No fitting,
likelihood multiplication, or prospective-prediction promotion is performed.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
from pathlib import Path

from . import endpoint_audit as endpoint
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins

ROOT = Path(__file__).resolve().parents[2]
PROJECTION = "code/cosmology/tilt_posterior/act_dr6_ns_weight.json.gz"
PROJECTION_SHA256 = "13dab540b652300bd000cf1ac3b197af55c7c411de215682482ea55993f40485"
ARCHIVE_SHA256 = "8f2801a69663ef1142bb05660e0118ef7b1e0bdf3d05afba434a0569023b98ab"
INHERITED_P = F("1.6309682094039593")
OUTPUT = ROOT / "model/predictions/shared_parameter_results.json"


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def digest(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def tilt_interval(p_lower, p_upper):
    """Exact image of the decreasing affine map; numerical errors only."""
    if not isinstance(p_lower, F) or not isinstance(p_upper, F):
        raise ValueError("exact rational P endpoints required")
    if not 0 < p_lower <= p_upper < 48:
        raise ValueError("ordered positive P interval below 48 required")
    return 1 - p_upper / 48, 1 - p_lower / 48


def frozen_forward_prediction(*, root=ROOT):
    # This stage reads source certificates only, never cosmology or CODATA data.
    certificate, mode, alpha_lo, alpha_hi = endpoint._read_sources(Path(root))
    box = mode["certified_enclosure"]["P"]
    p_lo, p_hi = F(box["lo"]), F(box["hi"])
    ns_lo, ns_hi = tilt_interval(p_lo, p_hi)
    return {
        "mode": endpoint.MODE,
        "P_interval": [str(p_lo), str(p_hi)],
        "inverse_alpha_interval": [str(alpha_lo), str(alpha_hi)],
        "scalar_tilt_interval": [str(ns_lo), str(ns_hi)],
        "tilt_map": "n_s=1-P/48",
        "coordinate_identification": "conditional identification of the primordial P with the forward coupling P",
        "inherited_forward_map_scope": certificate["claim_boundary"],
        "numerical_enclosure_is_physical_theory_error": False,
    }


def weighted_positions(rows, lower, upper, inherited):
    """Compare original decimal tokens exactly, including endpoint equality."""
    if not lower <= upper:
        raise ValueError("reversed scalar-tilt interval")
    total = lower_count = upper_count = inherited_count = 0
    for row in rows:
        if not isinstance(row, list) or len(row) != 2 or any(type(v) is not str for v in row):
            raise ValueError("original weight and n_s string tokens required")
        weight, ns = map(F, row)
        if weight <= 0 or weight.denominator != 1 or not 0 < ns < 2:
            raise ValueError("positive integer multiplicity and finite scalar tilt required")
        total += int(weight)
        lower_count += int(weight) * (ns <= lower)
        upper_count += int(weight) * (ns <= upper)
        inherited_count += int(weight) * (ns <= inherited)
    if not total:
        raise ValueError("empty retained chain")
    return {"multiplicity": total, "forward_tail_weight_interval": [lower_count, upper_count],
            "inherited_tail_weight": inherited_count}


def compare_projection(prediction, *, root=ROOT):
    raw = (Path(root) / PROJECTION).read_bytes()
    if hashlib.sha256(raw).hexdigest() != PROJECTION_SHA256:
        raise ValueError("ACT projection hash mismatch")
    payload = json.loads(gzip.decompress(raw))
    if (payload["schema"] != "oph.act_dr6.ns_projection.v1"
            or payload["archive_sha256"] != ARCHIVE_SHA256
            or str(payload["burn_in_fraction"]) != "0.5"
            or len(payload["chains"]) != 4
            or len({c["member"] for c in payload["chains"]}) != 4):
        raise ValueError("official four-chain burn-in contract changed")
    lower, upper = map(F, prediction["scalar_tilt_interval"])
    inherited = 1 - INHERITED_P / 48
    per_chain = []
    for chain in payload["chains"]:
        drop = len(chain["rows"]) // 2
        retained = chain["rows"][drop:]
        per_chain.append({"member": chain["member"], "source_sha256": bound_pins(chain["source_sha256"]),
                          "burn_in_rows": drop, "retained_rows": len(retained),
                          **weighted_positions(retained, lower, upper, inherited)})
    total = sum(c["multiplicity"] for c in per_chain)
    tail = [sum(c["forward_tail_weight_interval"][i] for c in per_chain) for i in (0, 1)]
    old_tail = sum(c["inherited_tail_weight"] for c in per_chain)
    return {
        "projection": PROJECTION, "projection_sha256": PROJECTION_SHA256,
        "archive_sha256": ARCHIVE_SHA256, "per_chain": per_chain,
        "retained_rows": sum(c["retained_rows"] for c in per_chain),
        "multiplicity": total,
        "forward_tail_weight_interval": tail,
        "forward_posterior_mass_interval": [str(F(value, total)) for value in tail],
        "inherited_P": str(INHERITED_P), "inherited_scalar_tilt": str(inherited),
        "inherited_tail_weight": old_tail,
        "forward_minus_inherited_tilt_interval": [str(lower-inherited), str(upper-inherited)],
        "forward_and_inherited_have_same_tail_weight": tail == [old_tail, old_tail],
        "posterior_mass_is_p_value": False,
        "tail_effective_sample_size_or_confidence_certificate": False,
    }


def build(*, root=ROOT):
    prediction = frozen_forward_prediction(root=root)
    seal = digest(prediction)
    cosmology = compare_projection(prediction, root=root)
    coupling = endpoint.audit(root=root)
    if digest(prediction) != seal:
        raise ValueError("forward prediction changed after comparison")
    return {
        "schema": "oph.predictions.shared_forward_parameter.v1",
        "classification": "retrospective conditional same-parameter comparison",
        "prediction": prediction, "prediction_sha256_before_comparison": seal,
        "source_pins": bound_pins(dict(endpoint.SOURCE_PINS)),
        "coupling": {"target": coupling["target"],
                     "relative_inverse_alpha_deficit_ppm": coupling["arithmetic"]["relative_inverse_alpha_deficit_ppm"],
                     "zero_missing_correction_verdict": coupling["verdicts"]["literal_complete_zero_correction_hypothesis"]},
        "cosmology": cosmology,
        "coverage": {
            "one_forward_P_used_in_both_maps": True,
            "target_fit_or_recalibration": False,
            "full_common_physical_construction_established": False,
            "unresolved_physical_maps": ["complete low-energy coupling transport", "native curvature perturbation and primordial tilt map", "CMB background and radiation transfer"],
            "observations_used": ["CODATA 2022 inverse fine-structure constant", "ACT DR6 plus Planck, lensing and DESI DR1 conditional scalar-tilt posterior"],
            "untested_by_this_comparison": ["particle spectra and holonomy defects", "physical observer clocks", "quantum intervention records", "gravitational wave and local gravity data"],
            "joint_likelihood": None, "FZ15_verdict": None,
            "all_public_datasets_reproduced": False,
            "within_stated_uncertainties_across_scales_established": False,
            "interpretation": "The zero-correction coupling endpoint fails its declared band; the CMB result is a conditional posterior position. No joint OPH-wide acceptance or exclusion follows.",
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = build()
    if args.check:
        if json.loads(OUTPUT.read_text(encoding="utf-8")) != result:
            raise SystemExit("shared-parameter comparison drift")
    else:
        with OUTPUT.open("x", encoding="utf-8") as stream:
            stream.write(canonical(result))
    print(canonical({"coupling": result["coupling"]["zero_missing_correction_verdict"],
                     "CMB_tail_weight_interval": result["cosmology"]["forward_tail_weight_interval"],
                     "CMB_total_multiplicity": result["cosmology"]["multiplicity"],
                     "same_P_in_both_maps": True}), end="")


if __name__ == "__main__":
    main()
