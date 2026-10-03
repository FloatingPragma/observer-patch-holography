#!/usr/bin/env python3
"""Retrospective comparison of both fixed conditional mass branches.

The producer never imports this module or either target fixture. Quoted PDG
ranges are retained as ranges, without interpreting 90% limits as one sigma.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from decimal import Decimal, localcontext
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PACKET = ROOT / "code/particles/runs/flavor/conditional_quark_mass_replay.json"
PDG = ROOT / "code/particles/data/pdg_2025_conditional_quark_comparison.json"
FLAG = ROOT / "code/particles/data/flag_2024_light_quark_ratio_fixture.json"
SPEC = Path(__file__).with_name("conditional_quark_mass_spec.json")
VERIFIER = Path(__file__).with_name("verify_conditional_quark_mass_replay.py")
TOP_PACKET = ROOT / "code/particles/runs/calibration/conditional_top_fixed_p.json"
TOP_VERIFIER = ROOT / "code/particles/calibration/verify_conditional_top_fixed_p.py"
OUTPUT = ROOT / "code/particles/runs/flavor/conditional_quark_mass_comparison.json"
BRANCHES = ("full_rscc", "lower_order_control")
D = Decimal


def quoted_range(value, row):
    central, minus, plus = (D(row[k]) for k in ("central", "minus", "plus"))
    low, high = central - minus, central + plus
    return {
        "computed": str(value), "reference": row,
        "quoted_lower": str(low), "quoted_upper": str(high),
        "inside_quoted_marginal_range": low <= value <= high,
        "signed_difference": str(value - central),
    }


def compare(packet, pdg, flag):
    spec = json.loads(SPEC.read_text())
    if packet["guards"] != spec["guards"]:
        raise ValueError("All frozen scientific guards must be retained")
    if packet["model"]["charts"] != spec["charts"]:
        raise ValueError("The declared native mass charts must be retained")
    expected_interpretation = {
        "joint_covariance_available": False, "pdg_and_flag_independent": False,
        "theory_uncertainty_supplied": False,
        "no_gaussian_sigma_conversion_of_quoted_ranges": True,
        "mass_chart_attachment_is_conditional": True,
        "six_marginal_matches_are_not_joint_model_acceptance": True,
        "all_scheduled_light_comparisons_must_be_retained": True,
        "target_informed_formula_discovery_remains_target_informed": True,
        "direct_reconstruction_top_mass_not_used_as_pole_mass": True,
        "new_independent_evidence_claimed": False,
    }
    if pdg["interpretation"] != expected_interpretation:
        raise ValueError("The comparison interpretation cannot be promoted")
    if [r["quark"] for r in pdg["masses"]] != list("udscbt"):
        raise ValueError("All six PDG comparisons required in fixed order")
    if [r["quantity"] for r in pdg["light_quantities"]] != ["mu_over_md", "ms_over_mud", "mud_gev"]:
        raise ValueError("All three light comparisons required")
    if [r["nf"] for r in flag["averages"]] != ["2+1+1", "2+1"]:
        raise ValueError("Both FLAG flavor theories must be retained")
    top_reference = pdg["criticality_top_comparison"]
    if (top_reference["scheme"], top_reference["scale"], top_reference["unit"]) != (
        "MSbar", "m_t(m_t), cross-section extraction", "GeV"
    ):
        raise ValueError("Criticality comparison must retain its running top chart")
    for row in pdg["masses"]:
        q = row["quark"]
        scheme = "pole_extraction" if q == "t" else "MSbar"
        scale = "2 GeV" if q in "uds" else f"m_{q}(m_{q})"
        if q == "t":
            scale = "pole from cross-section measurements"
        if row["scheme"] != scheme or row["scale"] != scale or row["unit"] != "GeV":
            raise ValueError("Comparison chart changed")
    with localcontext() as ctx:
        ctx.prec = 32
        branches = {}
        for name in BRANCHES:
            masses = {q: D(v) for q, v in packet[name]["native_mass_chart_GeV"].items()}
            if set(masses) != set("udscbt") or any(not m.is_finite() or m <= 0 for m in masses.values()):
                raise ValueError("Six finite positive mass coordinates required")
            values = {
                "mu_over_md": masses["u"] / masses["d"],
                "ms_over_mud": 2 * masses["s"] / (masses["u"] + masses["d"]),
                "mud_gev": (masses["u"] + masses["d"]) / 2,
            }
            flag_rows = []
            for average in flag["averages"]:
                for quantity in ("mu_over_md", "ms_over_mud"):
                    reference = average[quantity]
                    value = values[quantity]
                    flag_rows.append({
                        "nf": average["nf"], "quantity": quantity,
                        "computed": str(value), "reference": reference,
                        "marginal_standardized_residual": str(
                            (value - D(reference["value"])) / D(reference["standard_uncertainty"])),
                        "includes_theory_uncertainty": False,
                        "joint_or_independent_significance": False,
                    })
            branches[name] = {
                "pdg_individual_masses": [quoted_range(masses[r["quark"]], r) for r in pdg["masses"]],
                "pdg_light_quantities": [quoted_range(values[r["quantity"]], r) for r in pdg["light_quantities"]],
                "flag_light_ratios": flag_rows,
            }
        return {
            "schema": "oph.conditional_quark_mass_comparison.v1",
            "status": "RETROSPECTIVE_CONDITIONAL_COMPARISON",
            "primary_branch": "full_rscc",
            "control_role": "fixed_historical_lower_order_ablation_not_selected_by_agreement",
            "branches": branches,
            "interpretation": pdg["interpretation"],
            "joint_likelihood": None,
            "theory_uncertainty": None,
            "sources": {"pdg": pdg["source"], "flag": flag["source"]},
        }


def build():
    packet, pdg, flag = (json.loads(p.read_text()) for p in (PACKET, PDG, FLAG))
    module_spec = importlib.util.spec_from_file_location("independent_mass_verifier", VERIFIER)
    verifier = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(verifier)
    verifier.verify_payload(packet)
    result = compare(packet, pdg, flag)
    top = json.loads(TOP_PACKET.read_text())
    top_spec = importlib.util.spec_from_file_location("independent_top_verifier", TOP_VERIFIER)
    top_verifier = importlib.util.module_from_spec(top_spec)
    top_spec.loader.exec_module(top_verifier)
    top_verifier.verify_payload(top)
    if top["source_scales_GeV"]["E_star"] != float(json.loads(SPEC.read_text())["calibrations"]["E_star_GeV"]["value"]):
        raise ValueError("Both calculations must use the same external unit calibration")
    top_rows = {
        "one_loop": top["hybrid_rows"]["1"][-1],
        "two_loop_yukawa_quartic_one_loop_gauges": top["hybrid_rows"]["2"][-1],
        "coupled_two_loop_top_only": top["fully_coupled_two_loop_same_low_gauge_anchor"][-1],
    }
    result["criticality_top_running_coordinate_comparisons"] = {
        name: quoted_range(D(str(row["top_running_mass_coordinate_GeV"])), pdg["criticality_top_comparison"])
        for name, row in top_rows.items()
    }
    result["criticality_top_boundary"] = {
        "physical_electroweak_vev_matching_assumed": True,
        "finite_electroweak_and_complex_pole_matching_supplied": False,
        "QCD_pole_proxy_promoted_to_physical_pole": False,
        "truncation_spread_is_a_confidence_interval": False,
        "branches_selected_by_comparison": False,
    }
    result["source_pins"] = {
        p.relative_to(ROOT).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (PACKET, PDG, FLAG, SPEC, VERIFIER, TOP_PACKET, TOP_VERIFIER, Path(__file__))
    }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    rendered = json.dumps(build(), indent=2, sort_keys=True) + "\n"
    if args.check:
        if OUTPUT.read_text() != rendered:
            raise SystemExit("Conditional comparison receipt differs; regenerate it")
        print("Conditional comparison receipt verified")
    else:
        OUTPUT.write_text(rendered)
        print(OUTPUT)


if __name__ == "__main__":
    main()
