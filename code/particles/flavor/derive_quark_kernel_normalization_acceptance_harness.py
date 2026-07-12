#!/usr/bin/env python3
"""Executable acceptance harness for the issue-377 kernel program (K1-K3).

The closure route posted on #377 requires an OPH-derived transport kernel
whose normalized spectrum emits the two quark spread moduli. This module
turns the acceptance tests of that program into a mechanical gate: feed it
a candidate-kernel payload and it either blocks with the first failed gate
or reports acceptance (which unlocks promotion review; it is not itself a
promotion).

Candidate payload schema (JSON):

    {
      "candidate_id": "<name>",
      "kernel": {"refinements": [{"level": 0, "hermitian_descendant":
                 {"real": [[...]], "imag": [[...]]}}, ...]},
      "normalization_clause": {"kind": "<budget|modular|...>",
                               "sigma_u_per_side": <float>,
                               "sigma_d_per_side": <float>,
                               "statement": "<one line>"},
      "g_ch_replacement": {"value": <float>, "statement": "..."},
      "ancestry": {"artifacts": ["<artifact names>"],
                   "attestations": {
                     "quark_reference_values_consumed": false,
                     "fitted_spreads_consumed": false,
                     "numerical_flavor_template_consumed": false}}
    }

Gates:

  G0 schema: required fields present.
  G1 ancestry: every attestation false, including template-consumption,
     and no forbidden ancestor (the forbidden list is imported from the
     nonidentifiability obstruction module so the two artifacts cannot
     drift apart).
  G2 spectrum: the latest-level centered hermitian descendant has a simple
     spectrum; rho_ord and x2 are recomputed from it.
  G3 shape band with scheme rider: rho_ord must land in the quoted
     convention band of the common-scheme scan; the alternative bands are
     recorded because the band is convention dependent.
  G4 normalization clause: the payload must emit positive per-side spans.
     Absence of this clause is the operator-side restatement of the
     nonidentifiability theorem once ancestry has passed.
  G5 span agreement: emitted spans must match the version-B calibrated
     spans of the held-out artifact within the declared tolerance
     (compare-only acceptance check, not a solve input).
  G6 forward dominance: the six-mass forward map built from the candidate
     spectrum (rays and mean-law coefficients recomputed from the
     candidate rho_ord and x2, shared scale from the payload or the
     current writeback) must beat the zero-input version-C table in
     maximum absolute log deviation.

A passing run still leaves the formal step of rerunning the
nonidentifiability obstruction with the new kernel ancestor (the fiber
must collapse), which is recorded as the remaining burden rather than
executed here.

Run without arguments to evaluate the fail-closed baseline (the current
template kernel, which must block at G1):

    python3 code/particles/flavor/derive_quark_kernel_normalization_acceptance_harness.py
writes code/particles/runs/flavor/quark_kernel_normalization_acceptance_current.json.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys
from datetime import datetime, timezone

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNS = ROOT / "particles" / "runs" / "flavor"
DEFAULT_OUT = RUNS / "quark_kernel_normalization_acceptance_current.json"
KERNEL_PATH = RUNS / "family_transport_kernel.json"
WRITEBACK_PATH = RUNS / "charged_shared_absolute_scale_writeback.json"
HELDOUT_PATH = RUNS / "quark_shape_law_heldout_prediction_test.json"
SCAN_PATH = RUNS / "quark_common_scheme_shape_law_scan.json"
REFERENCE_PATH = ROOT / "particles" / "data" / "particle_reference_values.json"

sys.path.insert(0, str(HERE))
from derive_quark_sigma_source_nonidentifiability_obstruction import (  # noqa: E402
    FORBIDDEN_ANCESTORS,
)

SPAN_AGREEMENT_TOLERANCE = 0.02
UP_KEYS = ("up_quark", "charm_quark", "top_quark")
DOWN_KEYS = ("down_quark", "strange_quark", "bottom_quark")


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _decode(payload) -> np.ndarray:
    return np.asarray(payload["real"], dtype=float) + 1j * np.asarray(payload["imag"], dtype=float)


def rays(rho: float) -> tuple[np.ndarray, np.ndarray]:
    denom = 3.0 * (1.0 + rho)
    v_u = np.asarray([-(2.0 * rho + 1.0), rho - 1.0, rho + 2.0]) / denom
    v_d = np.asarray([-(rho + 2.0), 1.0 - rho, 2.0 * rho + 1.0]) / denom
    return v_u, v_d


def baseline_template_candidate() -> dict:
    """The current template kernel, wrapped as a candidate payload."""
    kernel = json.loads(KERNEL_PATH.read_text(encoding="utf-8"))
    return {
        "candidate_id": "current_template_kernel_baseline",
        "kernel": {"refinements": [
            {"level": item["level"],
             "hermitian_descendant": item["hermitian_descendant"]}
            for item in kernel["refinements"]]},
        "ancestry": {
            "artifacts": [str(kernel.get("artifact"))],
            "attestations": {
                "quark_reference_values_consumed": False,
                "fitted_spreads_consumed": False,
                "numerical_flavor_template_consumed": True,
            },
        },
        "baseline_note": (
            "the handwritten numerical flavor template is declared in ancestry, "
            "so the expected outcome is an immediate block at G1"
        ),
    }


def evaluate(candidate: dict) -> dict:
    gates = {}
    first_failed = None

    def record(name: str, passed: bool, detail: dict) -> bool:
        nonlocal first_failed
        gates[name] = {"passed": bool(passed), **detail}
        if not passed and first_failed is None:
            first_failed = name
        return bool(passed)

    # G0 schema
    has_kernel = bool(candidate.get("kernel", {}).get("refinements"))
    has_ancestry = "ancestry" in candidate
    record("G0_schema", has_kernel and has_ancestry, {
        "kernel_refinements_present": has_kernel,
        "ancestry_present": has_ancestry,
    })
    if first_failed:
        return {"gates": gates, "first_failed_gate": first_failed,
                "acceptance_passed": False}

    # G1 ancestry
    ancestry = candidate["ancestry"]
    attestations = dict(ancestry.get("attestations", {}))
    target_attestations = dict(attestations)
    clean_attestations = not any(bool(v) for v in target_attestations.values())
    names = [str(a) for a in ancestry.get("artifacts", [])]
    forbidden_hits = sorted(
        {name for name in names
         for bad in FORBIDDEN_ANCESTORS if bad in name})
    record("G1_ancestry", clean_attestations and not forbidden_hits, {
        "target_attestations": target_attestations,
        "template_ancestry_declared": bool(
            attestations.get("numerical_flavor_template_consumed", False)),
        "forbidden_ancestor_hits": forbidden_hits,
        "forbidden_list_source":
            "derive_quark_sigma_source_nonidentifiability_obstruction",
    })
    if first_failed:
        return {
            "gates": gates,
            "first_failed_gate": first_failed,
            "acceptance_passed": False,
            "status": "blocked_fail_closed",
            "remaining_formal_burden": (
                "replace the handwritten numerical flavor template with a source-emitted "
                "kernel before any spectrum or normalization gate is evaluated"
            ),
        }

    # G2 spectrum
    latest = candidate["kernel"]["refinements"][-1]
    hermitian = _decode(latest["hermitian_descendant"])
    centered = hermitian - (np.trace(hermitian) / 3.0) * np.eye(3, dtype=complex)
    evals = np.linalg.eigvalsh(centered)
    g21, g32 = float(evals[1] - evals[0]), float(evals[2] - evals[1])
    simple = min(g21, g32) > 1.0e-12
    rho = 3.0 * g32 / (2.0 * g32 + g21) if simple else float("nan")
    x2 = (2.0 * g21 / float(evals[2] - evals[0]) - 1.0) if simple else float("nan")
    record("G2_spectrum", simple, {
        "eigenvalues": [float(v) for v in evals.tolist()],
        "rho_ord": rho,
        "x2": x2,
        "simple_spectrum": simple,
    })

    # G3 shape band with scheme rider
    scan = json.loads(SCAN_PATH.read_text(encoding="utf-8"))
    bands = scan["findings"]["a2_band_by_convention"]
    quoted_band = bands["quoted_mixed"]
    inside = bool(simple and quoted_band[0] <= rho <= quoted_band[1])
    record("G3_shape_band_scheme_rider", inside, {
        "declared_convention": "quoted_mixed",
        "band": quoted_band,
        "alternative_bands": {k: v for k, v in bands.items()
                              if k != "quoted_mixed"},
        "rho_ord": rho,
        "scheme_rider": "the band is convention dependent; see the "
                        "common-scheme scan artifact",
    })

    # G4 normalization clause
    clause = candidate.get("normalization_clause") or {}
    sigma_u = clause.get("sigma_u_per_side")
    sigma_d = clause.get("sigma_d_per_side")
    clause_ok = (isinstance(sigma_u, (int, float)) and sigma_u > 0.0
                 and isinstance(sigma_d, (int, float)) and sigma_d > 0.0)
    record("G4_normalization_clause_present", clause_ok, {
        "clause_kind": clause.get("kind"),
        "sigma_u_per_side": sigma_u,
        "sigma_d_per_side": sigma_d,
        "reading_if_failed": "no normalization clause means the kernel "
                             "fixes only scale-invariant shape data, which "
                             "is the operator-side restatement of the "
                             "nonidentifiability theorem",
    })

    heldout = json.loads(HELDOUT_PATH.read_text(encoding="utf-8"))
    if clause_ok:
        # G5 span agreement against the version-B calibrated spans
        cal = heldout["version_b_two_in_four_out"]["calibrated"]
        dev_u = sigma_u / float(cal["sigma_u"]) - 1.0
        dev_d = sigma_d / float(cal["sigma_d"]) - 1.0
        record("G5_span_agreement",
               abs(dev_u) <= SPAN_AGREEMENT_TOLERANCE
               and abs(dev_d) <= SPAN_AGREEMENT_TOLERANCE, {
                   "reference": "version_b_two_in_four_out.calibrated "
                                "(compare-only acceptance check)",
                   "sigma_u_relative_deviation": dev_u,
                   "sigma_d_relative_deviation": dev_d,
                   "tolerance": SPAN_AGREEMENT_TOLERANCE,
               })

        # G6 forward dominance over the zero-input version-C table
        a_ud = 1.0 / (2.0 * (1.0 + rho - x2 * x2))
        b_ud = 1.0 / (2.0 * (1.0 - x2 * x2 - x2 * x2 / (1.0 + rho)))
        g_ch_block = candidate.get("g_ch_replacement") or {}
        if "value" in g_ch_block:
            g_ch = float(g_ch_block["value"])
            g_ch_source = "candidate_replacement"
        else:
            writeback = json.loads(WRITEBACK_PATH.read_text(encoding="utf-8"))
            g_ch = float(writeback["stored_shared_absolute_scale"])
            g_ch_source = "current_family_writeback (conditional layer "
            g_ch_source += "retained)"
        seed = 0.5 * (sigma_u + sigma_d)
        eta = 0.5 * (sigma_u - sigma_d)
        g_u = g_ch * math.exp(-(a_ud * seed - b_ud * eta))
        g_d = g_ch * math.exp(-(a_ud * seed + b_ud * eta))
        v_u, v_d = rays(rho)
        refs = json.loads(REFERENCE_PATH.read_text(encoding="utf-8"))["entries"]
        predicted = {}
        for keys, g, sigma, ray in ((UP_KEYS, g_u, sigma_u, v_u),
                                    (DOWN_KEYS, g_d, sigma_d, v_d)):
            for key, v in zip(keys, ray, strict=True):
                predicted[key] = g * math.exp(2.0 * sigma * float(v))
        dlogs = {key: abs(math.log(predicted[key]
                                   / float(refs[key]["value_gev"])))
                 for key in predicted}
        version_c_max = max(
            row["abs_dlog"] for row in
            heldout["version_c_zero_quark_inputs"]["comparison"].values())
        candidate_max = max(dlogs.values())
        record("G6_forward_dominance", candidate_max < version_c_max, {
            "g_ch_source": g_ch_source,
            "candidate_max_abs_dlog": candidate_max,
            "version_c_max_abs_dlog": version_c_max,
            "candidate_abs_dlog_by_mass": dlogs,
            "candidate_predicted_gev": predicted,
        })

    passed = first_failed is None
    return {
        "gates": gates,
        "first_failed_gate": first_failed,
        "acceptance_passed": passed,
        "status": ("acceptance_passed_promotion_review_unlocked" if passed
                   else "blocked_fail_closed"),
        "remaining_formal_burden": (
            "rerun derive_quark_sigma_source_nonidentifiability_obstruction "
            "with the accepted kernel as ancestor; the spread fiber must "
            "collapse from (R_>0)^2 to a unique tuple before any promotion"),
    }


def build(candidate: dict) -> dict:
    result = evaluate(candidate)
    return {
        "artifact": "oph_quark_kernel_normalization_acceptance",
        "generated_utc": _timestamp(),
        "github_issues": [377, 379, 380],
        "row_class": "acceptance_harness_record",
        "guards": {
            "public_promotion_allowed": False,
            "acceptance_is_not_promotion": True,
            "legacy_reciprocal_ray_subfamily_only": True,
            "template_ancestry_is_fatal": True,
            "physical_yukawa_closure_possible_from_this_harness": False,
            "measured_values_role": "compare-only acceptance checks "
                                    "(G5, G6); never solve inputs",
        },
        "candidate_id": candidate.get("candidate_id"),
        "candidate_baseline_note": candidate.get("baseline_note"),
        "span_agreement_tolerance": SPAN_AGREEMENT_TOLERANCE,
        **result,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Evaluate a candidate kernel against the issue-377 "
                    "acceptance gates.")
    parser.add_argument("--candidate", default=None,
                        help="path to a candidate payload JSON; omit to "
                             "evaluate the fail-closed template baseline")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    if args.candidate:
        candidate = json.loads(
            pathlib.Path(args.candidate).read_text(encoding="utf-8"))
    else:
        candidate = baseline_template_candidate()

    report = build(candidate)
    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")

    print(f"candidate: {report['candidate_id']}")
    for name, gate in report["gates"].items():
        print(f"  {name:34s} passed={gate['passed']}")
    print(f"first failed gate: {report['first_failed_gate']}")
    print(f"status: {report['status']}")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
