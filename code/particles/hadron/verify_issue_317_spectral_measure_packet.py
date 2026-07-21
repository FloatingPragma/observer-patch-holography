#!/usr/bin/env python3
"""Verifier for OPH issue #317: Ward-projected hadronic spectral-measure proof packet.

The packet certifies the construction layer that the production hadron
backend must satisfy: the renormalized electromagnetic hadronic current,
its positive spectral measure with thresholds, resonances-and-continuum
typing, covariance, and uncertainty propagation, all derived from the OPH
QCD/matter sector with no measured HVP input.

The packet is upstream of production execution: the production backend
(issue #425) "exports the Ward-current spectral payload consumed by #317",
so this verifier fixes and machine-checks the contract that the export
must satisfy, proves the structural theorem parts on the declared
construction, and demonstrates the constructor end-to-end on committed
real lattice data at diagnostic scale. Physical promotion of any hadronic
number stays with the production lane.

Machine checks executed at emission time:

1. Gauge-invariance witness: the U(1) vector correlator and the plaquette
   are recomputed on a random SU(3) gauge transform of a rough interacting
   background and must match the untransformed values.
2. Ward-identity witness: the exact backward-difference divergence of the
   conserved point-split current vanishes off the source up to solver
   tolerance on the same rough background.
3. Normalization witness: the conserved-local estimator
   Z_V^eff = C_CL/(2 kappa C_LL) plateaus at 1 on the free field.
4. Positivity/covariance demonstrator: a machine-readable single-atom
   spectral measure is emitted from the committed real quenched ensemble
   (hybrid_ir_ensembleA_2026-07-16.npz): ground vector level, Ward-projected
   residue, jackknife covariance with a positive-semidefiniteness check,
   and a typed budget ledger. Bare parameters only; no measured hadronic
   input of any kind enters.
5. Acceptance gate: a schema-conformant synthetic production payload is
   accepted and a battery of negative controls (measured-HVP guard flips,
   surrogate flag, missing budget, quenched branch, negative residue,
   empty level support, forbidden-target leak) is rejected, so the
   machine-readable emission fails closed exactly as the claim requires.
6. Empirical typing and blinding: the source-side comparison manifest is
   empty (blinding intact, no comparison event), the production schema and
   the executed gate exclude every measured-HVP path from the source
   construction, the declared-empirical companion schema pins the
   comparison row class non-source by constants, and the public-use
   property of that empirical surface is recorded as an explicit
   disclosure confined to the non-promotable row class.
7. Higher-point and radiative typing: the four-current receipt records the
   two-point measure as insufficient (higher-point half), and the
   Xi_same_scheme remainder ledger is a distinct typed artifact whose
   required terms include the QED and EW rows while the production schema
   excludes QED from the two-point branch (radiative half); the Xi source
   certificate itself is a disclosed production-level gap.

The certificate separates two verdicts. ``accepted`` covers the packet
construction/contract layer only. ``issue_closure_condition`` computes the
author-stated closure bar (backend export bundle, unquenching execution,
runtime receipt, uncertainty ledger) from the live production-lane
artifacts and currently reports ``met_locally: false``; the production
export bundle, populated base measure, and Ward-current source certificate
remain open in the production lane.

Run:
    python3 code/particles/hadron/verify_issue_317_spectral_measure_packet.py \
        --check --output \
        code/particles/runs/hadron/ward_projected_spectral_measure_proof_packet.json
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from lattice_backend.core import (  # noqa: E402
    average_plaquette,
    cold_start,
    gauge_transform,
    random_su3_near_identity,
    sweep,
)
from lattice_backend.dirac import WilsonClover, point_propagator  # noqa: E402
from lattice_backend.conserved_vector import (  # noqa: E402
    conserved_local_correlator,
    transverse_vector_correlator,
    ward_divergence_offsource_max,
    zv_effective,
)
from lattice_backend.vector_correlator import fold_correlator  # noqa: E402

ROOT = HERE.parents[1]
SCHEMA_PATH = HERE / "ward_projected_spectral_measure.schema.json"
EMPIRICAL_SCHEMA_PATH = HERE / "empirical_ward_projected_spectral_measure.schema.json"
ENSEMBLE_NPZ = ROOT / "particles" / "runs" / "hadron" / "hybrid_ir_ensembleA_2026-07-16.npz"
DEFAULT_OUT = ROOT / "particles" / "runs" / "hadron" / "ward_projected_spectral_measure_proof_packet.json"
BACKEND_DIR = ROOT / "particles" / "runs" / "qcd" / "hadron_source_backend"
COMPARISON_MANIFEST = BACKEND_DIR / "controls" / "comparison_data_manifest.json"
Q4_RECEIPT = BACKEND_DIR / "higher_point" / "Q4_HLbL_receipt.json"
XI_LEDGER = BACKEND_DIR / "endpoint" / "Xi_same_scheme.json"
BASE_MEASURE = BACKEND_DIR / "qcd_ensemble" / "base_measure.json"
WARD_CURRENT_DEFINITION = BACKEND_DIR / "currents" / "ward_current_definition.json"
READINESS_REPORT = ROOT / "particles" / "runs" / "hadron" / "hadron_production_readiness_report.json"

FORBIDDEN_TARGETS = [
    "CODATA_ALPHA",
    "MUON_G_MINUS_2",
    "EE_TO_HADRONS",
    "RARE_DECAY_DATA",
    "HADRON_MASS_TARGETS",
    "PDG_QCD_FITS",
]

GAUGE_INVARIANCE_TOL = 1e-8
WARD_RELATIVE_TOL = 1e-9
ZV_FREE_FIELD_TOL = 5e-3
EFFMASS_WINDOW_START = 3


def _now_utc() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _rel(path: Path) -> str:
    return path.relative_to(ROOT.parent).as_posix()


# ---------------------------------------------------------------------------
# Live artifact checks: typing criteria and the issue-closure condition.
# ---------------------------------------------------------------------------


def empirical_typing_check(gate: dict[str, Any]) -> dict[str, Any]:
    """Criterion: empirical data are used only for a blinded comparison.

    Machine-checked in three layers against live artifacts:

    1. Blinding of the source construction: no empirical dataset is attached
       to the source lane (comparison manifest empty), and the source export
       cannot admit empirical input (production schema pins the source
       guards as constants and the semantic gate fails closed on
       EE_TO_HADRONS; executed in the acceptance-gate witness).
    2. Comparison channel typing: the only surface where empirical data may
       ever attach to this lane is the declared-empirical companion row
       class, pinned non-source by schema constants.
    3. Disclosure: the empirical companion surface allows public final
       values (usable_for_public_final_values const true). That use is
       confined to the disclosed oph_plus_empirical_hadron_closure row
       class, which by the same constants can never be promoted as an OPH
       source theorem and never satisfies the production constructive
       artifact, so the source measure itself only ever meets empirical
       data as a comparison. The check records the disclosure explicitly
       rather than hiding it.

    No comparison event has occurred (the manifest is empty); the
    unblinding event is production-scope and must attach through the
    comparison manifest after a source payload exists.
    """
    manifest = _load_json(COMPARISON_MANIFEST)
    schema = _load_json(EMPIRICAL_SCHEMA_PATH)
    production_schema = _load_json(SCHEMA_PATH)
    guards = schema["properties"]["guards"]["properties"]
    production_guards = production_schema["properties"]["guards"]["properties"]
    gate_rejections = {
        row["control_id"]: row["rejected"] for row in gate["negative_controls"]
    }
    checks = {
        "comparison_manifest_empty": manifest.get("comparison_data") == [],
        "comparison_manifest_status_no_data": manifest.get("status") == "NO_COMPARISON_DATA_ATTACHED",
        "production_schema_rejects_compare_only_const_false": (
            production_guards["compare_only_external_endpoint"].get("const") is False
        ),
        "production_schema_rejects_surrogate_const_false": (
            production_guards["surrogate_hadron_artifact"].get("const") is False
        ),
        "gate_rejects_measured_hvp_comparison_endpoint": bool(
            gate_rejections.get("measured_hvp_comparison_endpoint")
        ),
        "gate_rejects_ee_to_hadrons_target_leak": bool(
            gate_rejections.get("forbidden_target_leak_ee_to_hadrons")
        ),
        "row_class_const_empirical_closure": (
            schema["properties"]["row_class"].get("const") == "oph_plus_empirical_hadron_closure"
        ),
        "empirical_surface_not_source_only_const": guards["source_only"].get("const") is False,
        "empirical_use_disclosed_const_true": (
            guards["external_cross_section_data_used"].get("const") is True
        ),
        "promotable_as_oph_source_theorem_const_false": (
            guards["promotable_as_oph_source_theorem"].get("const") is False
        ),
        "satisfies_production_constructive_next_artifact_const_false": (
            guards["satisfies_production_constructive_next_artifact"].get("const") is False
        ),
    }
    public_use_disclosure = {
        "usable_for_public_final_values_const": guards["usable_for_public_final_values"].get("const"),
        "confinement": (
            "public use is a property of the disclosed empirical-closure row class "
            "only; the same schema constants forbid promoting that surface as an "
            "OPH source theorem or substituting it for the production artifact, "
            "so with respect to the source spectral measure the empirical data "
            "remain comparison-only"
        ),
    }
    blinding_status = {
        "comparison_event_performed": bool(manifest.get("comparison_data")),
        "statement": (
            "no empirical dataset is attached to the source lane, so no comparison "
            "has occurred and there is nothing unblinded; the comparison event is "
            "production-scope and must append to the comparison manifest after the "
            "source payload is emitted"
        ),
    }
    return {
        "checked_artifacts": {
            "comparison_manifest": _rel(COMPARISON_MANIFEST),
            "empirical_companion_schema": _rel(EMPIRICAL_SCHEMA_PATH),
            "production_export_schema": _rel(SCHEMA_PATH),
        },
        "checks": checks,
        "public_use_disclosure": public_use_disclosure,
        "blinding_status": blinding_status,
        "passed": bool(all(checks.values())),
    }


def higher_point_typing_check() -> dict[str, Any]:
    """Criterion: higher-point and radiative corrections stay separately typed.

    Machine-checked for both halves against live artifacts:

    - higher-point half: the four-current receipt records the two-point
      measure as insufficient for HLbL-class objects, non-promoting, with no
      external targets;
    - radiative half: the same-scheme remainder ledger Xi_same_scheme is a
      distinct typed artifact whose required term rows include the QED and
      EW corrections, non-promoting, with no external targets; the
      production export schema additionally excludes QED from the two-point
      branch (branch.qed limited to off/explicitly_superseded).

    The Xi ledger's source certificate is itself still missing
    (MISSING_SOURCE_CERTIFICATE); that is a production-level content gap
    recorded under production_open_items. The typing separation checked
    here is exactly that the radiative rows live in that distinct ledger
    rather than inside the two-point measure.
    """
    receipt = _load_json(Q4_RECEIPT)
    xi = _load_json(XI_LEDGER)
    production_schema = _load_json(SCHEMA_PATH)
    qed_enum = (
        production_schema["properties"]["branch"]["properties"]["qed"].get("enum") or []
    )
    checks = {
        "q4_status_two_point_measure_insufficient": (
            receipt.get("status") == "TWO_POINT_MEASURE_INSUFFICIENT"
        ),
        "q4_promotion_not_allowed": receipt.get("promotion_allowed") is False,
        "q4_no_external_targets": receipt.get("external_targets_used") == [],
        "xi_ledger_distinct_artifact": xi.get("artifact") == "Xi_same_scheme",
        "xi_required_terms_include_qed_and_ew": (
            {"QED", "EW"} <= set(xi.get("required_terms", []))
        ),
        "xi_promotion_not_allowed": xi.get("promotion_allowed") is False,
        "xi_no_external_targets": xi.get("external_targets_used") == [],
        "production_branch_excludes_qed": "on" not in qed_enum and "off" in qed_enum,
    }
    return {
        "checked_artifacts": {
            "q4_hlbl_receipt": _rel(Q4_RECEIPT),
            "xi_same_scheme_ledger": _rel(XI_LEDGER),
            "production_export_schema": _rel(SCHEMA_PATH),
        },
        "checks": checks,
        "xi_ledger_source_certificate_status": xi.get("status"),
        "passed": bool(all(checks.values())),
    }


def issue_closure_condition() -> dict[str, Any]:
    """Author-stated closure bar, computed from live repo state.

    The issue author's close condition requires the backend export bundle,
    unquenching execution, runtime receipt, and uncertainty ledger. This
    packet certifies the construction/contract layer; the closure objects
    live in the production lane, and their current statuses are read from
    the repository rather than asserted.
    """
    readiness = _load_json(READINESS_REPORT)
    base_measure = _load_json(BASE_MEASURE)
    ward_current = _load_json(WARD_CURRENT_DEFINITION)
    remaining = readiness["exact_remaining_runtime_object"]
    closure = readiness["closure_status"]
    checked = {
        "production_backend_export_bundle": {
            "path": _rel(READINESS_REPORT),
            "status": remaining.get("status"),
            "closure_grade": closure.get("closure_grade"),
            "public_unsuppression_ready": closure.get("public_unsuppression_ready"),
        },
        "qcd_base_measure": {
            "path": _rel(BASE_MEASURE),
            "status": base_measure.get("status"),
        },
        "ward_current_source_certificate": {
            "path": _rel(WARD_CURRENT_DEFINITION),
            "status": ward_current.get("status"),
        },
    }
    met = bool(
        remaining.get("status") != "open"
        and closure.get("closure_grade") not in (None, "execution_incomplete")
        and base_measure.get("status") != "REQUIRED_NOT_POPULATED"
        and ward_current.get("status") != "MISSING_SOURCE_CERTIFICATE"
    )
    return {
        "author_condition": (
            "close only when the backend export bundle, unquenching execution, "
            "runtime receipt, and uncertainty ledger are present"
        ),
        "checked_objects": checked,
        "met_locally": met,
        "resolution_path": (
            "production backend export lane (dependency_note); the export must "
            "conform to the schema and acceptance gate certified in this packet"
        ),
    }


# ---------------------------------------------------------------------------
# Executed lattice witnesses (small, seeded, deterministic).
# ---------------------------------------------------------------------------


def gauge_and_ward_witness(
    shape: tuple[int, int, int, int] = (8, 2, 2, 2),
    beta: float = 5.5,
    kappa: float = 0.12,
    n_sweeps: int = 6,
    cg_tol: float = 1e-12,
    seed: int = 3170,
) -> dict[str, Any]:
    """Gauge invariance and exact Ward identity on one rough background."""
    rng = np.random.default_rng(seed)
    u = cold_start(shape)
    for _ in range(n_sweeps):
        sweep(rng, u, beta, n_or=1)
    op = WilsonClover(u, kappa=kappa, c_sw=1.0)
    prop, _ = point_propagator(op, shape, tol=cg_tol)
    g_ll = transverse_vector_correlator(prop)
    plaq = average_plaquette(u)

    g_field = random_su3_near_identity(rng, shape, eps=1.0)
    u_rot = gauge_transform(u, g_field)
    op_rot = WilsonClover(u_rot, kappa=kappa, c_sw=1.0)
    prop_rot, _ = point_propagator(op_rot, shape, tol=cg_tol)
    g_ll_rot = transverse_vector_correlator(prop_rot)
    plaq_rot = average_plaquette(u_rot)

    corr_scale = float(np.max(np.abs(g_ll)))
    corr_defect = float(np.max(np.abs(g_ll_rot - g_ll))) / corr_scale
    plaq_defect = abs(plaq_rot - plaq)

    ward_defect, ward_scale = ward_divergence_offsource_max(prop, op.ubc, kappa, nu=1)
    ward_relative = ward_defect / ward_scale

    return {
        "background": {
            "lattice_shape_TXYZ": list(shape),
            "beta": beta,
            "kappa": kappa,
            "n_sweeps": n_sweeps,
            "cg_tol": cg_tol,
            "rng_seed": seed,
            "plaquette": plaq,
        },
        "gauge_invariance": {
            "statement": (
                "the transverse U(1) vector correlator and the plaquette are invariant "
                "under U_mu(x) -> g(x) U_mu(x) g(x+mu)^dag for random SU(3) g"
            ),
            "correlator_relative_defect": corr_defect,
            "plaquette_absolute_defect": plaq_defect,
            "tolerance": GAUGE_INVARIANCE_TOL,
            "passed": bool(corr_defect < GAUGE_INVARIANCE_TOL and plaq_defect < GAUGE_INVARIANCE_TOL),
        },
        "ward_identity": {
            "statement": (
                "the conserved point-split current V_mu satisfies the exact backward-difference "
                "Ward identity sum_mu [V_mu(x) - V_mu(x-mu)] = 0 off the contact point"
            ),
            "max_offsource_divergence": ward_defect,
            "correlator_scale": ward_scale,
            "relative_defect": ward_relative,
            "tolerance": WARD_RELATIVE_TOL,
            "passed": bool(ward_relative < WARD_RELATIVE_TOL),
        },
    }


def zv_free_field_witness(
    shape: tuple[int, int, int, int] = (16, 2, 2, 2),
    kappa: float = 0.13,
    cg_tol: float = 1e-11,
) -> dict[str, Any]:
    """Free-field anchor: Z_V^eff = C_CL/(2 kappa C_LL) plateaus at exactly 1."""
    u = cold_start(shape)
    op = WilsonClover(u, kappa=kappa, c_sw=1.0)
    prop, _ = point_propagator(op, shape, tol=cg_tol)
    g_ll = transverse_vector_correlator(prop)
    g_cl = conserved_local_correlator(prop, op.ubc, kappa)
    zv = zv_effective(g_cl, g_ll, kappa)
    plateau = zv[4:8]
    max_dev = float(np.max(np.abs(plateau - 1.0)))
    return {
        "lattice_shape_TXYZ": list(shape),
        "kappa": kappa,
        "statement": (
            "on the free field the local-current renormalization estimator "
            "Z_V^eff(t) = C_CL(t)/(2 kappa C_LL(t)) equals 1, anchoring sign and normalization"
        ),
        "plateau_window_t": [4, 7],
        "plateau_values": [float(v) for v in plateau],
        "max_abs_deviation_from_one": max_dev,
        "tolerance": ZV_FREE_FIELD_TOL,
        "passed": bool(max_dev < ZV_FREE_FIELD_TOL),
    }


# ---------------------------------------------------------------------------
# Demonstrator measure from the committed real quenched ensemble.
# ---------------------------------------------------------------------------


def _cosh_effective_mass(t: int, ratio: float, half: int) -> float:
    """Solve cosh(m*(T/2 - t))/cosh(m*(T/2 - t - 1)) = ratio by bisection."""
    if ratio <= 1.0:
        return float("nan")
    lo, hi = 1e-8, 50.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        val = math.cosh(mid * (half - t)) / math.cosh(mid * (half - t - 1))
        if val > ratio:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def _effmass_window(folded: np.ndarray, start: int = EFFMASS_WINDOW_START) -> list[int]:
    half = len(folded) - 1
    ts: list[int] = []
    for t in range(start, half):
        if folded[t] > 0.0 and folded[t + 1] > 0.0 and folded[t] > folded[t + 1]:
            ts.append(t)
        elif ts:
            break
    return ts


def _windowed_mass(folded: np.ndarray, start: int = EFFMASS_WINDOW_START) -> float:
    half = len(folded) - 1
    vals = [
        _cosh_effective_mass(t, folded[t] / folded[t + 1], half)
        for t in _effmass_window(folded, start)
    ]
    vals = [v for v in vals if math.isfinite(v)]
    return float(np.mean(vals)) if vals else float("nan")


def _windowed_zv(folded_cl: np.ndarray, folded_ll: np.ndarray, kappa: float) -> float:
    half = len(folded_ll) - 1
    vals = [
        float(folded_cl[t] / (2.0 * kappa * folded_ll[t]))
        for t in range(EFFMASS_WINDOW_START, half)
        if folded_ll[t] != 0.0
    ]
    return float(np.mean(vals))


def _atom_pipeline(
    ll_mean: np.ndarray, cl_mean: np.ndarray, kappa: float, t_extent: int
) -> dict[str, float]:
    """Single-atom constructor: level, Z_V, Ward-projected residue weight.

    Single-state model for the folded correlator, f(t) = A (e^{-m t} +
    e^{-m (T-t)}).  For an R-normalized spectral atom rho_R = w delta(s - m^2)
    the amplitude is A = w m / (24 pi^2)
    (``vector_correlator.synthetic_atom_correlator``), so w_hop = 24 pi^2 A / m.
    The correlator is built from hopping-normalized propagators; the physical
    correlator is (2 kappa)^2 times the hopping one (psi_phys =
    sqrt(2 kappa) psi_hop; ``conserved_vector.py`` docstring: "the analysis
    layer applies that factor"), and the local current renormalizes with
    Z_V^2: w_phys = Z_V^2 (2 kappa)^2 w_hop.
    """
    folded_ll = fold_correlator(ll_mean)
    folded_cl = fold_correlator(cl_mean)
    m_v = _windowed_mass(folded_ll)
    z_v = _windowed_zv(folded_cl, folded_ll, kappa)
    window = _effmass_window(folded_ll)
    amps = [
        folded_ll[t] / (math.exp(-m_v * t) + math.exp(-m_v * (t_extent - t)))
        for t in window
    ]
    amp = float(np.mean(amps))
    w_hop = 24.0 * math.pi**2 * amp / m_v
    w_phys = z_v * z_v * (2.0 * kappa) ** 2 * w_hop
    return {"am_v": m_v, "z_v": z_v, "amplitude_hop": amp, "weight_phys": w_phys}


def demonstrator_measure(npz_path: Path = ENSEMBLE_NPZ) -> dict[str, Any]:
    """Machine-readable single-atom measure from the committed real ensemble."""
    d = np.load(npz_path)
    kappa = float(d["kappa"])
    t_extent = int(d["shape"][0])
    g_ll = np.asarray(d["g_ll"])
    g_cl = np.asarray(d["g_cl"])
    n_cfg = g_ll.shape[0]

    central = _atom_pipeline(g_ll.mean(axis=0), g_cl.mean(axis=0), kappa, t_extent)

    samples: dict[str, list[float]] = {k: [] for k in central}
    folded_samples = []
    for i in range(n_cfg):
        ll_i = np.delete(g_ll, i, axis=0).mean(axis=0)
        cl_i = np.delete(g_cl, i, axis=0).mean(axis=0)
        res = _atom_pipeline(ll_i, cl_i, kappa, t_extent)
        for key, value in res.items():
            samples[key].append(value)
        folded_samples.append(fold_correlator(ll_i))
    errors = {}
    for key, vals in samples.items():
        arr = np.array(vals)
        center = arr.mean()
        errors[key] = float(np.sqrt((n_cfg - 1) / n_cfg * np.sum((arr - center) ** 2)))

    folded_arr = np.array(folded_samples)
    folded_center = folded_arr.mean(axis=0)
    dev = folded_arr - folded_center
    covariance = (n_cfg - 1) / n_cfg * dev.T @ dev
    eigvals = np.linalg.eigvalsh(covariance)
    cov_scale = float(np.max(np.abs(covariance)))
    cov_psd = bool(float(eigvals.min()) >= -1e-12 * cov_scale)

    # Positivity is a statement about the exact expectation value; the
    # finite-sample estimate is checked statistically: every timeslice must be
    # positive where significant (|G| > 2 sigma) and consistent with zero
    # where not. A significantly negative value fails the witness.
    folded_mean = fold_correlator(g_ll.mean(axis=0))
    folded_err = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    significant = np.abs(folded_mean) > 2.0 * folded_err
    positive_where_significant = bool(np.all(folded_mean[significant] > 0.0))
    negative_only_within_noise = bool(
        np.all((folded_mean >= 0.0) | (np.abs(folded_mean) <= 2.0 * folded_err))
    )
    correlator_positive = positive_where_significant and negative_only_within_noise
    weight_positive = bool(central["weight_phys"] > 0.0)
    positivity_report = {
        "rule": (
            "Kallen-Lehmann positivity holds for the exact expectation value; the "
            "finite-sample witness requires G(t) > 0 on every statistically "
            "significant timeslice (|G| > 2 sigma) and |G| <= 2 sigma wherever the "
            "central value is negative"
        ),
        "timeslices": [
            {
                "t": int(t),
                "value": float(folded_mean[t]),
                "jackknife_error": float(folded_err[t]),
                "significant": bool(significant[t]),
            }
            for t in range(len(folded_mean))
        ],
        "positive_where_significant": positive_where_significant,
        "negative_values_consistent_with_zero_at_2_sigma": negative_only_within_noise,
        "noise_note": (
            "the vector channel enters noise at t >= 6 on this volume (declared in "
            "the hybrid IR bracket note); negative central values there are "
            "statistical fluctuations around zero, not spectral-positivity violations"
        ),
    }

    # quadrature/window budget: shift the extraction window start by +1
    window_shift_mass = abs(
        _windowed_mass(folded_mean, EFFMASS_WINDOW_START + 1) - central["am_v"]
    )

    unquantified = {
        "status": "not_quantified_at_diagnostic_scale",
        "production_scope": "production backend export budgets (dependency_note)",
    }
    am_v = central["am_v"]

    return {
        "artifact": "oph_diagnostic_ward_projected_spectral_measure_demonstrator",
        "role": (
            "executed witness that the packet constructor emits a machine-readable "
            "positive spectral measure from source-only lattice data; diagnostic "
            "scale, non-promoting, not a production payload"
        ),
        "ensemble": {
            "npz": npz_path.relative_to(ROOT).as_posix(),
            "ensemble_id": str(d["ensemble_id"]),
            "lattice_shape_TXYZ": [int(v) for v in d["shape"]],
            "beta": float(d["beta"]),
            "kappa": kappa,
            "n_configs": n_cfg,
            "rng_seed": int(d["seed"]),
            "flavors": "quenched_nf0",
            "external_inputs_used": False,
            "input_note": (
                "bare parameters (beta, kappa) only; the extracted level and residue "
                "are whatever the lattice dynamics produce; no measured HVP, R(s), "
                "alpha, or hadron-mass input enters the constructor"
            ),
        },
        "units": "lattice units a = 1; no scale setting at diagnostic scale (declared)",
        "charge_convention": (
            "single flavor, unit charge (R-normalization N_c Q^2 = 3); the U(1)_Q "
            "flavor charge weighting is a declared production-side factor"
        ),
        "finite_volume_levels": [
            {
                "ensemble_id": str(d["ensemble_id"]),
                "channel": "vector_transverse_unit_charge",
                "levels": [
                    {
                        "level_id": "vector_ground_A",
                        "s": am_v * am_v,
                        "energy": am_v,
                        "energy_jackknife_error": errors["am_v"],
                        "weight": central["weight_phys"],
                        "weight_jackknife_error": errors["weight_phys"],
                    }
                ],
            }
        ],
        "ward_projected_residues": [
            {
                "level_id": "vector_ground_A",
                "residue": central["weight_phys"],
                "current_normalization": (
                    "Z_V^eff = C_CL/(2 kappa C_LL) from the exactly conserved "
                    f"point-split current; Z_V = {central['z_v']:.6f} "
                    f"+- {errors['z_v']:.6f}"
                ),
            }
        ],
        "rho_had_or_measure": {
            "representation": "primitive_spectral_measure",
            "support_variable": "s",
            "pushforward_rule": (
                "rho(s) = w delta(s - (a m_V)^2) from the single-state model of the "
                "folded transverse correlator; the production pushforward to "
                "rho_Q(s;P) is the declared Luscher-class finite-volume map "
                "(rho_scattering/ scaffolding, production backend export lane)"
            ),
            "positivity_status": (
                "verified_nonnegative: correlator positive on all significant "
                "timeslices, noise-level tail consistent with zero, atom weight positive"
                if correlator_positive and weight_positive
                else "POSITIVITY_VIOLATION"
            ),
            "positivity_report": positivity_report,
        },
        "covariance": {
            "object": "jackknife covariance of the folded transverse correlator",
            "dimension": int(covariance.shape[0]),
            "matrix": [[float(x) for x in row] for row in covariance],
            "min_eigenvalue": float(eigvals.min()),
            "max_eigenvalue": float(eigvals.max()),
            "positive_semidefinite": cov_psd,
        },
        "extraction": {
            "am_vector": am_v,
            "am_vector_jackknife_error": errors["am_v"],
            "z_v": central["z_v"],
            "z_v_jackknife_error": errors["z_v"],
            "amplitude_hop": central["amplitude_hop"],
            "weight_physical_units": central["weight_phys"],
            "weight_jackknife_error": errors["weight_phys"],
            "window_rule": "cosh effective mass over the positive decaying window from t = 3",
        },
        "systematics": {
            "statistical_budget": {
                "method": "leave-one-out jackknife through the full constructor",
                "relative_error_weight": errors["weight_phys"] / central["weight_phys"],
                "relative_error_energy": errors["am_v"] / central["am_v"],
            },
            "quadrature_budget": {
                "method": "extraction-window start shifted by +1",
                "abs_mass_shift": window_shift_mass,
            },
            "continuum_budget": unquantified,
            "finite_volume_budget": unquantified,
            "chiral_budget": unquantified,
            "current_matching_budget": unquantified,
            "endpoint_remainder_budget": unquantified,
        },
        "guards": {
            "execution_class": "real_lattice_diagnostic_toy_scale",
            "promotion_allowed": False,
            "production_schema_conformant": False,
            "production_schema_nonconformance_reason": (
                "branch.flavors is quenched_nf0; the production schema requires the "
                "seeded 2+1 family, so this demonstrator is rejected by the same "
                "acceptance gate it demonstrates"
            ),
            "surrogate_hadron_artifact": False,
            "target_anchored": False,
            "external_inputs_used": False,
        },
        "checks_passed": bool(
            correlator_positive
            and weight_positive
            and cov_psd
            and math.isfinite(am_v)
            and am_v > 0.0
        ),
    }


# ---------------------------------------------------------------------------
# Acceptance gate: schema plus semantic checks, with negative controls.
# ---------------------------------------------------------------------------


def build_conformant_payload() -> dict[str, Any]:
    """Synthetic production payload conforming to the export schema."""
    return {
        "artifact": "oph_qcd_ward_projected_hadronic_spectral_measure",
        "format_version": 1,
        "profile_id": "issue_317_acceptance_gate_synthetic_conformant_payload",
        "branch": {"flavors": "2+1", "qed": "off"},
        "projection": {
            "lane": "U(1)_Q",
            "ward_projected": True,
            "external_kinematics": "zero-momentum rest-frame finite-volume levels",
            "zero_momentum_endpoint_compatible": True,
        },
        "backend": {
            "family": "synthetic_gate_check",
            "name": "acceptance_gate",
            "version": "1",
            "git_commit": "0000000000",
            "run_id": "gate",
            "build_id": "gate",
            "machine": "gate",
        },
        "finite_volume_levels": [
            {
                "ensemble_id": "gate_ens",
                "channel": "U(1)_Q_vector",
                "levels": [
                    {"level_id": "L0", "s": 1.0, "energy": 1.0, "weight": 0.5}
                ],
            }
        ],
        "ward_projected_residues": [
            {"level_id": "L0", "residue": 0.5, "current_normalization": "Z_V ledger"}
        ],
        "rho_had_or_measure": {
            "representation": "rho_had",
            "support_variable": "s",
            "pushforward_rule": "declared finite-volume to continuum map",
            "positivity_status": "certified_positive",
        },
        "systematics": {
            "statistical_budget": {"bound": 0.01},
            "continuum_budget": {"bound": 0.01},
            "finite_volume_budget": {"bound": 0.01},
            "chiral_budget": {"bound": 0.01},
            "current_matching_budget": {"bound": 0.01},
            "quadrature_budget": {"bound": 0.01},
            "endpoint_remainder_budget": {"bound": 0.01},
        },
        "guards": {
            "stable_channel_only": False,
            "surrogate_hadron_artifact": False,
            "compare_only_external_endpoint": False,
        },
    }


def accept_payload(payload: dict[str, Any], schema: dict[str, Any]) -> tuple[bool, str]:
    """Schema plus semantic acceptance; returns (accepted, reason)."""
    try:
        import jsonschema
    except ImportError:  # pragma: no cover - jsonschema is a pinned dependency
        return False, "jsonschema_unavailable"
    try:
        jsonschema.validate(instance=payload, schema=schema)
    except jsonschema.ValidationError as exc:
        return False, f"schema_violation: {exc.message}"
    for target in payload.get("external_targets_used", []):
        if target in FORBIDDEN_TARGETS:
            return False, f"TARGET_LEAK_DETECTED: {target}"
    for block in payload["finite_volume_levels"]:
        for level in block["levels"]:
            if level["weight"] < 0.0:
                return False, f"negative_level_weight: {level['level_id']}"
            if level["energy"] <= 0.0:
                return False, f"nonpositive_energy: {level['level_id']}"
    for row in payload["ward_projected_residues"]:
        if row["residue"] < 0.0:
            return False, f"negative_residue: {row['level_id']}"
    for name, budget in payload["systematics"].items():
        if not isinstance(budget, dict) or not budget:
            return False, f"empty_budget: {name}"
    return True, "accepted"


def build_negative_controls() -> list[dict[str, Any]]:
    """Payload mutations that the gate must reject, with the criterion they guard."""
    controls = []

    def add(control_id: str, guards_criterion: str, mutate) -> None:
        payload = build_conformant_payload()
        mutate(payload)
        controls.append(
            {"control_id": control_id, "guards_criterion": guards_criterion, "payload": payload}
        )

    add(
        "measured_hvp_comparison_endpoint",
        "no measured HVP input / blinded comparison only",
        lambda p: p["guards"].__setitem__("compare_only_external_endpoint", True),
    )
    add(
        "forbidden_target_leak_ee_to_hadrons",
        "no measured HVP input",
        lambda p: p.__setitem__("external_targets_used", ["EE_TO_HADRONS"]),
    )
    add(
        "surrogate_hadron_artifact",
        "source-derived only",
        lambda p: p["guards"].__setitem__("surrogate_hadron_artifact", True),
    )
    add(
        "stable_channel_only_export",
        "full vector-channel spectral support required",
        lambda p: p["guards"].__setitem__("stable_channel_only", True),
    )
    add(
        "missing_chiral_budget",
        "uncertainty propagation ledger complete",
        lambda p: p["systematics"].pop("chiral_budget"),
    )
    add(
        "quenched_branch",
        "physical 2+1 matter sector required",
        lambda p: p["branch"].__setitem__("flavors", "quenched"),
    )
    add(
        "negative_residue",
        "positivity proved and enforced",
        lambda p: p["ward_projected_residues"][0].__setitem__("residue", -0.5),
    )
    add(
        "empty_level_support",
        "nonempty finite-volume spectral support",
        lambda p: p.__setitem__("finite_volume_levels", []),
    )
    add(
        "ward_projection_dropped",
        "gauge invariance / Ward projection required",
        lambda p: p["projection"].__setitem__("ward_projected", False),
    )
    return controls


def run_acceptance_gate() -> dict[str, Any]:
    schema = json.loads(SCHEMA_PATH.read_text(encoding="utf-8"))
    ok, reason = accept_payload(build_conformant_payload(), schema)
    control_rows = []
    all_rejected = True
    for control in build_negative_controls():
        rejected, why = accept_payload(control["payload"], schema)
        control_rows.append(
            {
                "control_id": control["control_id"],
                "guards_criterion": control["guards_criterion"],
                "rejected": not rejected,
                "rejection_reason": why,
            }
        )
        all_rejected = all_rejected and not rejected
    return {
        "schema": SCHEMA_PATH.relative_to(ROOT).as_posix(),
        "semantic_checks": [
            "forbidden-target leak fails closed (TARGET_LEAK_DETECTED)",
            "level weights and residues nonnegative",
            "level energies positive",
            "all seven budget rows present and nonempty",
        ],
        "conformant_payload_accepted": bool(ok),
        "conformant_payload_reason": reason,
        "negative_controls": control_rows,
        "all_negative_controls_rejected": bool(all_rejected),
        "passed": bool(ok and all_rejected),
    }


# ---------------------------------------------------------------------------
# Packet assembly.
# ---------------------------------------------------------------------------


def build_theorem() -> dict[str, Any]:
    return {
        "theorem_id": "WardProjectedHadronicSpectralMeasureConstruction_Q",
        "statement": (
            "On the OPH QCD/matter sector with a declared quotient ensemble, "
            "reflection-positive Euclidean vacuum, and transfer operator, the "
            "point-split U(1)_Q vector current is gauge invariant and exactly "
            "conserved, its Ward projection selects the transverse channel "
            "Pi_Q^{mu nu} = (q^mu q^nu - q^2 eta^{mu nu}) Pi_Q(q^2), and the "
            "zero-momentum transverse correlator carries the positive finite-volume "
            "spectral measure rho_Q(s) = sum_n Z_n delta(s - E_n^2) with Z_n >= 0, "
            "support starting at the lowest vector-channel state, normalization "
            "fixed by the Z_V ledger, scheme dependence confined to the declared "
            "once-subtracted transport kernel, and jackknife-through-pipeline "
            "covariance and budget-ledger uncertainty propagation."
        ),
        "parts": {
            "a_gauge_invariance_and_source_derivation": (
                "each hop term psibar(x+mu)(1 +- gamma_mu) U_mu(x)^(dag) psi(x) of the "
                "point-split current is a color trace of gauge-covariant bilinears with the "
                "compensating link insertion, hence invariant under U_mu(x) -> "
                "g(x) U_mu(x) g(x+mu)^dag, psi(x) -> g(x) psi(x); correlators are "
                "functionals of the source ensemble only (executed witness 1)"
            ),
            "b_ward_identity": (
                "the exact Noether current of the U(1) phase rotation of the "
                "hopping-normalized Wilson-clover action satisfies the backward-difference "
                "Ward identity off contact points; the clover term is phase-invariant, so "
                "c_SW does not modify the identity (executed witness 2; "
                "lattice_backend/conserved_vector.py)"
            ),
            "c_positivity_and_thresholds": (
                "by reflection positivity and the transfer-operator spectral decomposition, "
                "G(t) = sum_n w_n (e^{-E_n t} + e^{-E_n (T-t)}) with w_n proportional to "
                "|<0|J_Q|n>|^2 through a positive normalization, hence w_n >= 0; rho_Q is a "
                "positive locally finite measure supported on the vector-channel spectrum "
                "above the vacuum, with threshold at the lowest such state (two-pion "
                "channel in infinite volume)"
            ),
            "d_normalization": (
                "the exactly conserved point-split current needs no renormalization "
                "(Z_V = 1); the local current is renormalized by the conserved-local "
                "estimator Z_V^eff = C_CL/(2 kappa C_LL), anchored at exactly 1 on the "
                "free field (executed witness 3)"
            ),
            "e_scheme_dependence": (
                "the spectral measure is scheme-independent spectral data of the conserved "
                "current; scheme enters only through the declared once-subtracted transport "
                "kernel mZ^2/(3 pi s (s + mZ^2)) and the same-subtraction lock to a0(P) "
                "(SOURCE_SPECTRAL_THEOREM.md; the kernel identity "
                "Delta_had = 4 pi Pihat(mZ^2) is exact, lattice_backend/vector_correlator.py)"
            ),
            "f_resonances_and_continuum": (
                "unstable resonances and the continuum enter through the declared "
                "Luscher-class finite-volume pushforward from levels and residues to "
                "rho_Q(s;P); the scattering scaffolding is typed under "
                "particles/hadron/rho_scattering/ and its production execution "
                "belongs to the production backend export lane (dependency_note)"
            ),
            "g_covariance_and_uncertainty": (
                "spectral moments are linear functionals of the measure, so errors "
                "propagate through the kernel by leave-one-out jackknife over "
                "configurations; production payloads must carry the seven-row budget "
                "ledger (statistical, continuum, finite-volume, chiral, current-matching, "
                "quadrature, endpoint-remainder) enforced by the acceptance gate"
            ),
        },
    }


def build_derivation_chain() -> list[dict[str, Any]]:
    return [
        {
            "step": 1,
            "premise": "OPH QCD/matter sector and quotient ensemble",
            "source_artifact": "physics-problems/hadronic_precision_endpoint.md",
            "receipts": "particles/runs/qcd/hadron_source_backend/qcd_ensemble/",
            "conclusion": (
                "the source law mu_r^QCD is a declared positive Gibbs measure on the "
                "gauge quotient at zero density and theta = 0, with refinement maps and "
                "no-target-leak controls; a normal form is not a probability law "
                "(source-law non-substitution theorem)"
            ),
        },
        {
            "step": 2,
            "premise": "Euclidean vacuum, reflection positivity, transfer operator",
            "source_artifact": "physics-problems/hadronic_precision_endpoint.md",
            "receipts": "particles/runs/qcd/hadron_source_backend/vacuum/",
            "conclusion": (
                "the Euclidean slab with reflection positivity defines a transfer operator "
                "whose spectral decomposition supplies the Hilbert-space representation "
                "used for positivity; reflection positivity is a theorem for the "
                "unimproved Wilson action (Osterwalder-Seiler), while the "
                "clover-improved engine restores it only in the continuum limit, so at "
                "finite spacing the demonstrator checks positivity as an executed "
                "witness rather than assuming it"
            ),
        },
        {
            "step": 3,
            "premise": "renormalized electromagnetic current",
            "source_artifact": "code/particles/hadron/lattice_backend/conserved_vector.py",
            "conclusion": (
                "the point-split conserved U(1)_Q current V_mu(x) = kappa [psibar(x+mu) "
                "(1+gamma_mu) U_mu(x)^dag psi(x) - psibar(x)(1-gamma_mu) U_mu(x) psi(x+mu)] "
                "is gauge invariant and exactly conserved; the local current carries "
                "Z_V^eff = C_CL/(2 kappa C_LL) with the free-field anchor exactly 1"
            ),
        },
        {
            "step": 4,
            "premise": "Ward projection",
            "source_artifact": "code/P_derivation/SOURCE_SPECTRAL_THEOREM.md",
            "conclusion": (
                "current conservation gives q_mu Pi_Q^{mu nu} = 0, so the two-current "
                "correlator has the transverse scalar form "
                "Pi_Q^{mu nu} = (q^mu q^nu - q^2 eta^{mu nu}) Pi_Q(q^2)"
            ),
        },
        {
            "step": 5,
            "premise": "positive spectral measure with thresholds",
            "source_artifact": "physics-problems/hadronic_precision_endpoint.md",
            "conclusion": (
                "the smeared-current quadratic spectral measure of the invariant-mass "
                "operator is positive and locally finite; in finite volume it is the "
                "atomic measure rho_Q(s) = sum_n Z_n delta(s - E_n^2), Z_n >= 0, with "
                "support above the vacuum starting at the lowest vector-channel level"
            ),
        },
        {
            "step": 6,
            "premise": "resonances and continuum",
            "source_artifact": "code/particles/hadron/rho_scattering/",
            "conclusion": (
                "two-body channels map to infinite volume through the declared "
                "Luscher-class pushforward; unstable resonances and thresholds are "
                "treated explicitly at that stage (production backend export lane)"
            ),
        },
        {
            "step": 7,
            "premise": "scheme lock and transport kernel",
            "source_artifact": "code/particles/hadron/lattice_backend/vector_correlator.py",
            "conclusion": (
                "the once-subtracted spacelike kernel identity "
                "Delta_had = 4 pi Pihat(mZ^2) is exact; the measure itself is "
                "scheme-independent and the same-subtraction lock ties the transport "
                "to a0(P) in the declared convention"
            ),
        },
        {
            "step": 8,
            "premise": "covariance and uncertainty propagation",
            "source_artifact": "code/particles/hadron/verify_issue_317_spectral_measure_packet.py",
            "conclusion": (
                "leave-one-out jackknife through the full constructor gives level, "
                "residue, and correlator covariance (executed witness 4, with a "
                "positive-semidefiniteness check); production payloads carry the "
                "seven-row budget ledger"
            ),
        },
        {
            "step": 9,
            "premise": "machine-readable emission with fail-closed guards",
            "source_artifact": "code/particles/hadron/ward_projected_spectral_measure.schema.json",
            "conclusion": (
                "the export schema plus semantic gate accepts exactly the "
                "schema-conformant positive payloads and rejects measured-HVP guard "
                "flips, surrogate flags, missing budgets, quenched branches, negative "
                "residues, empty level support, and forbidden-target leaks "
                "(executed witness 5)"
            ),
        },
        {
            "step": 10,
            "premise": "blinded comparison and higher-point typing",
            "source_artifact": "code/particles/runs/qcd/hadron_source_backend/",
            "conclusion": (
                "empirical compilations are typed oph_plus_empirical_hadron_closure with "
                "promotable_as_oph_source_theorem = false and enter only through the "
                "comparison manifest; the two-current measure is a marginal of the "
                "hadronic precision functor, with four-current and transition objects "
                "separately typed (Q4_HLbL_receipt: TWO_POINT_MEASURE_INSUFFICIENT) and "
                "QED/EW radiative corrections confined to the Xi_Q remainder ledger"
            ),
        },
    ]


def build_packet() -> dict[str, Any]:
    witnesses = gauge_and_ward_witness()
    zv_witness = zv_free_field_witness()
    demonstrator = demonstrator_measure()
    gate = run_acceptance_gate()
    empirical_typing = empirical_typing_check(gate)
    higher_point_typing = higher_point_typing_check()
    closure = issue_closure_condition()

    criteria = {
        "current_and_measure_gauge_invariant_and_source_derived": {
            "packet_level_passed": bool(
                witnesses["gauge_invariance"]["passed"]
                and witnesses["ward_identity"]["passed"]
            ),
            "machine_checks": "executed witnesses 1 and 2 (machine_witnesses.gauge_and_ward)",
            "production_open_items": [
                "source QCD parameter map P -> (g_3, quark masses, scheme)",
                "qcd_ensemble/base_measure.json: REQUIRED_NOT_POPULATED",
                "currents/ward_current_definition.json: MISSING_SOURCE_CERTIFICATE",
                "the demonstrator runs at declared diagnostic bare parameters on a "
                "quenched ensemble, not on the derived physical branch",
            ],
        },
        "positivity_normalization_thresholds_scheme_dependence_proved": {
            "packet_level_passed": bool(
                zv_witness["passed"] and demonstrator["checks_passed"]
            ),
            "machine_checks": (
                "executed witnesses 3 and 4 (machine_witnesses.zv_free_field_anchor, "
                "machine_witnesses.demonstrator_measure)"
            ),
            "production_open_items": [
                "reflection-positivity and transfer-operator source certificates",
                "production thresholds and resonance pushforward execution",
                "continuum limit; five budget rows not quantified at diagnostic scale",
            ],
        },
        "machine_readable_measure_emitted_without_measured_hvp_input": {
            "packet_level_passed": bool(
                demonstrator["checks_passed"]
                and demonstrator["ensemble"]["external_inputs_used"] is False
                and gate["passed"]
            ),
            "machine_checks": (
                "executed witnesses 4 and 5 (machine_witnesses.demonstrator_measure, "
                "machine_witnesses.acceptance_gate)"
            ),
            "production_open_items": [
                "the production 2+1 payload; the diagnostic demonstrator is rejected "
                "by the production gate by construction",
            ],
        },
        "empirical_data_typed_blinded_comparison_only": {
            "packet_level_passed": bool(empirical_typing["passed"]),
            "machine_checks": (
                "live artifact check with gate-derived fail-closed witnesses "
                "(machine_witnesses.empirical_typing)"
            ),
            "production_open_items": [
                "the comparison event itself: no empirical dataset is attached, so "
                "the blinded comparison is executed only after the production "
                "payload exists, appending to the comparison manifest",
            ],
        },
        "higher_point_and_radiative_corrections_separately_typed": {
            "packet_level_passed": bool(higher_point_typing["passed"]),
            "machine_checks": "live artifact check (machine_witnesses.higher_point_typing)",
            "production_open_items": [
                "endpoint/Xi_same_scheme.json source certificate "
                "(MISSING_SOURCE_CERTIFICATE): the radiative rows are separately "
                "typed but their source content is production-scope",
            ],
        },
    }
    accepted = all(row["packet_level_passed"] for row in criteria.values())

    return {
        "issue": 317,
        "certificate_id": "issue-317-hadronic-spectral-measure-proof-packet-v1",
        "artifact": "oph_ward_projected_spectral_measure_proof_packet",
        "status": "proof_packet_construction_certified_issue_closure_condition_open",
        "accepted": bool(accepted),
        "acceptance_scope": (
            "'accepted' certifies the packet construction/contract layer only: the "
            "theorem parts, executed witnesses, typing checks, and the fail-closed "
            "emission gate. It does not assert the author-stated issue-closure "
            "condition, which is computed separately in issue_closure_condition "
            "and is currently not met locally."
        ),
        "theorem": build_theorem(),
        "derivation_chain": build_derivation_chain(),
        "machine_witnesses": {
            "gauge_and_ward": witnesses,
            "zv_free_field_anchor": zv_witness,
            "demonstrator_measure": demonstrator,
            "acceptance_gate": gate,
            "empirical_typing": empirical_typing,
            "higher_point_typing": higher_point_typing,
        },
        "acceptance_criteria_status": criteria,
        "issue_closure_condition": closure,
        "acceptance_criteria_mapping": {
            "current_and_measure_gauge_invariant_and_source_derived": (
                "theorem parts (a), (b); executed witnesses 1 and 2; source-law "
                "non-substitution theorem excludes non-source substitutes"
            ),
            "positivity_normalization_thresholds_scheme_dependence_proved": (
                "theorem parts (c), (d), (e); executed witnesses 3 and 4"
            ),
            "machine_readable_measure_emitted_without_measured_hvp_input": (
                "export schema plus semantic gate; executed witnesses 4 and 5; the "
                "demonstrator consumes bare lattice parameters only"
            ),
            "empirical_data_typed_blinded_comparison_only": (
                "machine-checked live in three layers (machine_witnesses."
                "empirical_typing): blinding of the source construction (comparison "
                "manifest empty, production-schema guard constants, gate fails closed "
                "on the measured-HVP controls), comparison-channel typing (empirical "
                "companion row class oph_plus_empirical_hadron_closure pinned "
                "non-source by schema constants), and explicit disclosure that public "
                "use of the empirical surface is confined to that non-promotable row "
                "class; no comparison event has occurred, so blinding is intact"
            ),
            "higher_point_and_radiative_corrections_separately_typed": (
                "machine-checked live for both halves (machine_witnesses."
                "higher_point_typing): higher-point half via "
                "higher_point/Q4_HLbL_receipt.json (TWO_POINT_MEASURE_INSUFFICIENT, "
                "non-promoting, no external targets); radiative half via "
                "endpoint/Xi_same_scheme.json (distinct typed ledger whose required "
                "terms include QED and EW, non-promoting) plus the production schema "
                "branch.qed exclusion; the Xi source certificate itself remains a "
                "disclosed production-level gap"
            ),
        },
        "claim_boundary": {
            "closed_here": (
                "the construction theorem for the renormalized U(1)_Q current and its "
                "positive spectral measure (gauge invariance, Ward identity, positivity, "
                "normalization, thresholds, scheme confinement, covariance and "
                "uncertainty-propagation rules), the machine-readable emission contract "
                "with fail-closed guards and negative controls, an executed diagnostic "
                "demonstrator on committed real lattice data, and the typing of "
                "empirical, higher-point, and radiative objects"
            ),
            "not_closed_here": [
                "the production 2+1 payload with real correlator arrays and quantified "
                "continuum/finite-volume/chiral/current-matching/endpoint budgets "
                "(production backend execution, dependency_note)",
                "production-scale resonance extraction through the Luscher-class map",
                "continuum-limit and confinement certificates (H4 tier)",
                "the source QCD parameter map P -> (g_3, quark masses, scheme)",
                "physical promotion of any hadronic number",
            ],
            "scope": (
                "This packet is the proof and contract layer of the hadronic spectral "
                "measure: it fixes what the production export must satisfy and proves "
                "the structural claims on the declared construction. It emits no "
                "physical hadronic prediction; the demonstrator measure is diagnostic, "
                "non-promoting, and is itself rejected by the production acceptance "
                "gate it demonstrates."
            ),
        },
        "dependency_note": {
            "direction": (
                "this proof packet is upstream of production execution: issue #425's "
                "acceptance list requires that the production backend 'exports the "
                "Ward-current spectral payload consumed by #317', i.e. #425 supplies a "
                "payload conforming to the contract certified here; no result of #425 "
                "is consumed by this packet"
            ),
            "upstream_of_425_not_downstream": True,
            "gauge_carrier_selection": (
                "the physical gauge/matter carrier selection (#590) gates #425's "
                "physical promotion, not this packet"
            ),
        },
        "dependency_artifacts": {
            "export_schema": "code/particles/hadron/ward_projected_spectral_measure.schema.json",
            "constructive_contract": "code/particles/runs/hadron/ward_projected_spectral_measure_contract.json",
            "nonidentifiability_obstruction": "code/particles/hadron/derive_ward_projected_spectral_measure_obstruction.py",
            "endpoint_reduction_theorem": "code/P_derivation/source_spectral_theorem.py (WardProjectedHadronicSpectralEmission_Q)",
            "backend_spec_note": "physics-problems/hadronic_precision_endpoint.md",
            "receipt_bundle": "code/particles/runs/qcd/hadron_source_backend/ (claim: SOURCE_PROTOTYPE_NOT_PROMOTED)",
            "lattice_engine": "code/particles/hadron/lattice_backend/",
            "demonstrator_ensemble": "code/particles/runs/hadron/hybrid_ir_ensembleA_2026-07-16.npz",
        },
        "consumer_artifacts": {
            "production_backend_export": (
                "production backend export bundle conforming to the schema and "
                "acceptance gate certified here (dependency_note)"
            ),
            "thomson_endpoint_lane": (
                "code/P_derivation/source_spectral_theorem.py consumes a conforming "
                "payload for the source-only Thomson endpoint interval certificate"
            ),
        },
        "forbidden_inputs": FORBIDDEN_TARGETS,
        "verifier_command": (
            "python3 code/particles/hadron/verify_issue_317_spectral_measure_packet.py "
            "--check --output "
            "code/particles/runs/hadron/ward_projected_spectral_measure_proof_packet.json"
        ),
        "volatile": {
            "generated_utc": _now_utc(),
            "note": "timestamp only; excluded from acceptance logic",
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify the issue #317 spectral-measure proof packet.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    parser.add_argument("--check", action="store_true", help="exit nonzero unless the packet is accepted")
    args = parser.parse_args()

    packet = build_packet()
    text = json.dumps(packet, indent=1) + "\n"
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(text, encoding="utf-8")

    witnesses = packet["machine_witnesses"]
    print(f"gauge invariance defect: {witnesses['gauge_and_ward']['gauge_invariance']['correlator_relative_defect']:.2e}")
    print(f"ward identity relative defect: {witnesses['gauge_and_ward']['ward_identity']['relative_defect']:.2e}")
    print(f"Z_V free-field max deviation: {witnesses['zv_free_field_anchor']['max_abs_deviation_from_one']:.2e}")
    demo = witnesses["demonstrator_measure"]["extraction"]
    print(f"demonstrator: a*m_V = {demo['am_vector']:.4f} +- {demo['am_vector_jackknife_error']:.4f}, "
          f"Z_V = {demo['z_v']:.4f}, weight = {demo['weight_physical_units']:.4f}")
    gate = witnesses["acceptance_gate"]
    print(f"acceptance gate: conformant accepted = {gate['conformant_payload_accepted']}, "
          f"negative controls rejected = {gate['all_negative_controls_rejected']} "
          f"({len(gate['negative_controls'])} controls)")
    print(f"empirical typing check passed: {witnesses['empirical_typing']['passed']}")
    print(f"higher-point typing check passed: {witnesses['higher_point_typing']['passed']}")
    print(f"packet accepted: {packet['accepted']}")
    print(f"issue closure condition met locally: {packet['issue_closure_condition']['met_locally']}")
    print(f"wrote {out_path}")
    if args.check and not packet["accepted"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
