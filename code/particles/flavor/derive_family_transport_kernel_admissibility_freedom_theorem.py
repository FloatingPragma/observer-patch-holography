#!/usr/bin/env python3
"""Kernel admissibility freedom theorem (#377): certificates fix no spectrum.

Theorem (kernel admissibility freedom). Let the admissibility battery of
the family-transport lane consist of the certificates the corpus records
for the on-disk kernel: two refinement levels with positive-semidefinite
hermitian descendants, an open three-cluster gap at every level, a simple
centered spectrum, a conjugacy-Riesz margin (defect supremum below half
the minimal gap), persistent projector labeling (every same-label overlap
amplitude exceeds every cross-label amplitude), and nondegenerate
overlap-edge amplitudes above the recorded floor. Then for every target
pair (r, s) with r > 0 and s > 0 there exists a kernel family passing the
full battery whose centered compressed branch generator has raw gap ratio
exactly r and spectral span exactly s. Consequently the battery constrains
the emitted invariants

    rho_ord = 3/(2 + r),   x2 = (r - 1)/(r + 1),   span = s

not at all: the map from certificate-passing kernels onto (r, s) is
surjective onto (0, infinity)^2.

Proof (constructive, executed below). Fix the trace-free target spectrum
lambda with gaps g21 = s r/(1+r), g32 = s/(1+r), shift it positive, and
set T_0 = Q_0 diag(sqrt(mu)) Q_0^dagger with Q_0 = exp(i eps_0 H) for a
fixed fully off-diagonal Hermitian mixer H. Then the hermitian descendant
T_0 T_0^dagger = Q_0 diag(mu) Q_0^dagger has exactly the prescribed
spectrum, since unitary conjugation preserves eigenvalues. Level one uses
a slightly rotated frame and a relative eigenvalue drift matched to the
on-disk kernel scale, shrunk geometrically until the Riesz margin
certificate passes, which terminates because defects scale linearly in
the shrink factor while the gap is fixed. Off-diagonal mixer entries make
every cross-label eigenline overlap nonzero, clearing the amplitude
floor, while the small frame difference keeps same-label overlaps
dominant. Every certificate is then checked numerically for a witness
grid that brackets the on-disk operating point by more than an order of
magnitude in both coordinates, plus the operating point itself.

Corollary (content of K1). The persistence certificates are a filter, not
a generator: deriving a kernel that merely passes them cannot emit
rho_ord, x2, or the spans. The load-bearing content of the issue-377
program is a selection principle for the transport operator itself, which
on the current corpus is exactly the statement of the axiom-level
nonidentifiability theorem transported to the operator side.

No quark reference value, fitted spread, or flavor template enters any
step; the construction is target-free by inspection and the witness grid
covers arbitrary (r, s), not preferred values.

Run:
    python3 code/particles/flavor/derive_family_transport_kernel_admissibility_freedom_theorem.py
writes code/particles/runs/flavor/family_transport_kernel_admissibility_freedom_theorem.json.
"""

from __future__ import annotations

import argparse
import json
import math
import pathlib
from datetime import datetime, timezone

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RUNS = ROOT / "particles" / "runs" / "flavor"
DEFAULT_OUT = RUNS / "family_transport_kernel_admissibility_freedom_theorem.json"
KERNEL_PATH = RUNS / "family_transport_kernel.json"

AMPLITUDE_FLOOR = 1.0e-12
MIXER = np.asarray([[0.0, 1.0, 1.0j],
                    [1.0, 0.0, 1.0],
                    [-1.0j, 1.0, 0.0]], dtype=complex)
EPS_FRAME_0 = 0.05
EPS_FRAME_DRIFT = 0.02
EIGENVALUE_DRIFT_REL = 0.003

R_GRID = (0.05, 0.317889, 1.0, 5.0)      # includes the on-disk raw ratio
S_GRID = (0.1, 1.253553, 10.0)           # includes the on-disk span


def _timestamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _expm_hermitian_i(h: np.ndarray, eps: float) -> np.ndarray:
    evals, evecs = np.linalg.eigh(h)
    return evecs @ np.diag(np.exp(1j * eps * evals)) @ evecs.conj().T


def _centered_spectrum(r: float, s: float) -> np.ndarray:
    g21 = s * r / (1.0 + r)
    g32 = s / (1.0 + r)
    lam1 = -(2.0 * g21 + g32) / 3.0
    return np.asarray([lam1, lam1 + g21, lam1 + g21 + g32])


def _eigenlines(hermitian: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    evals, evecs = np.linalg.eigh(hermitian)
    return evals, evecs


def build_witness(r: float, s: float) -> dict:
    lam = _centered_spectrum(r, s)
    shift = 1.0 + abs(float(lam[0]))
    mu = lam + shift
    gap_min = float(min(lam[1] - lam[0], lam[2] - lam[1]))

    q0 = _expm_hermitian_i(MIXER, EPS_FRAME_0)
    t0 = q0 @ np.diag(np.sqrt(mu)) @ q0.conj().T

    shrink, attempts = 1.0, 0
    while True:
        attempts += 1
        eps1 = EPS_FRAME_0 + EPS_FRAME_DRIFT * shrink
        drift = EIGENVALUE_DRIFT_REL * shrink
        q1 = _expm_hermitian_i(MIXER, eps1)
        mu1 = mu * (1.0 + drift)
        t1 = q1 @ np.diag(np.sqrt(mu1)) @ q1.conj().T
        defect = float(np.linalg.norm(t1 - t0, ord=2))
        if defect < 0.5 * gap_min or attempts > 60:
            break
        shrink *= 0.5

    certificates = {}
    descendants = []
    projector_frames = []
    for level, t in enumerate((t0, t1)):
        desc = t @ t.conj().T
        evals, evecs = _eigenlines(desc)
        centered = desc - (np.trace(desc) / 3.0) * np.eye(3, dtype=complex)
        cevals = np.linalg.eigvalsh(centered)
        g21 = float(cevals[1] - cevals[0])
        g32 = float(cevals[2] - cevals[1])
        descendants.append({
            "level": level,
            "psd": bool(np.min(np.linalg.eigvalsh(desc)) > 0.0),
            "centered_eigenvalues": [float(x) for x in cevals],
            "g21": g21,
            "g32": g32,
        })
        projector_frames.append(evecs)
    certificates["psd_descendants"] = all(d["psd"] for d in descendants)
    certificates["three_cluster_gap_open_all_levels"] = all(
        min(d["g21"], d["g32"]) > 0.0 for d in descendants)
    certificates["simple_centered_spectrum"] = all(
        min(d["g21"], d["g32"]) > 1.0e-12 for d in descendants)
    certificates["riesz_margin_defect_below_half_gap"] = bool(
        defect < 0.5 * gap_min)

    left, right = projector_frames
    overlaps = np.abs(left.conj().T @ right) ** 2
    same_label = [float(overlaps[i, i]) for i in range(3)]
    cross_label = [float(overlaps[i, j]) for i in range(3) for j in range(3)
                   if i != j]
    certificates["projector_labeling_persistent"] = bool(
        min(same_label) > max(cross_label))
    certificates["edge_amplitudes_above_floor"] = bool(
        min(cross_label) > AMPLITUDE_FLOOR)

    base = descendants[0]
    r_emitted = base["g21"] / base["g32"]
    s_emitted = base["g21"] + base["g32"]
    rho_emitted = 3.0 * base["g32"] / (2.0 * base["g32"] + base["g21"])
    x2_emitted = (r_emitted - 1.0) / (r_emitted + 1.0)

    return {
        "target": {"r": r, "s": s},
        "emitted": {
            "r": r_emitted,
            "s": s_emitted,
            "rho_ord": rho_emitted,
            "x2": x2_emitted,
        },
        "exactness": {
            "r_abs_error": abs(r_emitted - r),
            "s_abs_error": abs(s_emitted - s),
        },
        "riesz": {"defect_sup": defect, "half_gap": 0.5 * gap_min,
                  "shrink_attempts": attempts},
        "certificates": certificates,
        "all_certificates_pass": all(bool(v) for v in certificates.values()),
    }


def build() -> dict:
    kernel = json.loads(KERNEL_PATH.read_text(encoding="utf-8"))

    witnesses = []
    for r in R_GRID:
        for s in S_GRID:
            witnesses.append(build_witness(r, s))
    all_pass = all(w["all_certificates_pass"] for w in witnesses)
    max_r_err = max(w["exactness"]["r_abs_error"] for w in witnesses)
    max_s_err = max(w["exactness"]["s_abs_error"] for w in witnesses)
    if not all_pass:
        raise AssertionError("a freedom witness failed the certificate battery")
    if max(max_r_err, max_s_err) > 1.0e-9:
        raise AssertionError("a freedom witness missed its target invariants")

    span_r = (min(R_GRID), max(R_GRID))
    span_s = (min(S_GRID), max(S_GRID))

    return {
        "artifact": "oph_family_transport_kernel_admissibility_freedom_theorem",
        "generated_utc": _timestamp(),
        "github_issues": [377, 379, 380],
        "proof_status": "closed_constructive_freedom_theorem",
        "claim_tier": "certificate_battery_nonselection_obstruction",
        "row_class": "theorem_grade_obstruction",
        "guards": {
            "quark_reference_values_consumed": False,
            "fitted_spreads_consumed": False,
            "numerical_flavor_template_consumed": False,
            "public_promotion_allowed": False,
        },
        "theorem_statement": (
            "For every (r, s) in (R_>0)^2 there is a two-level transport "
            "kernel passing the full recorded admissibility battery (PSD "
            "hermitian descendants, open three-cluster gaps, simple "
            "centered spectrum, conjugacy-Riesz margin, persistent "
            "projector labeling, overlap-edge amplitudes above the floor) "
            "whose centered compressed branch generator has raw gap ratio "
            "exactly r and spectral span exactly s. The battery therefore "
            "places no constraint on (rho_ord, x2, span): certificate-"
            "passing kernels surject onto the full invariant space."
        ),
        "proof_kind": "explicit_construction_with_executed_certificates",
        "construction": {
            "spectrum_placement": "trace-free target spectrum with gaps "
                                  "(s r/(1+r), s/(1+r)), shifted positive; "
                                  "T = Q diag(sqrt(mu)) Q_dagger so the "
                                  "descendant spectrum is exact under "
                                  "unitary conjugation",
            "frames": "Q_level = exp(i eps_level H) with a fully "
                      "off-diagonal Hermitian mixer H; the frame "
                      "difference makes every cross-label overlap nonzero "
                      "while keeping same-label overlaps dominant",
            "riesz_termination": "frame drift and eigenvalue drift shrink "
                                 "geometrically until the defect supremum "
                                 "is below half the minimal gap; "
                                 "termination is guaranteed because the "
                                 "defect is linear in the shrink factor "
                                 "at fixed gap",
        },
        "witness_grid": {
            "r_values": list(R_GRID),
            "s_values": list(S_GRID),
            "includes_on_disk_operating_point": True,
            "r_span_orders_of_magnitude": math.log10(span_r[1] / span_r[0]),
            "s_span_orders_of_magnitude": math.log10(span_s[1] / span_s[0]),
            "all_certificates_pass": all_pass,
            "max_target_error": max(max_r_err, max_s_err),
        },
        "witnesses": witnesses,
        "identity_x2_of_r": {
            "statement": "x2 = (r - 1)/(r + 1) holds identically, so the "
                         "mean-law coordinate carries no information "
                         "beyond the gap ratio",
            "max_deviation_on_witnesses": max(
                abs(w["emitted"]["x2"]
                    - (w["emitted"]["r"] - 1.0) / (w["emitted"]["r"] + 1.0))
                for w in witnesses),
        },
        "on_disk_kernel_context": {
            "artifact": kernel.get("artifact"),
            "status": kernel.get("status"),
            "note": "the on-disk template is one point of this free family; "
                    "see the rho_ord sensitivity audit for its local "
                    "elasticity",
        },
        "corollary_for_issue_377": (
            "Deriving a kernel that passes the persistence certificates "
            "cannot close K1: the certificates are a filter with zero "
            "selective power over the emitted invariants. The remaining "
            "burden of the kernel program is a source selection principle "
            "for the transport operator itself, equivalently a derivation "
            "of the invariant triple named by the three-scalar interface "
            "theorem."
        ),
        "relation_to_nonidentifiability": (
            "This is the operator-side transport of the axiom-level "
            "nonidentifiability theorem: the free (R_>0)^2 spread fiber "
            "reappears as the free (r, s) invariant space of "
            "certificate-passing kernels."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Prove the kernel admissibility freedom theorem "
                    "constructively.")
    parser.add_argument("--output", default=str(DEFAULT_OUT))
    args = parser.parse_args()

    report = build()
    out_path = pathlib.Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")

    grid = report["witness_grid"]
    print(f"witness grid: r in {grid['r_values']}  s in {grid['s_values']}")
    print(f"all certificates pass: {grid['all_certificates_pass']}")
    print(f"max target error: {grid['max_target_error']:.2e}")
    print(f"x2 = (r-1)/(r+1) max deviation: "
          f"{report['identity_x2_of_r']['max_deviation_on_witnesses']:.2e}")
    for w in report["witnesses"]:
        t, e = w["target"], w["emitted"]
        print(f"  r={t['r']:<9g} s={t['s']:<9g} -> rho={e['rho_ord']:.6f} "
              f"x2={e['x2']:+.6f}  pass={w['all_certificates_pass']}")
    print(f"proof status: {report['proof_status']}")
    print(f"saved: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
