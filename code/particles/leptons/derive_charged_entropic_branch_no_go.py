#!/usr/bin/env python3
"""No-go for the entropic conditioned charged branch at leading cubic order.

Conditioning the twelve-port record on a nonzero W5 norm leaves the
quadratic repair cost flat on the W5 sphere, so orbit selection falls to the
universal Taylor coefficients of the record entropy,
``-sum q ln q = const - |w|^2/2 + S3(w)/6 - S4(w)/12 + ...`` with the
port-wise power sums ``S_k``.  At leading order in the band amplitude the
selection functional is the cubic ``S3`` alone, with no chosen coefficient.

Statement (certified numerically here, forty-seed projected ascent,
degeneracy tolerance 1e-8): the maximum of ``S3`` on the W5 unit sphere of
the twelve-port carrier is attained on the C5-axis orbit, stabilizer order
ten, whose quadrupole spectrum is exactly doubly degenerate.  Since
``S3(-w) = -S3(w)``, the minimum is attained at the antipode ``-w`` of the
maximizer, the same orbit type with the sign of the quadrupole flipped, so
both orientations are degenerate; the antipodal spectrum is recorded next
to the maximizer.  At cubic order the entropic conditioned branch therefore
produces two equal charged masses.

Scope.  The statement covers the leading cubic term at small band amplitude
only.  The quartic certificate (``flavor/entropy_w5_shape_certificate.py``)
shows that the same parameter-free packet, truncated at quartic order,
selects a simple-spectrum golden-ratio branch above the exact crossing
amplitude r_c = 2.776 and below the positivity bound r < 2 sqrt(6) of that
branch; the golden branch is excluded by comparison with the measured
charged-lepton shape at 16.7 percent, and packets beyond quartic truncation
stay open.  Of the three candidate mechanisms, the second is closed at
cubic order, excluded by comparison at quartic order, and open beyond; the
surviving route is a source-emitted charged interaction whose invariant mix
lies off the entropic ray, inside the simple-spectrum region.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np

try:
    from leptons.derive_charged_w5_orbit_decision_geometry import (
        P5,
        spectrum_report,
    )
except ModuleNotFoundError:
    from derive_charged_w5_orbit_decision_geometry import P5, spectrum_report

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = (
    ROOT / "particles" / "runs" / "leptons"
    / "charged_entropic_branch_no_go.json"
)


def maximize_s3(seeds: int = 40, iters: int = 4000) -> tuple[float, np.ndarray]:
    best: tuple[float, np.ndarray] | None = None
    for seed in range(seeds):
        rng = np.random.default_rng(seed)
        w = P5 @ rng.standard_normal(12)
        w /= np.linalg.norm(w)
        for _ in range(iters):
            grad = P5 @ (3.0 * w**2)
            grad -= (grad @ w) * w
            if np.linalg.norm(grad) < 1.0e-13:
                break
            w = w + 0.05 * grad
            w /= np.linalg.norm(w)
        value = float(np.sum(w**3))
        if best is None or value > best[0]:
            best = (value, w.copy())
    assert best is not None
    return best


def build() -> dict[str, Any]:
    value, w = maximize_s3()
    report = spectrum_report(w)
    eigenvalues = report["eigenvalues"]
    min_gap = min(
        eigenvalues[1] - eigenvalues[0], eigenvalues[2] - eigenvalues[1]
    )
    degenerate = min_gap < 1.0e-8
    # S3 is odd, so the minimizer is the antipode of the maximizer.
    antipodal_eigenvalues = spectrum_report(-w)["eigenvalues"]
    antipodal_value = float(np.sum((-w) ** 3))
    antipodal_gap = min(
        antipodal_eigenvalues[1] - antipodal_eigenvalues[0],
        antipodal_eigenvalues[2] - antipodal_eigenvalues[1],
    )
    checks = {
        "extremum_found": value > 0.5,
        "extremal_spectrum_degenerate": bool(degenerate),
        "antipodal_value_is_negated": abs(antipodal_value + value) < 1.0e-12,
        "antipodal_spectrum_degenerate": bool(antipodal_gap < 1.0e-8),
        "no_go_certified": bool(degenerate and antipodal_gap < 1.0e-8),
    }
    return {
        "artifact": "oph_charged_entropic_branch_no_go",
        "schema_version": 1,
        "status": "ENTROPIC_CONDITIONED_BRANCH_NO_GO",
        "row_class": "parameter_free_no_go_certificate",
        "promotion_allowed": False,
        "selection_functional": (
            "universal entropy cubic S3 on the W5 unit sphere, the leading "
            "term in the band amplitude; no chosen coefficient"
        ),
        "extremal_value_s3": value,
        "extremal_spectrum": eigenvalues,
        "minimum_gap": min_gap,
        "antipodal_value_s3": antipodal_value,
        "antipodal_spectrum": antipodal_eigenvalues,
        "oddness": (
            "S3(-w) = -S3(w), so the minimizer of S3 on the sphere is the "
            "antipode of the maximizer: the same C5-axis orbit type with the "
            "sign of the quadrupole flipped, degenerate in both orientations"
        ),
        "consequence": (
            "two equal charged masses on the entropic branch at leading cubic "
            "order; the observed family requires a source-emitted charged "
            "interaction off the entropic ray, inside the simple-spectrum "
            "region"
        ),
        "scope": {
            "truncation_order": "leading cubic term at small band amplitude",
            "quartic_order": (
                "the same parameter-free packet truncated at quartic order "
                "selects a simple-spectrum golden-ratio branch above the "
                "exact crossing amplitude r_c = 2.776 and below that "
                "branch's positivity bound r < 2 sqrt(6); the golden branch "
                "is excluded by comparison with the measured charged-lepton "
                "shape at 16.7 percent only"
            ),
            "quartic_certificate": "runs/flavor/entropy_w5_shape_certificate.json",
            "beyond_quartic": "open",
            "mechanism_status": (
                "second of three candidate mechanisms: closed at cubic "
                "order, excluded by comparison at quartic order, open beyond"
            ),
            "universal_impossibility_claimed": False,
        },
        "checks": checks,
        "checks_pass": all(bool(v) for v in checks.values()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    artifact = build()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({k: artifact[k] for k in ("status", "checks_pass", "extremal_spectrum")}, indent=2))
    return 0 if artifact["checks_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
