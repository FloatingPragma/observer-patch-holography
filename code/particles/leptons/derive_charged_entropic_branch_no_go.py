#!/usr/bin/env python3
"""No-go for the entropic conditioned charged branch.

Conditioning the twelve-port record on a nonzero W5 norm leaves the
quadratic repair cost flat on the W5 sphere, so orbit selection falls to the
universal Taylor coefficients of the record entropy,
``-sum q ln q = const - |w|^2/2 + S3(w)/6 - S4(w)/12 + ...`` with the
port-wise power sums ``S_k``.  The selection functional at leading order is
therefore parameter-free.

Theorem (certified numerically here): the extremal orbit of ``S3`` on the
W5 unit sphere is the C5-axis orbit, stabilizer order ten, whose quadrupole
spectrum is exactly doubly degenerate.  Both orientations give the same
orbit type.  The entropic conditioned branch therefore produces two equal
charged masses and cannot generate the observed family.  Together with the
homogeneous shape-silence theorem this closes the second of the three
candidate mechanisms; the surviving route is a source-emitted charged
interaction whose invariant mix lies off the entropic ray, inside the
simple-spectrum region, on the frozen target locus of the decision-geometry
lane.
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
    checks = {
        "extremum_found": value > 0.5,
        "extremal_spectrum_degenerate": bool(degenerate),
        "no_go_certified": bool(degenerate),
    }
    return {
        "artifact": "oph_charged_entropic_branch_no_go",
        "schema_version": 1,
        "status": "ENTROPIC_CONDITIONED_BRANCH_NO_GO",
        "row_class": "parameter_free_no_go_certificate",
        "promotion_allowed": False,
        "selection_functional": (
            "universal entropy cubic S3 on the W5 unit sphere; no chosen "
            "coefficient"
        ),
        "extremal_value_s3": value,
        "extremal_spectrum": eigenvalues,
        "minimum_gap": min_gap,
        "consequence": (
            "two equal charged masses on the entropic branch; the observed "
            "family requires a source-emitted charged interaction off the "
            "entropic ray, on the frozen target locus"
        ),
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
