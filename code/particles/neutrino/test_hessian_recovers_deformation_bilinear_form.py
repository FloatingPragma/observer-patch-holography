#!/usr/bin/env python3
"""Ensure the upstream Majorana Hessian artifact reproduces the current deformation class."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_HESSIAN = ROOT / "particles" / "runs" / "neutrino" / "majorana_overlap_defect_hessian.json"
DEFAULT_DEFORMATION = ROOT / "particles" / "runs" / "neutrino" / "majorana_deformation_bilinear_form.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate Hessian recovery of the deformation bilinear form.")
    parser.add_argument("--hessian", default=str(DEFAULT_HESSIAN))
    parser.add_argument("--deformation", default=str(DEFAULT_DEFORMATION))
    args = parser.parse_args()

    hessian = json.loads(pathlib.Path(args.hessian).read_text(encoding="utf-8"))
    deformation = json.loads(pathlib.Path(args.deformation).read_text(encoding="utf-8"))
    h_residual = np.asarray(hessian.get("residual_hessian_2x2"), dtype=float)
    d_residual = np.asarray(deformation.get("residual_metric_class_2x2"), dtype=float)
    if h_residual.shape != (2, 2) or d_residual.shape != (2, 2):
        print("missing 2x2 residual matrices", file=sys.stderr)
        return 1
    if not np.allclose(h_residual, d_residual, atol=1.0e-12, rtol=1.0e-12):
        print("overlap-defect Hessian does not recover the deformation bilinear-form class", file=sys.stderr)
        return 1
    print("Majorana Hessian recovery guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
