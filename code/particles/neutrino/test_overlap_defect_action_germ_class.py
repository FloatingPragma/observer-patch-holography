#!/usr/bin/env python3
"""Validate the local quadratic Majorana action-germ class."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys

import numpy as np


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "neutrino" / "majorana_overlap_defect_action_germ.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the Majorana overlap-defect action germ.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT))
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    if str(payload.get("proof_status", "")) != "local_quadratic_germ_closed":
        print("action germ is not marked as locally closed", file=sys.stderr)
        return 1
    template = np.asarray(payload.get("hessian_class_residual_2x2"), dtype=float)
    target = np.asarray([[2.0, 1.0], [1.0, 2.0]], dtype=float)
    if template.shape != (2, 2) or not np.allclose(template, target, atol=1.0e-12, rtol=1.0e-12):
        print("action germ does not carry the expected residual Hessian class", file=sys.stderr)
        return 1
    print("Majorana action-germ class guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
