#!/usr/bin/env python3
"""CLI for the paper-math P/alpha closure experiment."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from paper_math import build_report


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Solve the OPH P/alpha closure directly from the paper equations.")
    parser.add_argument(
        "--mode",
        choices=("thomson_structured_running", "thomson_structured_running_asymptotic", "mz_anchor"),
        default="thomson_structured_running",
        help="Which alpha readout to feed into P = phi + alpha*sqrt(pi).",
    )
    parser.add_argument("--precision", type=int, default=40, help="Decimal precision for the paper-math solver.")
    parser.add_argument("--su2-cutoff", type=int, default=120, help="Representation cutoff for the SU(2) edge sum.")
    parser.add_argument("--su3-cutoff", type=int, default=90, help="Representation cutoff for the SU(3) edge sum.")
    parser.add_argument("--scan-points", type=int, default=60, help="Alpha-space scan points used to bracket closure.")
    parser.add_argument("--max-iterations", type=int, default=20, help="Maximum outer fixed-point iterations.")
    parser.add_argument("--json", action="store_true", help="Print the full report as JSON.")
    parser.add_argument("--output", help="Optional path for the JSON report.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = build_report(
        precision=args.precision,
        mode=args.mode,
        su2_cutoff=args.su2_cutoff,
        su3_cutoff=args.su3_cutoff,
        scan_points=args.scan_points,
        max_iterations=args.max_iterations,
    )

    if args.output:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    if args.json:
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0

    print(f"mode                     = {report['mode']}")
    print(f"precision                = {report['precision']}")
    print(f"alpha                    = {report['alpha']}")
    print(f"alpha^-1                 = {report['alpha_inv']}")
    print(f"P                        = {report['p']}")
    print(f"source anchor a0(P)      = {report['source_anchor_alpha_inv']}")
    print(f"alpha fixed-point resid. = {report['alpha_fixed_point_residual']}")
    print("")
    print("D10 point at the closure solution")
    print(f"  alpha_U                = {report['d10']['alpha_u']}")
    print(f"  M_U                    = {report['d10']['mu_u']}")
    print(f"  mZ_run                 = {report['d10']['mz_run']}")
    print(f"  v                      = {report['d10']['v']}")
    print(f"  alpha1(mZ)             = {report['d10']['alpha1_mz']}")
    print(f"  alpha2(mZ)             = {report['d10']['alpha2_mz']}")
    print(f"  alpha3(mZ)             = {report['d10']['alpha3_mz']}")
    if report["structured_running"] is not None:
        print("")
        print("Internal structured Thomson running")
        print(f"  mass source            = {report['structured_running']['mass_source']}")
        print(f"  transport kernel       = {report['structured_running']['transport_kernel']}")
        print(f"  lepton Delta alpha^-1  = {report['structured_running']['lepton_delta_alpha_inv']}")
        print(f"  quark naive            = {report['structured_running']['quark_delta_alpha_inv_naive']}")
        print(f"  quark screened         = {report['structured_running']['quark_delta_alpha_inv_screened']}")
        print(f"  total Delta alpha^-1   = {report['structured_running']['total_delta_alpha_inv']}")
        if "kernel_upgrade_delta_alpha_inv" in report["structured_running"]:
            print(f"  exact-vs-asymptotic    = {report['structured_running']['kernel_upgrade_delta_alpha_inv']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
