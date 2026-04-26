#!/usr/bin/env python3
"""Emit a numerical witness for the OPH P/alpha fixed-point closure."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from paper_math import build_fixed_point_witness


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Sample the OPH alpha -> alpha fixed-point map and emit a compare-only witness."
    )
    parser.add_argument(
        "--mode",
        choices=("thomson_structured_running", "thomson_structured_running_asymptotic", "mz_anchor"),
        default="thomson_structured_running",
        help="Which alpha readout to feed into P = phi + alpha*sqrt(pi).",
    )
    parser.add_argument("--precision", type=int, default=40, help="Decimal precision for the solver.")
    parser.add_argument("--su2-cutoff", type=int, default=120, help="Representation cutoff for the SU(2) edge sum.")
    parser.add_argument("--su3-cutoff", type=int, default=90, help="Representation cutoff for the SU(3) edge sum.")
    parser.add_argument("--scan-points", type=int, default=60, help="Alpha-space scan points used to bracket closure.")
    parser.add_argument("--max-iterations", type=int, default=20, help="Maximum outer fixed-point iterations.")
    parser.add_argument(
        "--derivative-step",
        default="0.000001",
        help="Alpha-space finite-difference step used for local slope diagnostics.",
    )
    parser.add_argument("--sample-points", type=int, default=5, help="Number of alpha probes around the fixed point.")
    parser.add_argument("--output", help="Optional path for the JSON witness.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    witness = build_fixed_point_witness(
        precision=args.precision,
        mode=args.mode,
        su2_cutoff=args.su2_cutoff,
        su3_cutoff=args.su3_cutoff,
        scan_points=args.scan_points,
        max_iterations=args.max_iterations,
        derivative_step=args.derivative_step,
        sample_points=args.sample_points,
    )

    text = json.dumps(witness, indent=2, sort_keys=True) + "\n"
    if args.output:
        out_path = Path(args.output)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(text, encoding="utf-8")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
