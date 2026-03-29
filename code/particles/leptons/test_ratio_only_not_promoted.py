#!/usr/bin/env python3
"""Ensure the downstream completion surface does not promote a non-absolute lepton artifact."""

from __future__ import annotations

import argparse
import json
import pathlib
import subprocess
import sys
import tempfile


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "leptons" / "forward_charged_leptons.json"
COMPLETION_SCRIPT = ROOT / "particles" / "leptons" / "run_qed_ew_completion.py"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate ratio-only promotion guards for the charged-lepton lane.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input forward charged-lepton artifact.")
    args = parser.parse_args()

    input_path = pathlib.Path(args.input)
    payload = json.loads(input_path.read_text(encoding="utf-8"))
    closure_state = str(payload.get("closure_state", payload.get("theorem_status", "open")))

    with tempfile.TemporaryDirectory(prefix="oph_lepton_ratio_guard_") as tmpdir:
        output_path = pathlib.Path(tmpdir) / "completion.json"
        result = subprocess.run(
            [sys.executable, str(COMPLETION_SCRIPT), "--input", str(input_path), "--output", str(output_path)],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            print(result.stderr.strip() or "completion command failed", file=sys.stderr)
            return 1
        completion = json.loads(output_path.read_text(encoding="utf-8"))

    blocked = bool(completion.get("promotion_blocked", False))
    reported_masses = completion.get("reported_masses_gev")

    if closure_state == "absolute_scale_closed":
        if blocked or reported_masses is None:
            print("absolute-scale-closed artifact did not complete cleanly", file=sys.stderr)
            return 1
        print("absolute-scale completion path is open")
        return 0

    if not blocked or reported_masses is not None:
        print("non-absolute artifact was promoted by the completion surface", file=sys.stderr)
        return 1

    print("ratio/open artifact remains blocked downstream")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
