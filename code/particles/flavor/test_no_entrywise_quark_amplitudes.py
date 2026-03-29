#!/usr/bin/env python3
"""Reject dense entrywise amplitude imports in certified quark mode."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "forward_yukawas.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the no-entrywise-quark-amplitudes guard.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input forward Yukawa artifact.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    if bool(payload.get("forward_certified", False)) and bool(payload.get("dense_entrywise_amplitude_used", False)):
        print("forward-certified quark artifact still uses dense entrywise amplitudes", file=sys.stderr)
        return 1
    print("no-entrywise-quark-amplitudes guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
