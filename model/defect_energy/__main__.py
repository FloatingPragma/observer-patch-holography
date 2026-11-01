"""Replay the exact same-defect repair/action receipt."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from .energy import report


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).with_name("results.json")
# Direct calculation inputs, and the checked implementation behind the raw
# repair/memory interface.  The wider source closure belongs to its own packet.
SOURCES = (
    "model/defect_energy/__init__.py", "model/defect_energy/__main__.py",
    "model/defect_energy/energy.py", "model/defect_energy/test_energy.py",
    "model/defect_energy/README.md", "model/defect_protection/protection.py",
    "model/defect_protection/service.py", "model/exchange_sector/statistics.py",
    "model/exchange_sector/statistics_service.py", "model/observer/native.py",
    "model/observer/record_recovery.py", "model/common_source/record_code.py",
    "model/common_source/native_routing.py", "model/observer/source_scalar/compiler.py",
    "Lean/Geometry/InternalEnergyInertia.lean", "Lean/Geometry/ProperTimeInternalAction.lean",
)


def build():
    return {**report(), "source_sha256": {
        name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in SOURCES}}


def check(actual, expected):
    if json.dumps(actual, sort_keys=True) != json.dumps(expected, sort_keys=True):
        raise ValueError("same-defect repair/action receipt differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--write", nargs="?", const=OUTPUT, type=Path)
    mode.add_argument("--check", nargs="?", const=OUTPUT, type=Path)
    args = parser.parse_args()
    result = build()
    if args.write:
        args.write.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    elif args.check:
        check(strict_json.load(args.check), result)
    else:
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
