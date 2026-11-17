"""Replay the finite charged-energy apparatus and every source pin."""
import argparse
import json
from pathlib import Path

from .certificate import certificate, validate_pins


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--write", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    path = Path(__file__).with_name("results.json")
    if args.check:
        saved = json.loads(path.read_text())
        validate_pins(saved["input_sha256"])
    result = certificate()
    if args.write:
        path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    elif result != saved:
        raise ValueError("complete finite charged-energy receipt differs")
    print("PASS: full-source charged phase readout, all outcomes and frozen science comparison")


if __name__ == "__main__":
    main()
