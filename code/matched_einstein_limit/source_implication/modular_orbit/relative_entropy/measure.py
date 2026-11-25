"""Compute the protected word and stationary angular relative-entropy records."""

import argparse
import json
from pathlib import Path

from . import angular, written


def run():
    return dict(written=written.receipt(), angular=angular.receipt())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    out = json.dumps(run(), allow_nan=False, indent=2) + "\n"
    if args.write:
        Path(__file__).with_name("results.json").write_text(out)
    else:
        print(out, end="")
