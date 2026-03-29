#!/usr/bin/env python3
"""Ensure the neutrino lane keeps only the declared symmetric-diagonal residual."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "sector_transport_pushforward.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate neutrino residual factorization.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input sector-response artifact.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    nu = dict(payload.get("sector_response_object", {}).get("nu", {}))
    if not nu:
        print("missing neutrino sector response", file=sys.stderr)
        return 1

    if nu.get("normalization_class") != "symmetric_diagonal":
        print("neutrino normalization class drifted from symmetric_diagonal", file=sys.stderr)
        return 1

    certificate = dict(nu.get("residual_factorization_certificate", {}))
    if certificate.get("entrywise_amplitude_free", True):
        print("neutrino residual factorization allows a free entrywise amplitude", file=sys.stderr)
        return 1

    if "K_core_majorana_sym" not in nu:
        print("missing explicit majorana symmetric kernel", file=sys.stderr)
        return 1

    print("neutrino residual factorization is explicit and bounded")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
