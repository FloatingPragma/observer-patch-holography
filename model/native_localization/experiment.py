"""Reproduce and source-bind the native localization successor packet."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from tools import strict_json
from . import localization
from ..observer.scalar_response.experiment import compare_receipt

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE/"results.json"


def bindings():
    result = localization.clock.source_sha256()
    paths = ["model/native_localization/__init__.py",
             "model/native_localization/localization.py",
             "model/native_localization/experiment.py",
             "model/native_localization/test_localization.py",
             "model/native_localization/README.md",
             "model/common_source/record_code.py",
             "model/observer/adaptive_record_recovery.py",
             "model/federated_quantum_sector/construction.py",
             "model/federated_quantum_sector/README.md",
             "model/observer/QUANTUM_MIXING.md"]
    result.update({name: hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in paths})
    return dict(sorted(result.items()))


def report():
    return {**localization.report(), "source_sha256": bindings()}


def verify(retained, current=None):
    """Compare every retained field; observable floats have fixed 2e-10 tolerance.

    No source budget, record probability, mean, spectrum, or entropy field is
    dropped. Structure, types, exact fractions, counts and source hashes are
    exact. The shared comparator rejects every nonfinite floating value.
    """
    current = report() if current is None else current
    compare_receipt(current, retained)
    for value in (current, retained):
        diagnostics = value["source_execution"]
        if (diagnostics["floating_diagnostic_gate"] != 1e-8
                or not 0 <= diagnostics["full_history_error_sum_frobenius"] < 1e-8
                or abs(diagnostics["complete_history_mass"]-1) > 1e-8
                or diagnostics["complete_raw_history_count"] != 128
                or value["canonical_physical_premises_discharged"] != 0
                or value["CP03_accepted"] is not False):
            raise ValueError("localization diagnostic or physical boundary changed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", nargs="?", const=OUTPUT, type=Path)
    parser.add_argument("--check", nargs="?", const=OUTPUT, type=Path)
    args = parser.parse_args()
    if (args.write is None) == (args.check is None):
        parser.error("choose exactly one of --write or --check")
    current = report()
    if args.write is not None:
        args.write.write_text(json.dumps(current, indent=2, sort_keys=True, allow_nan=False)+"\n")
        print(f"Wrote {args.write}")
    else:
        retained = strict_json.load(args.check)
        verify(retained, current)
        print("Native localization exact certificate, source replay, and bindings verified")


if __name__ == "__main__":
    main()
