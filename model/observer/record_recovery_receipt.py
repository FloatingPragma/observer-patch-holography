"""Reproduce complete-record recovery, encoded memory and identifying schedules."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from . import record_recovery as recovery
from . import record_recovery_schedule as schedule
from . import record_recovery_source as source


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = (
    "model/observer/native.py",
    "model/observer/record_recovery.py",
    "model/observer/record_recovery_schedule.py",
    "model/observer/record_recovery_source.py",
    "model/observer/source_discrimination/circuit.py",
    "model/observer/source_scalar/compiler.py",
    "code/source_selection_model/response.json",
    "code/m1_source_realization/receipt.json",
    "model/requirements.txt",
    "model/common_source/slow_geometry.py",
    "evidence/native_geometric_source_20260925/spec.json",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py",
)


def canonical(value):
    return json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"


def read_receipt(path):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate receipt key: {key}")
            result[key] = value
        return result

    def invalid(value):
        raise ValueError(f"nonfinite receipt value: {value}")

    return json.loads(Path(path).read_text(encoding="utf-8"),
                      object_pairs_hook=unique, parse_constant=invalid)


def compare_receipt(recorded, actual):
    # Canonical bytes bind every field and distinguish bool/int/float.
    if canonical(recorded) != canonical(actual):
        raise ValueError("exact record recovery receipt mismatch")


def run():
    return {
        "schema": "oph.observer.record_recovery.v1",
        "source_sha256": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                          for name in SOURCE_PATHS},
        "existing_observer": {
            domain: {law: recovery.observer_report(law, binary=binary) for law in ("C", "E")}
            for domain, binary in (("full", False), ("binary", True))},
        "uniform_independent_schedule": {
            domain: {law: recovery.iid_report(law, binary=binary) for law in ("C", "E")}
            for domain, binary in (("full", False), ("binary", True))},
        "one_event_encoded_qubit": {law: recovery.one_event_code(law) for law in ("C", "E")},
        "identifying_schedules": schedule.audit(),
        "source_recovery": source.audit(),
        "contract": {
            "candidate_revision": "f0494b54b8e2f695e46065bcd83802a8bc8f6940",
            "carrier": "native total-two twelve-port source: 78 full or 66 binary states",
            "record": "complete classical seam, fair coin and native current/contrast",
            "between_events": "native repair and classical record-only seam selection",
            "recovery": "arbitrary final CPTP decoder of final source and complete classical tape",
            "supplied_services": "blank scratch, classical archive, coherent code controls and readout",
            "time": "repair attempts; no physical clock calibration",
        },
        "scope": {
            "classical_initial_loads_exactly_recoverable_with_final_source": True,
            "unknown_quantum_input_generally_recoverable": False,
            "one_event_encoded_qubit_recoverable": True,
            "all_quantum_codes_or_active_error_correction_excluded": False,
            "coherent_uncommitted_tape_excluded": False,
            "new_physical_law_or_observations": False,
            "programme_wide_no_go": False,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    files = parser.add_mutually_exclusive_group()
    files.add_argument("--write", type=Path)
    files.add_argument("--check", type=Path)
    args = parser.parse_args()
    result = run()
    encoded = canonical(result)
    if args.write:
        with args.write.open("x", encoding="utf-8", newline="\n") as stream:
            stream.write(encoded)
    if args.check:
        compare_receipt(read_receipt(args.check), result)
    print(json.dumps({"exact_replay": True,
                      "result_sha256": hashlib.sha256(encoded.encode()).hexdigest()}, sort_keys=True))


if __name__ == "__main__":
    main()
