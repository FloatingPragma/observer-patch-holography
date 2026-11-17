"""Reproduce the exact native quantum mixing and memory certificates.

No numerical eigenvalues, sampled counts or fitted convergence thresholds enter
this receipt. Its scope is the independent uniform-seam forgotten-record law.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from . import quantum_mixing, quantum_memory
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins


ROOT = Path(__file__).resolve().parents[2]
SOURCE_PATHS = (
    "model/observer/native.py",
    "model/observer/quantum_mixing.py",
    "model/observer/quantum_memory.py",
    "model/observer/source_discrimination/circuit.py",
    "model/observer/source_discrimination/certificate.py",
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

    def invalid_constant(value):
        raise ValueError(f"nonfinite receipt constant: {value}")

    return json.loads(Path(path).read_text(encoding="utf-8"),
                      object_pairs_hook=unique, parse_constant=invalid_constant)


def compare_receipt(recorded, actual):
    # Canonical JSON comparison distinguishes bool/int and int/float, and checks
    # every key, list entry, rational, source pin and scope field exactly.
    if canonical(recorded) != canonical(actual):
        raise ValueError("exact quantum mixing receipt mismatch")


def carrier_spectrum():
    """Check the actual incidence spectrum with an exact integer polynomial."""
    import sympy as sp

    laplacian = sp.zeros(12)
    for a, b in quantum_mixing.seams():
        laplacian[a, a] += 1
        laplacian[b, b] += 1
        laplacian[a, b] -= 1
        laplacian[b, a] -= 1
    variable = sp.Symbol("x")
    polynomial = laplacian.charpoly(variable).as_expr()
    target = variable * (variable - 6)**5 * (variable**2 - 10*variable + 20)**3
    if sp.expand(polynomial - target) != 0:
        raise ArithmeticError("the bound source incidence has the wrong exact spectrum")
    return {"characteristic_polynomial": "x*(x-6)^5*(x^2-10*x+20)^3",
            "integer_coefficients_descending": [int(x) for x in sp.Poly(polynomial, variable).all_coeffs()],
            "positive_gap": "5-sqrt(5)",
            "classical_exclusion_rate_factor": "1/60",
            "single_event_C_E_quantum_gap": "(5-sqrt(5))/60"}


def run():
    return {
        "schema": "oph.observer.quantum_mixing.v1",
        "source_sha256": bound_pins({name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                          for name in SOURCE_PATHS}),
        "sectors": quantum_mixing.audit(),
        "carrier_spectrum": carrier_spectrum(),
        "memory": quantum_memory.audit(),
        "contract": {
            "candidates": ["C", "E"],
            "carrier": "fixed twelve-port thirty-seam source, total load two",
            "schedule": "independent uniform seam choices and fair placement coins",
            "channel": "trace out each committed repair record; no feedback or postselection",
            "time": "number of repair attempts, with no physical clock identification",
            "classical_gap_import": "Caputo, Liggett, Richthammer (2010), section 4.1.1",
            "separability_import": "Gurvits, Barnum (2002), Corollary 3",
        },
        "scope": {
            "binary_and_full_total_two_quantum_gap": "(5-sqrt(5))/60",
            "binary_fixed_observables": "scalar multiples of the identity",
            "full_stationary_state": "uniform binary-shell density, rank 66 on dimension 78",
            "binary_eventual_entanglement_breaking": True,
            "full_78_finite_time_entanglement_breaking": False,
            "retained_record_instrument_replaced_by_forgotten_channel": False,
            "adaptive_controller_or_protected_response_registers_covered": False,
            "source_controls_select_C_over_E": False,
            "physical_data_or_programme_wide_exclusion": False,
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
                      "quantum_gap": result["scope"]["binary_and_full_total_two_quantum_gap"],
                      "result_sha256": hashlib.sha256(encoded.encode()).hexdigest()}, sort_keys=True))


if __name__ == "__main__":
    main()
