"""Lift the exact source/bit reference to the actual full-record observer.

This constructs a stationary finite source/FIFO marginal. The exported tape
keeps growing and is not stationary. No modular implementation is supplied.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

from .family import FamilyParameters, binding
from . import observer, stationary_reference
try:
    from tools.pin_scope import bound_pins
except ModuleNotFoundError:  # run as a script: put the repository root on sys.path
    import sys as _sys
    _sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parents[2]))
    from tools.pin_scope import bound_pins


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def source_bit_projection(blocks):
    result = defaultdict(F)
    for memory, density in blocks.items():
        bit = int(memory[-1].current > 0)
        for (i, a, j, b), value in density.items():
            if (i, a) != (j, b) or a != 0:
                raise ValueError("diagonal source, no reference factor required for this projection")
            result[observer.STATES[i], bit] += value
    return {key: value for key, value in result.items() if value}


def seed_blocks(stationary_joint):
    """Only the last-current bit affects the next instrument. These two
    representative records embed that quotient for one step; they are NOT
    asserted to be the stationary full-record distribution.
    """
    index = {state: i for i, state in enumerate(observer.STATES)}
    blocks = {}
    for row in stationary_joint:
        probability = F(row["probability"])
        if not probability:
            continue
        bit = row["memory"]
        record = observer.Record("repair", "feedback", observer.native.E, "C", 0, bit, bit)
        i = index[tuple(row["state"])]
        blocks.setdefault((record,), {})[i, 0, i, 0] = probability
    return blocks


def lift(stationary_joint, c_weight):
    params = FamilyParameters(c_weight=c_weight, repair_probability=F(1),
                              exploration_probability=F(1, 2), memory_capacity=1,
                              overflow="export_oldest")
    seed = seed_blocks(stationary_joint)
    result = observer.advance_marginal(seed, params)
    expected = {(tuple(row["state"]), row["memory"]): F(row["probability"])
                for row in stationary_joint if F(row["probability"])}
    if source_bit_projection(result) != expected:
        raise ArithmeticError("full-record reference does not project to the solved bit reference")
    if observer.advance_marginal(result, params) != result:
        raise ArithmeticError("actual full-record observer reference is not stationary")
    if observer.trace(observer.summed_density(result)) != 1:
        raise ArithmeticError("full-record reference is not normalized")
    return result


def serialize_record(record):
    return [record.kind, record.selector, list(record.seam) if record.seam else None,
            record.law, record.coin, record.label, record.current]


def certify(blocks, c_weight):
    params = FamilyParameters(c_weight=c_weight, repair_probability=F(1),
                              exploration_probability=F(1, 2), memory_capacity=1,
                              overflow="export_oldest")
    if observer.trace(observer.summed_density(blocks)) != 1:
        raise ArithmeticError("reference is not normalized")
    if observer.advance_marginal(blocks, params) != blocks:
        raise ArithmeticError("reference is not stationary under the declared full-record observer")
    rows, supports = [], []
    for memory, density in sorted(blocks.items(), key=lambda item: repr(item[0])):
        record = serialize_record(memory[0])
        support = 0
        for (i, a, j, b), probability in sorted(density.items()):
            if (i, a) != (j, b) or a != 0 or probability <= 0:
                raise ArithmeticError("lift ceased to be a positive diagonal reference")
            rows.append([record, i, str(probability)])
            support += 1
        supports.append({"record": record, "positive_source_dimensions": support,
                         "record_probability": str(observer.trace(density))})
    payload = json.dumps(rows, separators=(",", ":"), sort_keys=True).encode()
    return {"normalization": "1", "stationarity_residual": "0",
            "positive_full_record_blocks": len(blocks),
            "positive_joint_diagonal_dimensions": len(rows),
            "canonical_exact_reference_sha256": hashlib.sha256(payload).hexdigest(),
            "reference_encoding": "sorted [serialized_record,source_index,rational_probability] rows; reconstructed exactly by lift()",
            "record_block_supports": supports,
            "faithful_on_full_M78_for_each_record": False,
            "faithful_on_declared_direct_sum_of_positive_supports": True}


def report():
    quotient = stationary_reference.run()
    endpoint_certificates = {}
    for law, c in (("C", F(1)), ("E", F(0))):
        blocks = lift(quotient["stationary_joint"], c)
        endpoint_certificates[law] = certify(blocks, c)
        # lift explicitly verified this equality against the solved quotient.
        endpoint_certificates[law]["source_bit_projection_equals_exact_stationary_law"] = True
    files = [Path(__file__), HERE/"observer.py", HERE/"stationary_reference.py",
             HERE/"adaptive_reference.py", HERE/"family.py"]
    return {"schema": "oph.completion.observer_reference.v1", "gate": "CP-02",
            "status": "EXACT_FINITE_REFERENCE_JOIN_GATE_OPEN", "campaign": binding(),
            "source_sha256": bound_pins({path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
                              for path in files}),
            "quotient_reference": "stationary_reference.run(): exact129-support source/bit stationary law",
            "observer_parameters": {"repair_probability": "1", "exploration_probability": "1/2",
                                    "memory_capacity": 1, "overflow": "export_oldest"},
            "endpoint_certificates": endpoint_certificates,
            "boundaries": [
                "Stationarity belongs to the finite source and last complete record marginal, with all previous records exported and traced out in this marginal.",
                "The complete history has growing record storage and is not asserted to be stationary.",
                "The old bit law is lifted through the actual C/E observer instrument; it is not relabelled as a full-record law.",
                "Each positive record block is faithful only after explicitly restricting the source algebra to its positive diagonal support.",
                "No physically available full support algebra, preparation accuracy, modular action, clock comparison or energy normalization is constructed.",
                "The lift covers the r=1 C and E endpoints; a hold writes current0 in this observer and needs its own joint stationary solve.",
                "No integrated CP-01 or CP-02 gate is accepted."
            ]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    value = report()
    path = HERE/"observer_reference_results.json"
    if args.check:
        from tools import strict_json
        if strict_json.load(path) != value:
            raise SystemExit("observer reference receipt differs")
    if args.write:
        path.write_text(json.dumps(value, indent=2, sort_keys=True)+"\n")
    print(json.dumps({law: {key: row[key] for key in ("positive_full_record_blocks",
        "positive_joint_diagonal_dimensions", "stationarity_residual")}
        for law, row in value["endpoint_certificates"].items()}, indent=2))


if __name__ == "__main__":
    main()
