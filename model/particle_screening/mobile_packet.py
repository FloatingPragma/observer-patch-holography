"""Replay the complete moving-memory contract and its source implementation."""
from __future__ import annotations

import argparse
import ast
from fractions import Fraction as F
import hashlib
import json
import math
from pathlib import Path

from tools import strict_json
from ..screened_recovery import experiment as predecessor
from . import mobile, mobile_service


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "mobile_results.json"
CONTRACT_SHA256 = "a1f15e50b8af96928d327ae16e0159d2242b9052fc62b185035205f77378052e"
CAMPAIGN = "mobile-screened-memory-20261029-v1"
OWN_SOURCES = (
    "model/__init__.py",
    "model/observer/scalar_response/__init__.py",
    "model/particle_screening/mobile.py",
    "model/particle_screening/mobile_service.py",
    "model/particle_screening/mobile_packet.py",
    "model/particle_screening/mobile_contract.json",
    "model/particle_screening/test_mobile.py",
    "model/particle_screening/test_mobile_service.py",
    "model/particle_screening/test_mobile_packet.py",
    "evidence/native_geometric_source_20260925/spec.json",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py",
)


def contract():
    path = HERE / "mobile_contract.json"
    if hashlib.sha256(path.read_bytes()).hexdigest() != CONTRACT_SHA256:
        raise ValueError("mobile memory contract changed")
    return strict_json.load(path)


def source_bindings():
    # The geometry loader imports an archived package outside the predecessor's
    # model/tools closure. Include that package's eager initializers as well as
    # recursive imports, independently of the current process's module caches.
    # No predecessor receipt is accepted as a replacement for this execution.
    roots = ((mobile.slow_geometry.ARCHIVE, {"oph_exact", "oph_fpe"}),
             (ROOT, {"model", "tools"}))
    pending = [ROOT / name for name in (*predecessor.source_bindings(), *OWN_SOURCES)]
    seen = set()
    while pending:
        path = pending.pop()
        if path in seen:
            continue
        seen.add(path)
        if path.suffix != ".py":
            continue
        owner = next(root for root, _ in roots if path.is_relative_to(root))
        parts = path.relative_to(owner).with_suffix("").parts
        package = parts[:-1]
        # Importing a.b.c also executes a/__init__.py and a/b/__init__.py.
        pending.extend(init for i in range(1, len(package)+1)
                       if (init := owner.joinpath(*package[:i], "__init__.py")).is_file())
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            candidates = []
            if isinstance(node, ast.Import):
                candidates.extend(tuple(alias.name.split(".")) for alias in node.names)
            elif isinstance(node, ast.ImportFrom):
                base = package[:len(package)-node.level+1] if node.level else ()
                base += tuple(node.module.split(".")) if node.module else ()
                candidates.append(base)
                candidates.extend(base+tuple(alias.name.split("."))
                                  for alias in node.names if alias.name != "*")
            for candidate in candidates:
                for root, prefixes in roots:
                    if candidate and candidate[0] in prefixes:
                        module = root.joinpath(*candidate)
                        pending.extend(found for found in (module.with_suffix(".py"), module / "__init__.py")
                                       if found.is_file())
    return {path.relative_to(ROOT).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(seen)}


def admission_intervals():
    """Fixed gates, not fitted residuals or rigorous floating-error enclosures."""
    bounds = {}
    for key in ("maximum_encoder_operator_error", "maximum_decoder_operator_error"):
        bounds[("source_composition", key)] = F(1, 10**8)
    for law in ("C", "E"):
        prefix = ("source_composition", "laws", law)
        bounds[(*prefix, "maximum_attempt_log_error")] = F(1, 10**7)
        for n in (1, 2, 10, 100):
            bounds[(*prefix, "encoded_t_attempt_operator_error_upper", str(n))] = F(1, 10**5)
    for i in range(2):
        bounds[("source_replay", "source_executed_prototypes", i, "source_encoder_error")] = F(1, 10**8)
    for name in ("motion_then_forbidden_reset", "closed_core_exchange"):
        for key in ("maximum_terminal_source_operator_error", "selected_terminal_decoder_error"):
            bounds[("source_replay", "complete_stopped_histories", name, key)] = F(1, 10**8)
    return bounds


def admit(diagnostics):
    bounds, seen = admission_intervals(), set()

    def visit(value, path=()):
        if path in bounds:
            upper = bounds[path]
            if type(value) is not float or not math.isfinite(value) or not 0 <= F(value) <= upper:
                raise ArithmeticError(f"fresh mobile source diagnostic outside fixed gate at {path}")
            seen.add(path)
            return {"kind": "fresh_numerical_admission", "lower": "0", "upper": str(upper)}
        if type(value) is dict:
            return {key: visit(item, (*path, key)) for key, item in value.items()}
        if type(value) in (list, tuple):
            return [visit(item, (*path, i)) for i, item in enumerate(value)]
        if type(value) is float:
            raise ArithmeticError(f"unregistered mobile source diagnostic at {path}")
        return value

    result = visit(diagnostics)
    if seen != set(bounds):
        raise ArithmeticError("missing fresh mobile source diagnostic")
    return result


def report():
    frozen = contract()
    pins = source_bindings()
    exact = mobile.exact_certificate()
    executed = admit(mobile_service.report())
    if source_bindings() != pins:
        raise ArithmeticError("mobile source files changed during replay")
    return {
        "schema": "oph.mobile_screened_memory.receipt.v1",
        "campaign_id": CAMPAIGN,
        "contract_sha256": CONTRACT_SHA256,
        "scope": frozen["scope"],
        "exact": exact,
        "execution": executed,
        "numerical_admission_semantics": (
            "Fresh source execution and counted operator-composition diagnostics under fixed gates; "
            "not outward interval proofs, physical hardware/reset error bounds, or infinite-time accuracy."
        ),
        "integrated_gate_accepted": False,
        "physical_comparison_count_added": 0,
        "physical_particle_identification": False,
        "source_sha256": pins,
    }


def compare(actual, expected):
    # Canonical JSON distinguishes false from 0, rejects NaN, and preserves
    # the entire record/cost/scope tree rather than matching headline scalars.
    encode = lambda value: json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    if encode(actual) != encode(expected):
        raise ValueError("mobile memory receipt differs from complete source replay")


def verify(path=OUTPUT):
    compare(strict_json.load(Path(path)), report())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument("--write", nargs="?", const=OUTPUT, type=Path)
    modes.add_argument("--check", nargs="?", const=OUTPUT, type=Path)
    args = parser.parse_args()
    value = report()
    if args.write:
        args.write.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n",
                              encoding="utf-8", newline="\n")
    elif args.check:
        compare(strict_json.load(args.check), value)
    print("Verified moving memory: complete native flags, optimal core resets, and compiled recovery.")


if __name__ == "__main__":
    main()
