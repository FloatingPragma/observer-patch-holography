"""Replay the native exchange classification and its source-operated readouts."""
from __future__ import annotations

import argparse
import ast
import hashlib
import json
from pathlib import Path

from tools import strict_json
from ..common_source import slow_geometry
from . import statistics, statistics_service, statistics_topology


ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "statistics_results.json"
CAMPAIGN = "native-exchange-statistics-20261029-v1"
CONTRACT_SHA256 = 'fe979eb6e059f059d3196b785bdf578e71371bcc60fc6e3af68e6862aea1349b'
OWN_SOURCES = (
    "model/exchange_sector/statistics.py",
    "model/exchange_sector/statistics_service.py",
    "model/exchange_sector/statistics_topology.py",
    "model/exchange_sector/statistics_packet.py",
    "model/exchange_sector/statistics_contract.json",
    "model/exchange_sector/test_statistics.py",
    "model/exchange_sector/test_statistics_service.py",
    "model/exchange_sector/test_statistics_topology.py",
    "model/exchange_sector/test_statistics_packet.py",
    "model/exchange_sector/STATISTICS.md",
    "evidence/native_geometric_source_20260925/spec.json",
    "evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py",
    "code/source_selection_model/response.json",
    "code/m1_source_realization/receipt.json",
    "Lean/Screen/HolonomyInterference.lean",
)


def source_bindings():
    """Static import closure, including eager package initializers and geometry.

    Computing this from files makes custody independent of Python's module
    cache. The geometry's archived import root is resolved before the outer
    repository; its dynamically imported carrier is an explicit seed.
    """
    roots = ((slow_geometry.ARCHIVE, {"oph_exact", "oph_fpe"}),
             (ROOT, {"model", "tools"}))
    pending = [ROOT / name for name in OWN_SOURCES]
    seen = set()
    while pending:
        path = pending.pop()
        if path in seen:
            continue
        if not path.is_file() or not path.is_relative_to(ROOT):
            raise ValueError(f"existing repository source required: {path}")
        seen.add(path)
        if path.suffix != ".py":
            continue
        owner = next(root for root, _ in roots if path.is_relative_to(root))
        package = path.relative_to(owner).with_suffix("").parts[:-1]
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


def compare(actual, expected):
    """Compare every output, cost, scope and source leaf, with exact JSON types."""
    encode = lambda value: json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    if encode(actual) != encode(expected):
        raise ValueError("exchange statistics receipt differs from complete replay")


# A valid new manifest does not reload modules already held by Python.
# Canonical replay uses a fresh interpreter; a repeated API call must fail
# closed if its source tree has changed since this packet was imported.
_IMPORTED_SOURCE_BINDINGS = source_bindings()


def contract():
    path = HERE / "statistics_contract.json"
    if hashlib.sha256(path.read_bytes()).hexdigest() != CONTRACT_SHA256:
        raise ValueError("exchange statistics contract changed")
    return strict_json.load(path)


def report():
    frozen = contract()
    pins = source_bindings()
    # Importing the instrument already warms this geometry cache.  Recheck
    # the retained source/spec binding even when report() is called again
    # in that interpreter; hashing changed bytes is not their validation.
    slow_geometry.source.cache_clear()
    slow_geometry.source()
    if pins != _IMPORTED_SOURCE_BINDINGS:
        raise ValueError('exchange sources changed since import; use a fresh interpreter')
    topology = statistics_topology.certificate()
    instruments = statistics.certificate()
    execution = statistics_service.certificate()
    if source_bindings() != pins:
        raise ArithmeticError("exchange statistics source files changed during replay")
    return {
        "schema": "oph.native_exchange_statistics.receipt.v1",
        "campaign_id": CAMPAIGN,
        "contract_sha256": pins["model/exchange_sector/statistics_contract.json"],
        "scope": frozen["scope"],
        "topology": topology,
        "instruments": instruments,
        "execution": execution,
        "integrated_gate_accepted": False,
        "physical_comparison_count_added": 0,
        "physical_particle_identification": False,
        "source_sha256": pins,
    }


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
    print("Verified exchange topology, complete signed instruments and source readouts.")


if __name__ == "__main__":
    main()
