"""Run small scientific modules and record their complete local verdicts."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import time

RER_COMMIT = "3bcacf1276ddb3ead35be2b9adcd63b6643e2d35"
OPENWAVE_COMMIT = "5d686cb049b82cc55f8822210a254d5833be23a7"
MODULES = ("charge", "spin", "leptons", "conjugation", "dispersion", "construction", "quantum_records", "observer_protocol", "record_preparation", "record_energy", "interacting_records", "flux_controls", "regional_charge", "general_cut", "modular_native", "sparse_native", "uniform_native", "compact_native", "pinned_charge", "native_currents", "many_load_native", "fiber_preparation")


def validate_result(name: str, result: dict) -> None:
    """A scope declaration and non-vacuous, literal Boolean checks are required."""
    if not isinstance(result, dict) or not result.get("scope"):
        raise ValueError(f"{name}: scientific scope missing")
    checks = result.get("checks")
    if not isinstance(checks, dict) or not checks:
        raise ValueError(f"{name}: nonempty check map required")
    if any(value is not True for value in checks.values()):
        raise ValueError(f"{name}: scientific check failed")
    mutations = result.get("mutations")
    if not isinstance(mutations, dict) or not mutations:
        raise ValueError(f"{name}: broken-control checks missing")
    if any(value is not True for value in mutations.values()):
        raise ValueError(f"{name}: a broken control was not detected")
    if result.get("full_criterion_met") is not False:
        raise ValueError(f"{name}: unsupported full-criterion promotion")


def execute(names: list[str]) -> dict:
    rows = []
    for name in names:
        if name not in MODULES:
            raise ValueError(f"unknown computation: {name}")
        module = importlib.import_module(f"oph_openwave.{name}")
        started = time.monotonic()
        result = module.run()
        validate_result(name, result)
        rows.append({"task": name, "wall_seconds": time.monotonic() - started,
                     "command": f"python -m oph_openwave.runner {name}", "exit_code": 0,
                     "verdict": "PASS_SCOPED_COMPUTATION", "result": result})
    sources = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
               for p in sorted(Path(__file__).parent.glob("*.py"))}
    return {"schema": "oph.openwave_candidate.run.v1",
            "written_utc": datetime.now(timezone.utc).isoformat(),
            "machine": "workstation", "legacy_source_basis_commit": RER_COMMIT,
            "construction_source_identity": "candidate_python_sha256",
            "openwave_criteria_commit": OPENWAVE_COMMIT,
            "candidate_python_sha256": sources,
            "versions": {"python": platform.python_version(),
                         **{n: importlib.metadata.version(n) for n in ("numpy", "sympy", "mpmath")}},
            "jobs": rows, "official_model_id": "M10", "official_openwave_admission": True,
            "official_openwave_cells_awarded": 0,
            "verdict": "PASS_SCOPED_COMPUTATIONS_NO_PHYSICAL_PROMOTION"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("task", choices=(*MODULES, "all"))
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = execute(list(MODULES) if args.task == "all" else [args.task])
    payload = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload)
    else:
        sys.stdout.write(payload)


if __name__ == "__main__":
    main()
