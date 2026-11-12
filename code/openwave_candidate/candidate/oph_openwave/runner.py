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
OPENWAVE_COMMIT = "ef55175d04737e73c92f33f9f9cf3d17c593da54"
MODULES = ("charge", "spin", "leptons", "conjugation", "dispersion")


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
                     "command": f"python -m oph_openwave.{name}", "exit_code": 0,
                     "verdict": "PASS_SCOPED_COMPUTATION", "result": result})
    sources = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
               for p in sorted(Path(__file__).parent.glob("*.py"))}
    return {"schema": "oph.openwave_candidate.run.v1",
            "written_utc": datetime.now(timezone.utc).isoformat(),
            "machine": "workstation", "source_basis_commit": RER_COMMIT,
            "openwave_criteria_commit": OPENWAVE_COMMIT,
            "candidate_python_sha256": sources,
            "versions": {"python": platform.python_version(),
                         **{n: importlib.metadata.version(n) for n in ("numpy", "sympy", "mpmath")}},
            "jobs": rows, "official_openwave_admission": False,
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
