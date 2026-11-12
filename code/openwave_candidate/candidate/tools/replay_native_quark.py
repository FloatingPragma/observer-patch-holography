"""Recompute the actual pinned quark measures, preserving their negative result.

This optional source replay imports the public RER source at the explicit
--rer path. It is separate from the self-contained candidate installation.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys
import time

PIN = "3bcacf1276ddb3ead35be2b9adcd63b6643e2d35"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rer", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.rer.resolve()
    head = subprocess.check_output(["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip()
    if head != PIN:
        raise ValueError("pinned public source commit required")
    subprocess.run(["git", "-C", str(root), "diff", "--exit-code", PIN, "--", "model", "code"], check=True)
    sys.path.insert(0, str(root))
    from model.native_quark_response.response import positive_response
    from model.coherent_defect_wave.chiral_joint.charged_ground_exchange.resolvent import field_data
    import sympy as s
    start = time.monotonic()
    up = positive_response("uc")
    down = positive_response("dc")
    checks = {key: up[key] == down[key] for key in ("moments", "numerator", "denominator", "visible_support")}
    checks["thirteen_positive_records_per_channel"] = all(
        len(data["records"]) == 13 and all(row["weight_interval"][0] > 0 for row in data["records"])
        for data in (up, down))
    field = field_data()[2]
    checks["complete_success_weight_one_half"] = all(
        data['moments'][0] == field.convert(s.Rational(1, 2)) for data in (up, down))
    shift = s.Rational(1, 100)
    broken_numerator = down["numerator"].shift(shift)
    broken_denominator = down["denominator"].shift(shift)
    mutations = {"shifted_down_energy_breaks_measure_equality":
                 up["numerator"] * broken_denominator != broken_numerator * up["denominator"]}
    if not all(value is True for value in (*checks.values(), *mutations.values())):
        raise ArithmeticError("native quark degeneracy replay or mutation control failed")
    roots = [{"excitation_interval": list(map(str, row["difference_interval"])),
              "weight_interval": list(map(str, row["weight_interval"]))} for row in up["records"]]
    result = {"schema": "oph.openwave.native_quark_replay.v1", "source_commit": head,
              "written_utc": datetime.now(timezone.utc).isoformat(), "machine": "workstation",
              "versions": {"python": platform.python_version(), "sympy": importlib.metadata.version("sympy")},
              "command": "python candidate/tools/replay_native_quark.py --rer <pinned-rer-checkout> --output evidence/native_quark.json",
              "wall_seconds": time.monotonic() - start, "exit_code": 0,
              "checks": checks, "mutations": mutations, "complete_positive_records": roots,
              "success_total_weight_exact": "1/2", "pump_failure_weight_exact": "1/2",
              "outside_support_weight_exact": "0",
              "method_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "verdict": "CONFIRMED_EQUAL_SOURCE_MEASURES_PHYSICAL_SPLITTING_NOT_DERIVED",
              "scope": "Exact uc/dc ground-selected local energy measures of the pinned equal-coupling native source; no physical quark masses, propagation or confinement are inferred.",
              "full_criterion_met": False}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(result["verdict"])


if __name__ == "__main__":
    main()
