#!/usr/bin/env python3
"""Validate the OPH icosahedral screen-sieve certificate."""

from __future__ import annotations

import json
import pathlib
import sys
from typing import Any


def validate(certificate: dict[str, Any]) -> dict[str, bool]:
    poly = {item["name"]: item for item in certificate.get("polyhedral_comparison", [])}
    checks = dict(certificate.get("checks", {}))
    checks.update(
        {
            "status_is_closed": certificate.get("status")
            == "closed_on_declared_triangulated_screen_branch",
            "all_examples_have_total_charge_12": all(
                item.get("total_charge") == 12 for item in poly.values()
            ),
            "costs_are_ordered": (
                poly.get("icosahedral", {}).get("defect_cost_sum_q2") == 12
                and poly.get("octahedral", {}).get("defect_cost_sum_q2") == 24
                and poly.get("tetrahedral", {}).get("defect_cost_sum_q2") == 36
            ),
            "orbit_size_is_12": certificate.get("orbit_stabilizer", {}).get("orbit_size") == 12,
            "projection_is_p_over_12": certificate.get(
                "capacity_electroweak_projection", {}
            ).get("gamma_EW")
            == "(P/12)*log(N/pi)",
        }
    )
    return checks


def main(path: str) -> int:
    cert_path = pathlib.Path(path)
    certificate = json.loads(cert_path.read_text(encoding="utf-8"))
    checks = validate(certificate)
    payload = {"checks": checks, "pass": all(checks.values())}
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0 if payload["pass"] else 1


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: validate_screen_sieve_certificate.py CERTIFICATE", file=sys.stderr)
        raise SystemExit(2)
    raise SystemExit(main(sys.argv[1]))
