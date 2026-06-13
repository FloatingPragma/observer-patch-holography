#!/usr/bin/env python3
"""Validate the OPH issue #332 RG/Higgs naturality certificate."""

from __future__ import annotations

import json
import pathlib
import sys


def main(path: str = "issue_332_rg_naturality_certificate.json") -> int:
    cert_path = pathlib.Path(path)
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    optional = cert.get("optional_upstream_resonance_check", {})
    checks = {
        "issue_is_332": cert.get("issue") == 332,
        "accepted": cert.get("accepted") is True,
        "epsilon_H_zero": cert.get("epsilon_H") == "0",
        "epsilon_interval_zero": cert.get("epsilon_H_interval") == ["0", "0"],
        "weak_scale_forbidden": any("measured weak scale" in item for item in cert.get("forbidden_calibrations", [])),
        "higgs_mass_forbidden": any("Higgs" in item for item in cert.get("forbidden_calibrations", [])),
        "diagnostic_N_not_source": optional.get("n_crc_source") != "provided",
        "strict_resonance_not_required": optional.get("strict_resonance") is False,
    }
    payload = {"checks": checks, "pass": all(checks.values())}
    print(json.dumps(payload, indent=2))
    return 0 if payload["pass"] else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "issue_332_rg_naturality_certificate.json"))
