#!/usr/bin/env python3
"""Guard the universal charged scalarization law against label-dependent closure."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_INPUT = ROOT / "particles" / "runs" / "flavor" / "charged_budget_transport.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate charged scalarization label blindness.")
    parser.add_argument("--input", default=str(DEFAULT_INPUT), help="Input charged-budget artifact.")
    args = parser.parse_args()

    payload = json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    certificate = dict(payload.get("charged_dirac_scalarization_certificate", {}))
    law_scope = str(certificate.get("law_scope", ""))
    label_blindness_status = str(certificate.get("label_blindness_status", ""))
    label_blindness_candidate = bool(certificate.get("label_blindness_candidate", False))

    if law_scope != "direct_sum_u_plus_d_plus_e_pre_normal_form":
        print("charged scalarization law scope is not the universal u+d+e direct sum", file=sys.stderr)
        return 1
    if label_blindness_status == "closed" and not label_blindness_candidate:
        print("closed label blindness claimed without a label-blind candidate certificate", file=sys.stderr)
        return 1

    print("scalarization label-blindness guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
