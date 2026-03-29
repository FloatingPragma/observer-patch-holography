#!/usr/bin/env python3
"""Keep the quark zero-odd-scalar corollary gated by the shared charged law."""

from __future__ import annotations

import argparse
import json
import pathlib
import sys


ROOT = pathlib.Path(__file__).resolve().parents[2]
DEFAULT_ODD_FORM = ROOT / "particles" / "runs" / "flavor" / "charged_dirac_odd_deformation_form.json"
DEFAULT_RESPONSE = ROOT / "particles" / "runs" / "flavor" / "quark_odd_response_law.json"


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate the quark zero-odd-scalar corollary gate.")
    parser.add_argument("--odd-form", default=str(DEFAULT_ODD_FORM))
    parser.add_argument("--response", default=str(DEFAULT_RESPONSE))
    args = parser.parse_args()

    odd_form = json.loads(pathlib.Path(args.odd_form).read_text(encoding="utf-8"))
    response = json.loads(pathlib.Path(args.response).read_text(encoding="utf-8"))
    corollary = str(odd_form.get("quark_zero_odd_scalar_corollary", "open"))
    status = str(response.get("delta_logg_q_status", ""))
    if corollary != "closed" and status == "closed_zero_corollary":
        print("quark odd response closed the zero-odd-scalar corollary before the charged law closed", file=sys.stderr)
        return 1
    if odd_form.get("odd_scalar_slot_present") is False and response.get("delta_logg_q") != 0.0:
        print("quark odd response reintroduced an odd scalar slot even though the odd codomain excludes it", file=sys.stderr)
        return 1
    print("quark zero-odd-scalar corollary guard passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
