#!/usr/bin/env python3
"""Validate the issue #320 collar Poisson witness certificate.

The validator recomputes every declared family from the certificate inputs
with the same exact arithmetic as the generator and checks that the stored
total-variation values, bounds, and pass flags match the recomputation.
"""

from __future__ import annotations

import json
import pathlib
import sys
from decimal import Decimal
from fractions import Fraction

SCRIPTS = pathlib.Path(__file__).resolve().parent
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))

import collar_poisson_witness as witness  # noqa: E402

RECOMPARE_TOLERANCE = Decimal("1e-35")


def parse_frac(text: str) -> Fraction:
    num, den = text.split("/")
    return Fraction(int(num), int(den))


def main(path: str | None = None) -> int:
    cert_path = pathlib.Path(
        path
        or SCRIPTS.parent / "certificates" / "collar_poisson_witness_certificate.json"
    )
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    recomputed = witness.compute()
    rec_families = {f["family_id"]: f for f in recomputed["families"]}

    checks: dict[str, bool] = {
        "issue_is_320": cert.get("issue") == 320,
        "certificate_id_v1": (
            cert.get("certificate_id") == "issue-320-collar-poisson-witness-v1"
        ),
        "all_families_flag_true": cert.get("all_families_within_bound") is True,
        "family_count_matches": (
            len(cert.get("families", [])) == len(recomputed["families"])
        ),
        "counting_step_theorem_grade": (
            cert.get("claim_boundary", {}).get("counting_step")
            == "theorem grade on the declared model class"
        ),
        "physical_realization_stays_conditional": (
            "conditional"
            in cert.get("claim_boundary", {}).get("physical_realization", "")
        ),
        "lambda_collar_not_certified_here": (
            cert.get("claim_boundary", {}).get("lambda_collar_value")
            == "branch data, not certified here"
        ),
        "paper_surface_named": (
            "thm:collar-poisson"
            in cert.get("theorem", {}).get("paper_surface", "")
        ),
        "simulator_proxy_named": (
            "STATIC_GALAXY_RAR_BTFR_RECEIPT"
            in cert.get("consumers", {}).get("simulator_galaxy_proxy", "")
        ),
        "sparc_consumer_named": (
            "sparc_rar_compare.py"
            in cert.get("consumers", {}).get("sparc_comparison", "")
        ),
    }

    for fam in cert.get("families", []):
        fid = fam.get("family_id", "unnamed")
        rec = rec_families.get(fid)
        if rec is None:
            checks[f"{fid}:known_family"] = False
            continue
        stored_tv = Decimal(fam["exact_tv"])
        rec_tv = Decimal(rec["exact_tv"])
        stored_bound = parse_frac(fam["bound_exact"])
        rec_bound = parse_frac(rec["bound_exact"])
        checks[f"{fid}:tv_matches_recomputation"] = (
            abs(stored_tv - rec_tv) <= RECOMPARE_TOLERANCE
        )
        checks[f"{fid}:bound_matches_recomputation"] = stored_bound == rec_bound
        checks[f"{fid}:tv_le_bound"] = (
            fam.get("tv_le_bound") is True
            and stored_tv + witness.NUMERIC_TOLERANCE
            <= witness.frac_to_decimal(stored_bound)
        )

    failures = [name for name, ok in checks.items() if not ok]
    print(f"checks: {len(checks)}, failures: {len(failures)}")
    for name in failures:
        print(f"FAIL {name}")
    if not failures:
        print("collar Poisson witness certificate: VALID")
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1] if len(sys.argv) > 1 else None))
