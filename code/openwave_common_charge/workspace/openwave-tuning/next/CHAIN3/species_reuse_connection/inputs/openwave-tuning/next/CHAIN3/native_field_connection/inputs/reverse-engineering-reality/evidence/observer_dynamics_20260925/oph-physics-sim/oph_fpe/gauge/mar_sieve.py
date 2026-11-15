from __future__ import annotations

from typing import Any

from oph_fpe.claims import CONTINUATION, with_claim_metadata


EXPECTED_SM_GATE = {
    "g_phys": "(SU(3)xSU(2)xU(1))/Z6",
    "hypercharge_lattice": "exact",
    "Nc": 3,
    "Ng": 3,
    "higgs_doublets": 1,
    "light_chiral_exotics": 0,
    "extra_low_scale_u1": 0,
    "xy_gauge_bosons": 0,
}


def standard_model_candidate_sieve(candidate: dict[str, Any]) -> dict[str, Any]:
    checks = {
        "sm_quotient": _normalize_group(candidate.get("g_phys", candidate.get("G_phys"))) == EXPECTED_SM_GATE["g_phys"],
        "exact_hypercharge_lattice": str(candidate.get("hypercharge_lattice", "")).lower() == "exact",
        "Nc_equals_3": _integer_equal(candidate.get("Nc", candidate.get("colors")), 3),
        "Ng_equals_3": _integer_equal(candidate.get("Ng", candidate.get("generations")), 3),
        "one_higgs_package": _integer_equal(
            candidate.get("higgs_doublets", candidate.get("higgs_package_count")), 1
        ),
        "no_light_chiral_exotics": _integer_equal(candidate.get("light_chiral_exotics"), 0),
        "no_extra_low_scale_u1": _integer_equal(candidate.get("extra_low_scale_u1"), 0),
        "no_xy_gauge_bosons": _integer_equal(
            candidate.get("xy_gauge_bosons", candidate.get("X_Y_gauge_bosons")), 0
        ),
    }
    matches_target = all(checks.values())
    report = {
        "mode": "standard_model_target_conformance_diagnostic",
        "SM_TARGET_CONFORMANCE_DIAGNOSTIC": bool(matches_target),
        "SM_QUOTIENT_GATE_RECEIPT": False,
        "PHYSICAL_STANDARD_MODEL_FROM_SCREEN_RECEIPT": False,
        "receipt": False,
        "checks": checks,
        "candidate": dict(candidate),
        "expected_gate": dict(EXPECTED_SM_GATE),
        "claim_boundary": (
            "Target-conformance diagnostic only. Caller-supplied final Standard Model labels cannot "
            "derive an A5 coefficient module, physical port-current algebra, global quotient, matter "
            "selection, family descent, or continuum QFT realization."
        ),
    }
    return with_claim_metadata(report, claim_level=CONTINUATION, receipt="SM_TARGET_CONFORMANCE_DIAGNOSTIC")


def _normalize_group(value: Any) -> str:
    return str(value).replace(" ", "").replace("*", "x")


def _integer_equal(value: Any, expected: int) -> bool:
    return bool(not isinstance(value, bool) and isinstance(value, int) and value == expected)
