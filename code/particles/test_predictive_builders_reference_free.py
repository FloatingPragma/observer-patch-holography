#!/usr/bin/env python3
"""Directory-wide guard against reference injection into particle builders."""

from __future__ import annotations

import pathlib


ROOT = pathlib.Path(__file__).resolve().parents[1]
SCAN_DIRS = [
    ROOT / "particles" / "calibration",
    ROOT / "particles" / "flavor",
    ROOT / "particles" / "leptons",
]

# These modules are explicitly target-fed audits, inverse adapters, or
# compare-only diagnostics. New reference consumers are rejected by default.
REFERENCE_CONSUMER_ALLOWLIST = {
    "calibration/derive_d10_ew_common_transport_gap_report.py",
    "calibration/derive_d10_ew_exactness_audit.py",
    "calibration/derive_d10_ew_target_emitter_candidate.py",
    "calibration/derive_d10_ew_target_free_repair_value_law.py",
    "calibration/derive_d10_ew_tau2_current_carrier_obstruction.py",
    "calibration/derive_d10_ew_w_anchor_neutral_shear_factorization.py",
    "calibration/derive_d11_criticality_comparison.py",
    "calibration/derive_d11_live_exact_higgs_promotion.py",
    "calibration/derive_d11_live_exact_top_promotion.py",
    "calibration/derive_d11_reference_exact_adapter.py",
    "calibration/derive_direct_top_bridge_contract.py",
    "flavor/derive_family_transport_kernel_rho_ord_sensitivity_audit.py",
    "flavor/derive_quark_common_scheme_shape_law_scan.py",
    "flavor/derive_quark_current_family_exact_readout.py",
    "flavor/derive_quark_current_family_exactness_audit.py",
    "flavor/derive_quark_current_family_quadratic_readout_theorem.py",
    "flavor/derive_quark_kernel_normalization_acceptance_harness.py",
    "flavor/derive_quark_shape_law_heldout_prediction_test.py",
    "leptons/audit_charged_stage5_frozen_candidate.py",
    "leptons/derive_charged_d12_continuation_followup.py",
    "leptons/derive_lepton_current_family_exact_readout.py",
    "leptons/derive_lepton_current_family_exactness_audit.py",
    "leptons/derive_lepton_current_family_quadratic_readout_theorem.py",
}
FORBIDDEN_SNIPPETS = ("particle_reference_values.json", "REFERENCE_JSON")


def _production_files() -> list[pathlib.Path]:
    return sorted(
        path
        for directory in SCAN_DIRS
        for path in directory.rglob("*.py")
        if not path.name.startswith("test_") and "__pycache__" not in path.parts
    )


def _relative(path: pathlib.Path) -> str:
    return path.relative_to(ROOT / "particles").as_posix()


def _failures() -> list[str]:
    failures: list[str] = []
    seen_consumers: set[str] = set()
    for path in _production_files():
        text = path.read_text(encoding="utf-8")
        hits = [snippet for snippet in FORBIDDEN_SNIPPETS if snippet in text]
        if not hits:
            continue
        relative = _relative(path)
        seen_consumers.add(relative)
        if relative not in REFERENCE_CONSUMER_ALLOWLIST:
            failures.append(f"{relative}: undeclared reference consumer ({', '.join(hits)})")
    stale = REFERENCE_CONSUMER_ALLOWLIST - seen_consumers
    failures.extend(f"{path}: stale reference-consumer exception" for path in sorted(stale))
    return failures


def test_particle_builder_directories_have_no_undeclared_reference_consumers() -> None:
    assert _failures() == []


def main() -> int:
    failures = _failures()
    if failures:
        print("\n".join(failures))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
