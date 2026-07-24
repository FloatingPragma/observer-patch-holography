"""Pytest policy for optional legacy surfaces."""

from __future__ import annotations

import os
from pathlib import Path

import pytest


LEGACY_D10_TESTS = {
    "calibration/test_d10_current_carrier_frontier_split.py",
    "calibration/test_d10_ew_common_transport_gap_report.py",
    "calibration/test_d10_ew_exact_closure_beyond_current_carrier.py",
    "calibration/test_d10_ew_exact_mass_pair_chart_current_carrier.py",
    "calibration/test_d10_ew_exactness_audit.py",
    "calibration/test_d10_ew_forward_transmutation_certificate.py",
    "calibration/test_d10_ew_source_transport_pair_artifact.py",
    "calibration/test_d10_ew_source_transport_readout_artifact.py",
    "calibration/test_d10_ew_transport_kernel_artifact.py",
    "calibration/test_d10_gravity_shared_edge_entropy_bridge.py",
    "calibration/test_d10_observable_family_artifact.py",
    "calibration/test_single_p_consistency.py",
    "flavor/test_p_driven_flavor_candidate_motion.py",
    "leptons/test_charged_p_to_affine_anchor_bridge_no_go.py",
    "leptons/test_charged_p_to_affine_anchor_reduction.py",
    "test_compute_current_output_table_runtime_surface.py",
}


def _legacy_d10_enabled() -> bool:
    return os.environ.get("OPH_RUN_LEGACY_D10") == "1"


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line(
        "markers", "slow: long-running physics validation (deselect with -m 'not slow')")


def pytest_collection_modifyitems(config: pytest.Config, items: list[pytest.Item]) -> None:
    if _legacy_d10_enabled():
        return
    # Anchor on this conftest's directory, not config.rootpath: the repo-root
    # pytest.ini makes rootpath the repository root, under which the
    # LEGACY_D10_TESTS entries (relative to code/particles/) would never match
    # and the opt-in gate would be silently inert.
    root = Path(__file__).resolve().parent
    skip_legacy = pytest.mark.skip(
        reason="legacy arXiv D10 helper is opt-in; set OPH_RUN_LEGACY_D10=1 and OPH_LEGACY_PARTICLE_DIR"
    )
    for item in items:
        path = Path(str(item.fspath))
        try:
            rel = path.relative_to(root).as_posix()
        except ValueError:
            continue
        if rel in LEGACY_D10_TESTS:
            item.add_marker(skip_legacy)
