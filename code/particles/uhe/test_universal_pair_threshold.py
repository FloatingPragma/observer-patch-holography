import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

import pytest


MODULE_PATH = Path(__file__).with_name("universal_pair_threshold.py")
RUNTIME_RECEIPT_PATH = Path(__file__).with_name("runtime") / "universal_pair_threshold.json"
SPEC = importlib.util.spec_from_file_location("universal_pair_threshold", MODULE_PATH)
assert SPEC and SPEC.loader
mod = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = mod
SPEC.loader.exec_module(mod)


def planck_delta() -> float:
    a_evinv = mod.PLANCK_LENGTH_M / mod.HBAR_C_EV_M
    return mod.delta_abs_from_a(a_evinv)


def test_planck_scale_transition_and_global_minimum() -> None:
    delta = planck_delta()
    kt = mod.transition_energy(mod.ELECTRON_MASS_EV, delta)
    kmin, epsmin = mod.global_minimum(mod.ELECTRON_MASS_EV, delta)
    assert kt == pytest.approx(2.538373e17, rel=3e-6)
    assert kmin == pytest.approx(1.928746e17, rel=3e-6)
    assert epsmin == pytest.approx(1.805109e-6, rel=3e-6)
    assert kmin < kt


def test_representative_radio_energy_has_no_window() -> None:
    result = mod.kinematic_window(4e-8, mod.ELECTRON_MASS_EV, planck_delta())
    assert result["status"] == "no_leading_kinematic_window"
    assert result["lower_root_ev"] is None
    assert result["upper_root_ev"] is None


@pytest.mark.parametrize(
    ("epsilon", "lower", "upper"),
    [
        (6.34e-4, 4.11861e14, 7.82222e19),
        (3.0e-3, 8.703998e13, 3.70137e20),
    ],
)
def test_planck_scale_finite_windows(epsilon: float, lower: float, upper: float) -> None:
    result = mod.kinematic_window(epsilon, mod.ELECTRON_MASS_EV, planck_delta())
    assert result["status"] == "finite_leading_kinematic_window"
    assert result["lower_root_ev"] == pytest.approx(lower, rel=5e-6)
    assert result["upper_root_ev"] == pytest.approx(upper, rel=5e-6)
    assert mod.universal_threshold_envelope(
        float(result["upper_root_ev"]), mod.ELECTRON_MASS_EV, planck_delta()
    ) == pytest.approx(epsilon, rel=2e-12)


def test_root_below_old_scan_is_not_mislabeled_absent() -> None:
    delta = mod.delta_abs_from_a(
        1e4 * mod.PLANCK_LENGTH_M / mod.HBAR_C_EV_M
    )
    result = mod.kinematic_window(6.34e-4, mod.ELECTRON_MASS_EV, delta)
    assert result["status"] == "finite_leading_kinematic_window"
    assert float(result["upper_root_ev"]) == pytest.approx(7.82222e15, rel=5e-6)


def test_interior_minimizer_saturates_am_gm_envelope() -> None:
    delta = planck_delta()
    k = 1e20
    u = mod.minimizing_u(k, mod.ELECTRON_MASS_EV, delta)
    explicit = mod.universal_threshold_at_share(k, mod.ELECTRON_MASS_EV, delta, u)
    closed = mod.unconstrained_am_gm_lower_bound(
        k, mod.ELECTRON_MASS_EV, delta
    )
    assert explicit == pytest.approx(closed, rel=2e-15)


def test_constrained_minimum_is_not_am_gm_bound_below_transition() -> None:
    delta = planck_delta()
    kt = mod.transition_energy(mod.ELECTRON_MASS_EV, delta)
    below = kt / 2.0
    above = kt * 2.0

    assert mod.minimizing_u(below, mod.ELECTRON_MASS_EV, delta) == 0.25
    assert mod.universal_threshold_envelope(
        below, mod.ELECTRON_MASS_EV, delta
    ) > mod.unconstrained_am_gm_lower_bound(
        below, mod.ELECTRON_MASS_EV, delta
    )
    assert mod.universal_threshold_envelope(
        above, mod.ELECTRON_MASS_EV, delta
    ) == pytest.approx(
        mod.unconstrained_am_gm_lower_bound(
            above, mod.ELECTRON_MASS_EV, delta
        ),
        rel=2e-15,
    )


def test_invalid_share_fails_closed() -> None:
    with pytest.raises(ValueError, match="finite and lie"):
        mod.universal_threshold_at_share(1e19, mod.ELECTRON_MASS_EV, planck_delta(), 0)


def test_nonfinite_inputs_fail_closed() -> None:
    with pytest.raises(ValueError, match="positive and finite"):
        mod.delta_abs_from_a(math.nan)
    with pytest.raises(ValueError, match="finite and lie"):
        mod.universal_threshold_at_share(
            1e19, mod.ELECTRON_MASS_EV, planck_delta(), math.nan
        )


def test_photon_only_planck_tail_root_corrects_grid_estimate() -> None:
    delta = planck_delta()
    result = mod.photon_only_window(3.0e-3, mod.ELECTRON_MASS_EV, delta)
    assert result["status"] == "finite_leading_kinematic_window"
    upper = float(result["upper_root_ev"])
    assert upper == pytest.approx(3.29497e18, rel=2e-6)
    assert mod.photon_only_threshold(
        upper, mod.ELECTRON_MASS_EV, delta
    ) == pytest.approx(3.0e-3, rel=2e-12)


def test_receipt_keeps_epistemic_firewall() -> None:
    receipt = mod.build_receipt()
    status = receipt["epistemic_status"]
    assert status["frozen_prediction_score"] is False
    assert status["auger_exclusion_recomputed_for_universal_sector"] is False
    assert status["opacity_or_flux_prediction"] is False


def test_receipt_binds_producer_and_lean_source() -> None:
    receipt = mod.build_receipt()
    binding = receipt["source_binding"]
    assert binding["producer_path"] == "code/particles/uhe/universal_pair_threshold.py"
    assert binding["lean_source_path"] == (
        "Lean/Screen/SeamCurrentPhotonLeptonThreshold.lean"
    )
    assert binding["producer_sha256"] == hashlib.sha256(
        mod.PRODUCER_PATH.read_bytes()
    ).hexdigest()
    assert binding["lean_source_sha256"] == hashlib.sha256(
        mod.LEAN_SOURCE_PATH.read_bytes()
    ).hexdigest()


def test_runtime_receipt_is_bound_to_current_sources() -> None:
    receipt = json.loads(RUNTIME_RECEIPT_PATH.read_text(encoding="utf-8"))
    binding = receipt["source_binding"]
    assert binding["producer_sha256"] == hashlib.sha256(
        mod.PRODUCER_PATH.read_bytes()
    ).hexdigest()
    assert binding["lean_source_sha256"] == hashlib.sha256(
        mod.LEAN_SOURCE_PATH.read_bytes()
    ).hexdigest()


def test_external_bound_translations_keep_sector_assumptions() -> None:
    bounds = mod.build_receipt()["conditional_external_bounds"]
    lhaaso = bounds["LHAASO_GRB221009A_quadratic_photon_TOF"]
    assert lhaaso["input_sub_luminal_EQG2_lower_bound_ev"] == 6.9e20
    assert lhaaso["a_max_over_lP"] == pytest.approx(7.9130242e7, rel=2e-8)
    assert "95% subluminal" in lhaaso["scope"]

    auger = bounds["Auger_photon_only_translation_not_universal"]
    assert auger["a_max_over_lP"] == pytest.approx(0.54599867, rel=2e-8)
    assert "subdominant proton component" in auger["scope"]
    assert "photon-only LIV" in auger["scope"]
