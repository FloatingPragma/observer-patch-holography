import importlib.util
import json
import sys
from decimal import Decimal
from pathlib import Path

import pytest


MODULE_PATH = Path(__file__).with_name("capacity_alpha_tangent.py")
SPEC = importlib.util.spec_from_file_location("capacity_alpha_tangent", MODULE_PATH)
assert SPEC and SPEC.loader
mod = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = mod
SPEC.loader.exec_module(mod)


def test_chain_rule_reproduces_independent_direct_probe() -> None:
    # Values from the deterministic committed-solver h=1e-4 probe.  This test
    # exercises the chain rule without paying the ~45 s solver cost on every
    # unit-test run.
    out = mod.jacobian_from_probe(
        p=Decimal("1.630968209403959324879279848"),
        alpha=Decimal(1) / Decimal("137.035999177"),
        alpha_u=Decimal("0.0411243361955857835710048675537109375"),
        d_alpha_u_d_p=Decimal("-0.02283674315549433231353759766"),
        pi=Decimal("3.1415926535897932384626433832795028841971693993751"),
        sqrt_pi=Decimal("1.7724538509055160272981674833411451827975494561224"),
    )
    assert out["dln_n_dln_alpha"] == pytest.approx(
        Decimal("-0.2101831101348"), rel=Decimal("2e-12")
    )
    assert out["dln_alpha_dln_n"] == pytest.approx(
        Decimal("-4.75775622198"), rel=Decimal("2e-12")
    )


def test_corrected_conservative_public_data_envelopes() -> None:
    diagnostics = mod.observational_diagnostics(Decimal("-0.21018311013484586"))

    clock = diagnostics["Filzinger_2023_clock"]
    assert Decimal(
        clock["conservative_95_percent_abs_rate_envelope_per_year"]
    ) == Decimal("6.70e-19")

    espresso = diagnostics["ESPRESSO_HE0515_4414"]
    espresso_envelope = Decimal(
        espresso["conservative_95_percent_abs_fractional_envelope"]
    )
    assert espresso_envelope == pytest.approx(
        Decimal("3.965888219712146e-6"), rel=Decimal("2e-15")
    )
    # The old 2*sigma-only value (2.72 ppm) omitted the nonzero central value.
    assert espresso_envelope > Decimal("2.8e-6")

    planck = diagnostics["Hart_Chluba_Planck_2018"]
    assert Decimal(planck["central_fractional_change"]) == Decimal("0.0005")
    assert Decimal(planck["one_sigma_fractional_change"]) == Decimal("0.0024")
    assert Decimal(
        planck["conservative_95_percent_abs_fractional_envelope"]
    ) == Decimal("0.005204")


def test_finite_fractional_change_is_converted_to_log_space() -> None:
    envelope = Decimal("0.005204")
    converted = mod.fractional_to_abs_log_envelope(envelope)
    assert converted == pytest.approx(
        -(Decimal(1) - envelope).ln(), rel=Decimal("1e-27")
    )
    assert converted > envelope


def test_large_integrated_capacity_extrapolation_fails_closed() -> None:
    derivative = Decimal("-4.757756221983946")
    small = mod.linearized_delta_ln_alpha_from_capacity(Decimal("0.001"), derivative)
    assert small == Decimal("-0.004757756221983946")

    with pytest.raises(ValueError, match="extrapolation forbidden"):
        mod.linearized_delta_ln_alpha_from_capacity(Decimal("0.01"), derivative)
    with pytest.raises(ValueError, match="exact epoch branch"):
        mod.linearized_delta_ln_alpha_from_capacity(Decimal("0.1"), derivative)


def test_all_three_physical_premises_are_explicit_and_undischarged() -> None:
    assert set(mod.PREMISES) == {
        "B1_same_capacity",
        "B2_epochwise_bridge",
        "B3_physical_solver_tangent",
    }
    assert all(
        str(premise["status"]).startswith("undischarged")
        for premise in mod.PREMISES.values()
    )
    assert (
        "local optical clocks"
        in mod.PREMISES["B3_physical_solver_tangent"]["statement"]
    )


def test_committed_receipt_binds_current_producer_and_solver() -> None:
    receipt_path = (
        MODULE_PATH.with_name("runtime") / "capacity_alpha_tangent_retrospective.json"
    )
    receipt = json.loads(receipt_path.read_text(encoding="utf-8"))
    assert receipt["solver"]["receipt_script_sha256"] == mod._sha256(MODULE_PATH)
    assert receipt["solver"]["paper_math_sha256"] == mod._sha256(mod.PAPER_MATH)
    assert receipt["solver"]["h_sweep_passed"] is True
