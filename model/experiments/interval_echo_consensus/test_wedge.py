import json
from fractions import Fraction

from . import __main__ as wedge
from . import phase
from . import collar_reference


def test_two_unit_channels_and_balanced_input_control():
    result = wedge.run()
    saved = json.loads((wedge.HERE / 'results.json').read_text())
    assert json.loads(json.dumps(result)) == saved
    assert [row['middle_double_order_trace_distance'] for row in result['priority_rows']] == [1, 0, 0]
    assert [r['initial'] for r in result['legacy_ready_on_balanced_blank']] == [[1, 1, 0], [0, 1, 1]]
    assert all(r['opposite_phases_with_transcript_distance'] == 0 for r in result['priority_rows'])


def test_addressed_window_preserves_within_winner_phase():
    result = phase.run()
    assert result == json.loads((wedge.HERE / 'phase.json').read_text())
    for row in result['priority_rows']:
        assert row['priority_probe_rail_probabilities'] == [0, 1]
        assert row['both_probes_rail_probabilities'] == [.5, .5]


def test_local_reference_forecasts_held_magnitude():
    observed = collar_reference.measure()
    assert observed["public_records"]["magnitude_two_ratio"] == Fraction(132651, 20000)
    assert observed["equal_parity_splitting_matches"]
    assert observed["all_walk_lengths_ratio_floor"] == Fraction(3456, 625)
    assert json.loads(json.dumps(observed, default=str)) == json.loads(
        (wedge.HERE / "collar_reference.json").read_text()
    )


def test_same_offer_unit_links_energy_to_physical_incidence():
    from .modular_standard import compile_native, energies, local_spectrum, selected_offers, standard_count
    from .offer_service import finite_service, reversible_diagnostics

    m = standard_count()
    depth = selected_offers(m, "common")
    common = finite_service(depth, m, reader="common")
    fair = finite_service(depth, m, reader="fair")
    biased = finite_service(depth, m, reader="common", word=(1, 1, 2))
    assert (m, depth) == (8, 17)
    assert energies(common)[1] == energies(biased)[1] == m
    assert abs(energies(fair)[1] - m / 2) < 1e-12
    assert all(abs(x-1) < Fraction(1,10**11) for x in common["held_ratios"])
    assert all(x > 2 for x in biased["held_ratios"])
    masked = finite_service(depth, m-1, reader="common")
    assert abs(local_spectrum(compile_native(masked))[0]["energy"]-(m-1)) < 1e-12
    capped = reversible_diagnostics(depth, m, cap=512)
    assert capped["modes"] == common["modes"]
    assert 0 < capped["gap"] < 1 and capped["detailed_balance_error"] < 1e-14


def test_stationarity_calibration_recovers_complete_existing_seam():
    import numpy as np
    from model.candidate_laws.contrast_memory_spherical_v1 import law
    from .modular_standard import stationary_operator

    recovered, source = stationary_operator()
    exchange, contrast = law.local_term_parts()
    np.testing.assert_allclose(recovered, exchange+8*contrast, rtol=0, atol=2e-13)
    assert source["native_energy_levels"] == [0,8,9,32,33]


def test_modular_receipt_preserves_rare_probabilities_and_resolved_residuals():
    import pytest
    from .modular_standard import check, quadratic_energy_residual, stable
    from .offer_service import finite_service

    common = finite_service(17, 8)
    residual = quadratic_energy_residual(common)
    assert -7.67e-13 < residual < -7.65e-13
    resolved_ratio = stable({"held_even_ratio": float(common["held_ratios"][0])})
    assert resolved_ratio["held_even_ratio"] > 1
    with pytest.raises(ValueError, match="Exact finite ratio changed"):
        check(resolved_ratio, {"held_even_ratio": 1.0})

    for expected, mutated in (
        ({"one_batch_branches": {"accepted": 6.6e-26}}, {"one_batch_branches": {"accepted": 0.0}}),
        ({"second_energy_minus_quadratic": residual}, {"second_energy_minus_quadratic": 0.0}),
        ({"held_odd_relative_error": 2.65e-13}, {"held_odd_relative_error": 0.0}),
    ):
        with pytest.raises(AssertionError):
            check(expected, mutated)
        check(expected, expected)
    check({"full_operator_max_error": 1e-14}, {"full_operator_max_error": 0.0})
