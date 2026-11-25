import json

from . import control, counts, measure, optical_certificate, screen_ward


def bounds(v):
    return json.loads(v)


def test_full_regional_modular_operator_and_boundary_entropy_flow():
    records = measure.run()["charges"]
    for q, sign in (("-1", -1), ("1", 1)):
        r = records[q]
        assert min(bounds(z)[0] for z in r["eigenvalues"]) > 0
        assert bounds(r["work_basis_coefficients"]["strain_work"])[0] > 0.67
        flux = r["flux_rates_nats_per_tick"]
        a, b = bounds(flux["head_onsite"])
        assert sign * a > 0 and sign * b > 0
        a, b = bounds(flux["total_rate"])
        assert sign * a > 0.2 and sign * b > 0.2


def test_equal_traction_at_identical_ranges_and_native_force():
    record = control.run()
    assert (
        record["constraint_rank"] == 22 and record["work_derivative_exact"] == "1/16788"
    )
    for v in record["charges"].values():
        assert bounds(v["finite_entropy_increase_at_0p01_lower"])[0] > 4e-7
    fair = record["fair_untagged_load_region"]
    assert len({tuple(w) for w in fair["load_words"]}) == 5
    assert bounds(fair["finite_entropy_change_at_0p01"])[0] > 4.93893e-7
    assert bounds(fair["signed_entropy_slope"])[1] < 4.95e-5
    observed = counts.run()
    for v in observed["charges"].values():
        assert bounds(v["signed_entropy_difference"])[0] > 0
        assert bounds(v["signed_work_difference"])[0] > 0
    assert observed["fair_untagged"]["source_charge_tag_retained"] is False
    assert bounds(observed["fair_untagged"]["signed_entropy_difference"])[0] > 0


def test_two_hop_material_ward_and_signed_held_focusing():
    numeric_record = screen_ward.run()
    proof_record = optical_certificate.run()
    numeric = numeric_record["charges"]
    proved = proof_record["charges"]
    for name, sgn in [("forward", -1), ("reverse", 1)]:
        a, b = bounds(
            proof_record["fair_untagged_images"][name]["source_focusing_slope"]
        )
        val = numeric_record["fair_untagged_images"][name][
            "material_terminal_ward_read"
        ][3]
        assert a <= val <= b and sgn * a > 0 and sgn * b > 0
    for q, forward, reverse in (("-1", -1, 1), ("1", 1, -1)):
        assert numeric[q]["phase_and_head_settings_per_receiver_max"] == [1, 6, 49]
        for ray, sign in (("forward", forward), ("reverse", reverse)):
            a, b = bounds(proved[q][ray]["source_focusing_slope"])
            assert (
                a
                <= numeric[q]["directional_images"][ray]["material_terminal_ward_read"][
                    3
                ]
                <= b
            )
            assert sign * a > 0 and sign * b > 0


def test_pooled_odd_focusing_at_identical_bridge_energy_current():
    from . import energy_odd, energy_phase

    data = energy_odd.run()
    assert data["fixed_energy_current_derivatives_by_charge"] == [0, 0]
    assert data["equal_work_control"]["constraint_rank"] == 24
    assert data["fixed_work_control"]["constraint_rank"] == 25
    assert data["fixed_work_control"]["work_derivative_exact"] == "0"
    required = []
    for key in (
        "equal_work_control",
        "fixed_work_control",
        "fixed_work_second_control",
    ):
        arm = data[key]
        assert bounds(arm["native_lottery_slack_lower"])[0] > 0.03
        screen = arm["screen"]
        a, b = bounds(screen["odd_forward_minus_reverse"]["source_focusing_slope"])
        assert a * b > 0
        required.append(bounds(screen["required_common_clock_acceleration"]))
    assert required[0][1] < required[1][0] or required[1][1] < required[0][0]
    assert required[2][1] < required[1][0] or required[1][1] < required[2][0]
    records = energy_phase.run(data)
    assert records["source_charge_tag_retained"] is False
    for report in records["controls"].values():
        a, b = report["plus_minus_contrast_99pct"]
        assert a < 0 < b
