"""Guard the independent flux trigger, relative phase and detector eligibility."""

import json
from pathlib import Path
import numpy as np

from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (
    tomography,
)

HERE = Path(__file__).resolve().parent


def test_herald_matches_archived_family_and_spin_can_change_fringe():
    report = json.loads((HERE / "results.json").read_text())
    archive = json.loads((HERE.parent / "dirac_read/results.json").read_text())[
        "physical_degree_six_path"
    ]
    for sector in report["sectors"]:
        for j, source in enumerate(archive["modes"]):
            comparison = {
                x["time"]: x["distal_charge_same_character"]
                for x in source["two_positive_energy_mode_source_reads"]
            }
            for row in sector["records"][f"family_{j}_positive"]:
                p = comparison[row["time"]]
                data = row["fringe"]
                np.testing.assert_allclose(
                    data["charge"]["max_over_axes"], [p, 0, 1 - p], atol=2e-12
                )
                np.testing.assert_allclose(
                    data["supplied_spin"]["min_over_axes"], [0, p, 1 - p], atol=2e-12
                )


def test_count_inversion_keeps_complex_coherences():
    # The Y-positive state has rho[0,1] = -i/2. Integer counts here have no sampling error.
    half = tomography.SHOTS // 2
    record = {
        "sharp_and_miss": [half, half, 0],
        "pair_counts": [[[half, half, 0, 0], [2 * half, 0, 0, 0]]],
    }
    expected = np.outer(tomography.INPUTS[-2], tomography.INPUTS[-2].conj())
    np.testing.assert_array_equal(
        tomography.invert(record, [(0, 1)]), expected.round(14)
    )


def test_detection_accepts_transmitted_quantum_information_and_rejects_misrouted_turns():
    records = json.loads((HERE / "circulation.json").read_text())["sectors"]
    counts = json.loads((HERE / "counts.json").read_text())["panels"]
    for model, panel in zip(records, counts, strict=True):
        assert model["sign"] == panel["sign"] and model["source_spur"] == panel["spur"]
        assert model["rank_preserving_half_turns"] == 0
        assert sum(n > 0 for n in model["closed_address_subspace_sizes"]) == 10
        assert not any(model["order_disagreements_on_closed_address_subspace"])
        read = model["panels"][0]["readings"][2]["source_qubit_to_receiver"]
        assert (
            abs(
                read["minimum_reference_partial_transpose_eigenvalue"]
                - panel["minimum_partial_transpose_eigenvalue"]
            )
            < panel["simultaneous_eigenvalue_error"]
        )
        assert panel["certifies_entanglement"] == bool(model["source_spur"])
        for src in panel["counts"]:
            assert sum(src["sharp_and_miss"]) == tomography.SHOTS
            assert all(
                sum(s) == tomography.SHOTS
                for settings in src["pair_counts"]
                for s in settings
            )


def test_received_winding_distinguishes_quantum_transport_from_spin():
    from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (
        winding_turn,
    )

    report = winding_turn.run()
    for sector in report["sectors"]:
        guide = sector["staged_positive_guide"]
        np.testing.assert_allclose(
            guide["bare_positive_rest_energy"],
            report["family_mass_effects"][0],
            atol=2e-14,
        )
        assert len(guide["guide_spectrum_max_errors"]) == 4
        assert all(
            r["face_captured"] < 1e-25
            for r in sector["staged_cut_after_preparation"]["reads"]
        )
        certified = [
            r["time"]
            for r in sector["positive_branch_face_tomography"]
            if r["certifies_entanglement"]
        ]
        assert certified == [4, 8]
        read = next(r for r in guide["reads"] if r["time"] == 8)
        one, three = read["interior_apertures"]
        assert (
            one["sharp_vertex_records"] == 1
            and three["reference_partial_transpose_min"] < 0
        )
        probe = next(
            r for r in sector["positive_branch_rotation_readings"] if r["time"] == 8
        )
        eligible = [panel for panel in probe["per_axis_pair"] if panel is not None]
        assert len(probe["per_axis_pair"]) == 15 and len(eligible) == 6
        for panel in eligible:
            for setting in panel.values():
                assert sum(setting["counts"]) == sector["rotation_shots_per_setting"]
            p = read["face_captured"]
            np.testing.assert_allclose(
                panel["turn_orders"]["probabilities"][:2], [p, 0], atol=2e-14
            )
            np.testing.assert_allclose(
                panel["supplied_spin"]["probabilities"][:2], [0, p], atol=2e-14
            )


def test_incidence_grade_returns_and_all_rested_endpoint_copies_have_the_same_interior_records():
    from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (
        incidence_turn,
    )

    report = incidence_turn.run()
    assert report["rest_energy_multiplicities_on_one_undressed_guide"] == [1, 1, 1]
    contact = report["single_contact"]["light_character"]
    assert contact["incidence_probabilities_at_laplacian_swap_wait"][2] > 0.4
    assert contact["incidence_return_probabilities_from_sharp_endpoint"][2] < 1e-28
    for sector in report["sectors"]:
        largest_region = sector["interior_rest_measurements"][-1]
        assert largest_region["scalar_matrix_units_checked"] == 100
        assert (
            largest_region["rest_copies"]
            == largest_region["exterior_words_at_capture"]
            == 4
        )
        for axes in sector["axes"]:
            light = axes["variants"]["light_character"]
            assert 1e-8 < light["destructive_given_capture"] < 2e-6
            assert (
                light["destructive_bound_for_every_state_with_commuting_ideal_orders"]
                < 0.00044
            )
            assert (
                axes["variants"]["zero_mass_control"]["destructive_given_capture"] == 0
            )


def test_stationary_family_partners_have_measured_charge_and_negative_frequency():
    from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (
        rest_partner,
    )

    board, route = rest_partner.guide.longest_unpinned_route()
    fixtures = sorted(rest_partner.geometries(board, route))
    assert len(fixtures) == 20
    patch = rest_partner.Patch(board, route, fixtures[0])
    masses = rest_partner.guide.register_effect()[1]
    measured = rest_partner.observer(patch, masses, route)
    for mass, record in zip(masses, measured, strict=True):
        prep = record["heralded_preparation"]
        assert (
            prep["success_from_vertex_source"] > 0.098
            and prep["bright_fidelity_given_success"] > 1 - 1e-10
        )
        assert (
            abs(
                record["relative_energy_read"]["estimated_target_minus_guide_frequency"]
                + 2 * mass
            )
            < 1e-7
        )
        assert record["settings"][-1]["dark_input_to_bright_state"] > 0.9998
        assert (
            max(record["free_charge_reads"][-1]["charge_outside_hub_by_input"]) < 1e-24
        )
        for read in prep["postpreparation_reads"]:
            assert (
                sum(read["counts_failure_rounds_then_click_then_no_click"])
                == prep["shots_per_setting"]
            )
    assert [
        r["positive_mass_level_multiplicity"]
        for r in rest_partner.guide_control(patch, route, masses)["controls"]
    ] == [4, 2, 1]


def test_one_anchor_ratios_reconstruct_from_every_enrolled_detector_count():
    from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (
        anchor_spectrum as spectrum,
    )

    result = json.loads((HERE / "anchor_spectrum.json").read_text())
    primary = result["discovery"]
    for record, expected_rank in (
        (primary, 3),
        (primary["light_filter"], 1),
        (result["scalar_control"], 1),
    ):
        survey = spectrum.unpack(record["survey_counts"])
        blocks = spectrum.unpack(record["refinement_counts"])
        for counts in (survey, blocks):
            assert (counts.sum(axis=-1) == spectrum.SHOTS).all()
        omega, _, _ = spectrum.infer(survey, blocks)
        assert len(omega) == expected_rank
        estimate = (
            record["fresh_frequencies"]
            if "fresh_frequencies" in record
            else record["frequencies_by_acquisition"][0]
        )
        np.testing.assert_allclose(omega, estimate, atol=1e-13, rtol=0)
    assert primary["largest_heldout_quadrature_error"] < primary["heldout_tolerance"]
    assert primary["light_filter"]["acceptance_per_original_enrollment"] < 0.034
    audit = result["scale_audit"]
    assert 9.83 < audit["muon_difference_ppm"] < 9.84
    assert 1.03 < audit["muon_excess_keV_exact_e_anchor"] < 1.05
    assert 1.25 < audit["muon_excess_keV_literal_511_anchor"] < 1.26


def test_finite_clock_visibility_changes_the_three_path_spectrum():
    from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (
        cadence_audit,
    )

    records = cadence_audit.run()["finite_register_comparisons"]
    assert records["phase_only"]["N_fitted_from_muon"] == 302
    assert 149.03 < records["all_records"]["at_N"]["302"]["middle_ratio"] < 149.05
    assert records["all_records"]["N_fitted_from_muon"] == 11959296


def test_scalar_position_impulse_preserves_parked_charge_and_changes_control_flux():
    from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (
        rest_impulse,
    )

    result = rest_impulse.run()
    for read in result["unlabeled_selected_source"]:
        hub, other, bridge, failed = read["counts_hub_elsewhere_outgoing_failed"]
        assert (
            hub + other + bridge + failed == result["shots"] and hub > 0 and failed > 0
        )
        assert other == bridge == 0
    controls = [
        r
        for r in result["acquisitions"]
        if r["character"] == 0 and r["read_wait"] == 0.5
    ]
    assert (
        abs(
            controls[0]["vertex_control_probabilities"][2]
            - controls[1]["vertex_control_probabilities"][2]
        )
        > 0.02
    )


def test_positive_recoil_masses_predict_the_held_mode_from_port_counts():
    from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (
        recoil,
    )

    report = json.loads((HERE / "recoil.json").read_text())
    for name, rank in (("family", 3), ("scalar_control", 1)):
        panel = report[name]
        positive = []
        for row in panel["reads"]:
            survey = recoil.spectrum.unpack(row["survey_counts"])
            blocks = recoil.spectrum.unpack(row["refinement_counts"])
            assert (survey >= 0).all() and survey.max() <= recoil.SHOTS
            assert (blocks >= 0).all() and blocks.max() <= recoil.SHOTS
            omega, _ = recoil.read(survey, blocks, row["normalization"])
            assert len(omega) == 2 * rank
            positive.append(omega[omega > 0][::-1])
        mass, stiffness, held = recoil.fit_first_two(
            len(report["route"]), *positive[:2]
        )
        np.testing.assert_allclose(mass, panel["masses_from_positive_gaps"], atol=1e-12)
        np.testing.assert_allclose(stiffness, 4, atol=1e-8)
        np.testing.assert_allclose(positive[2], held, atol=1e-9, rtol=0)
    control = report["terminal_cost_control"]
    assert control["0.0"]["held_half_rms_complex_fringe_error"] < 0.002
    assert control["1.0"]["held_half_rms_complex_fringe_error"] > 0.4


def test_localized_flight_infers_mass_and_wavelengths_from_arriving_charge():
    from model.candidate_laws.interval_echo.openwave.electron_field.lepton_family.arrival_spin import (
        localized_recoil as study,
    )

    report = json.loads((HERE / "localized_recoil.json").read_text())
    unpack = study.anchor_spectrum.unpack
    n = len(report["route"])

    def records(stored):
        counts = unpack(stored)
        assert np.all(counts.sum(axis=-1) == report["shots_per_basis"])
        return counts

    for row in report["studies"]["family"]["runs"]:
        measured = []
        for name in ("first_flight", "fresh_flight"):
            panel = row[name]
            measured.append(study.read(records(panel["survey_counts"]), records(panel["long_wait_counts"])))
        first, held = measured
        mass, _, predicted = study.recoil.fit_first_two(n, first[:1], first[1:2])
        assert mass[0] > 0 and abs(predicted[0] - held[2]) < 1e-5
        local = records(row["spatial_read"]["survey_counts_site_time_basis_bin"])
        scores = study.spatial_scores(local, first, study.SURVEY[:local.shape[1]])
        np.testing.assert_array_equal(np.argmin(scores, axis=1), [1, 2, 3])
        assert np.max(np.min(scores, axis=1)) < 0.01
        np.testing.assert_array_equal(records(row["release_cut_counts"])[:, :, :2], 0)
    scalar = report["studies"]["scalar"]["same_three_tag_settings"]
    assert np.ptp([row["dispersion_mass"] for row in scalar]) < 2e-7
