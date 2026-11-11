from fractions import Fraction as F
import sympy as s

from . import comparison, neutral, self_energy, sterile, sterile_records, transport


def test_full_quantum_sterile_interaction_and_majorana_selection():
    record = sterile.operator_certificate()
    assert record['all_unit_blank_quark_leakage_squared'] == '24'
    assert record['neutral_algebraic_contact_slice']['neutral_four_mode_rank'] == 4
    native = neutral.certificate()
    assert not native['native_Feshbach_elimination_can_generate_DeltaL_two']
    assert native['declared_unit_Majorana_architecture']['nonzero_normal_ordered_pair_terms'] == 6


def test_neutral_particle_and_hole_packets_reach_a_separate_site_in_full_source():
    for row in transport.calculate()['complete_addition_and_removal_transport_controls']:
        assert row['exact_power_supports'] == [2, 23, 1252, 15739]
        assert row['normalized_target_matrix_elements_of_H_powers'] == ['0', '0', 'I/3', '3*I/2']
        low, high = row['receiver_probability_bounds_numeric']
        assert 0 < low < (1e-4)**4/36 < high


def test_native_complete_response_has_normal_positive_variance_complement():
    row = self_energy.calculate()
    assert row['complete_original_addition_removal_records']['complete_joint_record_count'] == 44
    assert row['cyclic_Feshbach_complement_dimension'] == 8
    assert F(row['normalized_second_central_energy_moment']['exact_outward_interval'][0]) > 0


def test_entire_dirac_control_ground_and_exact_active_sterile_eom_spectral_relation():
    result = sterile_records.calculate(False)
    assert result['ground_energy_interval'] == ['-4', '-4']
    assert sum(row['rank'] for row in result['complete_original_preparation_energy_records']) == 80
    kinds = result['complete_neutral_probe_records']
    weights = {}
    for kind, record in kinds.items():
        rows = record['four_settings'][0]['complete_energy_records']
        assert len(rows) == 12
        weights[kind] = {F(row['energy_transfer_interval'][0]): F(row['exact_unscaled_probe_residue'])
                         for row in rows if not row['exact_zero_residue']}
    assert weights['sterile'] == {F(1): F(35, 72), F(5): F(1, 72)}
    assert weights['active_composite'] == {gap: gap*gap*weight for gap, weight in weights['sterile'].items()}
    z = s.Symbol('z')
    function = s.sympify(result['unit_Dirac_control_complete_normal_selfenergy']['complete_CAR_retarded_propagator'])
    assert s.limit(z*function, z, s.oo) == 1


def test_majorana_new_ground_and_all_conjugate_records_preserve_original_probability():
    result = sterile_records.calculate(True)
    assert F(result['ground_energy_interval'][1]) < -4
    assert len(result['complete_original_preparation_energy_records']) == 35
    for probe in result['complete_neutral_probe_records'].values():
        lower = upper = F(0)
        for setting in probe['four_settings']:
            rows = setting['complete_energy_records']
            assert len(rows) == 35 and sum(not row['exact_zero_residue'] for row in rows) == 12
            for row in rows:
                a, b = map(F, row['complete_preparation_joint_probability_interval'])
                lower += a; upper += b
            a, b = map(F, setting['joint_complete_preparation_failure_probability_interval'])
            lower += a; upper += b
        assert lower <= F(1, 80) <= upper
    assert sum(F(row['original_preparation_probability']) for row in result['terminal_preparation_records_before_ground_probe_refinement']) == F(79, 80)


def test_conditional_unit_preserves_gap_ratio_and_identifies_its_empirical_input():
    result = comparison.calculate()
    values = result['comparison_calibration']['conditional_GeV_values_exact']
    assert s.simplify(s.sympify(values[1])/s.sympify(values[0])) == 5
    assert not result['physical_rest_pole_attachment_proved']
    assert not result['comparison_calibration']['active_or_sterile_species_assigned_to_either_value']
