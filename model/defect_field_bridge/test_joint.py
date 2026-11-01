from fractions import Fraction as F
import math

import numpy as np
import pytest

from . import joint
from .dynamics_packet import compare_arrays


def test_complete_source_preparation_retains_charge_and_which_path_reference():
    state, resources = joint.preparation()
    assert state.shape == (16, 8)
    assert resources['all_preparation_addresses'] == 128
    reference = state.conj().T@state
    expected = np.zeros((8, 8))
    expected[4, 4] = expected[5, 5] = .5
    np.testing.assert_allclose(reference, expected, atol=2e-9)
    spectator_charge = np.diag(np.repeat([0., 1., -1., 0.], 2))
    np.testing.assert_allclose(joint.charge()@state+state@spectator_charge, 0, atol=2e-9)


@pytest.mark.parametrize('current', [-1, 0, 1])
def test_all16_source_field_columns_match_independent_full_CAR_product(current):
    actual, expected, record = joint.field_gate(current)
    np.testing.assert_allclose(actual, expected, atol=2e-9, rtol=0)
    np.testing.assert_allclose(actual.conj().T@actual, np.eye(16), atol=2e-9, rtol=0)
    parity = np.diag([(-1)**i.bit_count() for i in range(16)])
    np.testing.assert_allclose(actual@parity, parity@actual, atol=2e-9)
    assert F(record['effective_duration']) == 1+F(current, 4)
    assert record['resources']['native_pulses'] == len(record['source_word'])


def test_corrected_chart_is_needed_for_both_static_and_Higgs_terms():
    from ..defect_field_interactions import interaction
    R = joint.ideal_chart()
    energy = interaction.bare_component_chart((27, 78))['energy']
    for phi in (0., .37, -1.2):
        bare = np.zeros((4, 4), complex)
        bare[3, 0], bare[0, 3] = energy*np.exp(1j*phi), energy*np.exp(-1j*phi)
        np.testing.assert_allclose(R@interaction.pair_hamiltonian(energy, phi)@R.conj().T,
                                   bare, atol=2e-15, rtol=0)
    # Omitting the extra QP phase still diagonalizes the static mass, but
    # fails the driven phase. This control guards the actual joining issue.
    X = np.array([[0, 1], [1, 0]], complex)
    Y = np.array([[0, -1j], [1j, 0]], complex)
    from scipy.linalg import expm
    naive = expm(1j*math.pi/8*(np.kron(X, Y)+np.kron(Y, X)))
    correct = R@interaction.pair_hamiltonian(energy, math.pi/2)@R.conj().T
    wrong = naive@interaction.pair_hamiltonian(energy, math.pi/2)@naive.conj().T
    assert np.linalg.norm(correct-wrong) > energy


@pytest.fixture(scope='module')
def packet():
    return joint.build()


def test_every_native_history_and_clock_outcome_reconstructs_from_complete_factors(packet):
    report, arrays = packet
    for law in report['laws']:
        values = arrays[law+'_all_joint_probabilities']
        assert values.shape == (32, 32, 4)
        assert values.min() >= 0
        assert abs(values.sum()-1) < 2e-8
        for first in range(32):
            for second in range(32):
                field = joint.field_output(arrays, law, first, second)
                sequence = (joint.current(first % 16)+1)*3+joint.current(second % 16)+1
                expected = np.linalg.norm(field)**2*np.sum(abs(arrays['clock_outputs'][sequence])**2,
                                                          axis=(1, 2))
                np.testing.assert_allclose(values[first, second], expected, atol=2e-12, rtol=0)
        assert report['laws'][law]['total_complete_outcomes'] == 4096


def test_full_operator_completeness_reference_preservation_and_created_collisions(packet):
    report, arrays = packet
    original = np.kron(arrays['field_preparation'].conj().T@arrays['field_preparation'],
                       arrays['clock_preparation'].conj().T@arrays['clock_preparation'])
    for law in report['laws']:
        np.testing.assert_allclose(arrays[law+'_joint_effect'], np.eye(64), atol=2e-8, rtol=0)
        np.testing.assert_allclose(arrays[law+'_untouched_reference_marginal'], original,
                                   atol=2e-8, rtol=0)
        row = report['laws'][law]
        assert row['first_collision_probability'] < 1e-14
        assert row['created_second_collision_probability'] > 1e-8
        distribution = arrays[law+'_final_active_number_distribution']
        assert distribution[3] > 1e-8
        assert distribution[[0, 2, 4]].sum() < 1e-14
        if row['policy'] == 'coherent':
            assert row['second_collision_flag_probability'] < 1e-14
        else:
            assert abs(row['second_collision_flag_probability']
                       -row['created_second_collision_probability']) < 2e-8


def test_coherent_policy_and_flagged_control_both_complete_but_physically_distinct(packet):
    report, _ = packet
    assert report['primary_collision_policy'] == 'coherent existing-raw0 hold'
    distances = report['unread_final_field_trace_distance_between_policies']
    assert distances['C'] > 1e-7
    assert distances['E'] < 2e-8  # E already records binary odd endpoint occupancy


def test_clock_inverse_uses_records_and_encloses_each_shared_action_time(packet):
    report, _ = packet
    for row in report['clock_sequences']:
        lo, hi = map(F, row['inferred_time']['effective_time_GeV_inverse_interval'])
        target = F(row['specified_time_for_comparison_only'])
        assert lo <= target <= hi
        assert not row['inferred_time']['inverse_uses_matter_energies_or_program_duration']
        assert F(row['inverse_probability_grid']) == F(1, 10**9)
        assert F(row['inferred_time']['probability_error_bound']) == F(11, 10**9)


def test_quantum_phase_mutation_is_rejected_even_if_all_probabilities_are_unchanged(packet):
    _, arrays = packet
    changed = dict(arrays)
    changed['field_unitaries'] = -arrays['field_unitaries']
    with pytest.raises(ValueError):
        compare_arrays(changed, arrays)


@pytest.mark.parametrize('bad', [-1, 16, '1', True])
def test_invalid_raw_record_is_rejected(bad):
    with pytest.raises(ValueError):
        joint.current(bad)
