"""Operational phase readout, complete histories and its interpretation boundary."""
from fractions import Fraction as F
import sympy as s

from model.coherent_defect_wave.model import JointObserver, combine
from .readout import source_experiment, fourier_branches, filter_clock_branches, powers
from .readout import energy_branch_control, numerical_free_phase_bound, resolution_plan
from .certificate import certificate
from .background_readout import background_experiment


def test_every_original_preparation_filter_and_clock_outcome_is_retained():
    result = certificate()
    rows = result['complete_original_trial_records']
    assert len(rows) == 9
    assert sum(F(row['probability']) for row in rows) == 1
    assert rows[0]['probability'] == '1/2'
    assert sum(F(row['probability']) for row in rows
               if row.get('filter') == 'nonflat_filter') == F(17, 288)
    assert all(row['positive_constraint_preserved'] for row in rows[1:])
    assert result['scope']['physical_energy_branch_selected'] is False


def test_compiled_polynomial_matches_sequential_filter_then_clock():
    observer, _, _, source_powers, direct, _ = source_experiment()
    selected = combine(source_powers[0], source_powers[2], F(1, 2), F(-1, 2))
    sequential = fourier_branches(powers(observer, selected, maximum=3))
    for (_, a), (_, b) in zip(direct[4:], sequential):
        assert a == b
        assert a.gauss_positive()
    assert sum(state.norm_squared() for _, state in sequential) == F(17, 288)


def test_identity_clock_cannot_manufacture_the_active_phase_response():
    observer = JointObserver()
    state = observer.embedded_port(1, 0)
    rows = filter_clock_branches((state,)*6)
    assert [branch.norm_squared() for _, branch in rows] == [1, 0, 0, 0, 0, 0, 0, 0]


def test_nonflat_phase_bound_and_extended_clock_plan_retain_scope():
    *_, moments = source_experiment()
    bound = numerical_free_phase_bound(moments)
    assert bound['spectral_weight_lower_given_positive_Gauss_preparation'] == '8/135'
    assert bound['original_trial_weight_lower_after_Gauss_and_nonflat_filter'] == '1/36'
    assert bound['conditional_nonflat_spectral_weight_lower'] == '8/17'
    plan = resolution_plan()
    assert plan['clock_size'] == 1024
    assert plan['maximum_filter_plus_clock_U_power'] == 1025
    assert plan['complete_clock_tail_upper'] == '17/512'


def test_identical_tick_does_not_determine_a_positive_energy_gap():
    result = energy_branch_control()
    assert result['different_positive_gap_ratio'] == '5'
    H0, H1 = s.diag(0, s.pi/2), s.diag(0, 5*s.pi/2)
    assert H0 != H1
    for tick in (0, 1, 2, 5):
        assert (-s.I*tick*H0).exp() == (-s.I*tick*H1).exp()


def test_stationary_background_face_disturbances_have_complete_active_readout():
    result = background_experiment()
    rows = result['complete_original_background_trial_records']
    assert len(rows) == 40
    assert sum(F(row['probability']) for row in rows) == 1
    active = [row for row in rows if row['filter'] == 'nonflat_filter']
    assert sum(F(row['probability']) for row in active) == F(283,7200)
    assert all(row['orthogonal_to_stationary_reference'] for row in active)
    assert all(row['positive_constraint_preserved'] for row in rows)
    assert result['scope']['physical_ground_state_or_energy_gap_identified'] is False
    control = result['same_source_positive_unique_ground_branch_control']
    assert control['exact_gaps_in_pi_units'] == ['2','4']
    assert control['same_unique_ground_vector'] == 'Omega0'
    assert control['nonlocal_countercontrols_not_selected_physical_generators'] is True
