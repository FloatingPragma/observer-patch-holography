"""Literal C/E geometry and the exact scope of phase discrimination."""
import json
import pytest
import sympy as s
from . import phase_selection as a


@pytest.fixture(scope='module')
def result():
    return a.certificate()


def test_actual_CE_order_is_only_orientation_reversal(result):
    row=result['source_halfturns']
    assert row['actual_halfturn_group_order']==6
    assert row['actual_carrier_rotation_count']==60
    assert row['ordered_EC_is_inverse_of_ordered_CE']
    assert row['unordered_action_spectrum_preserved_by_orientation_reversal']


def test_literal_permutations_conjugate_the_cycle_without_changing_weights():
    d=a.native_geometry();c,e=d['gC'],d['gE']
    ec=tuple(e[c[i]] for i in range(12));ce=tuple(c[e[i]] for i in range(12))
    assert all(ec[ce[i]]==i for i in range(12))
    assert all(c[ec[c[i]]]==ce[i] for i in range(12))


def test_exact_two_distinct_phase_controls_have_same_reversed_spectrum():
    d=a.native_geometry();r=d['R'];reflection=d['C']
    for z in (s.Integer(1),(s.Integer(99)+20*s.I)/101):
        c=s.sqrt(2)*s.eye(3)+z*r+s.conjugate(z)*r*r
        assert a._zero(reflection*c*reflection-c.conjugate())
        assert a._zero(reflection*c*c*reflection-c.conjugate()*c.conjugate())


def test_reflection_forces_real_circulant_and_A5_scalar(result):
    row=result['source_halfturns']
    assert row['halfturn_invariant_Hermitian_commutant_real_dimension']==2
    assert row['reflection_invariant_balanced_positive_response_has_repeated_pair']
    assert row['full_A5_invariant_triplet_response_scalar']
    assert not row['full_A5_invariance_assumed_for_actual_reset_memory']


def test_duty_change_is_not_orientation_change(result):
    row=result['complete_base_memory_duties']
    assert row['fine_memory_rows']==1536
    assert [r['rank'] for r in row['duties']]==[3,3,2,2]
    assert [r['trace'] for r in row['duties']]==['1','1/12','1/6','1/72']
    assert not row['exchanging_all_C_E_record_roles_preserves_source_law']
    assert not row['duties'][0]['commutes_with_actual_C3']


def test_real_history_grade_does_not_erase_order_or_invent_unitarity(result):
    row=result['ordered_feedback']
    assert row['actual_real_projectors_checked']==12
    assert row['noncommuting_ordered_pulse_control_retained']
    assert row['Z_conjugation_Z_fixes_every_pulse_and_word']
    assert not row['finite_step_reduced_channel_is_unitary']
    assert not row['all_family_conjugates_or_mass_spectra_excluded']


def test_fixed_program_and_changed_control_are_distinct(result):
    assert not result['fixed_coherent_program_default']['readable_program_is_ambiguous']
    assert result['fixed_coherent_program_default']['action_ratios_to_smallest']==['1','1','34+24*sqrt(2)']
    assert not result['changed_coherent_control_is_same_complete_program']
    assert not result['empirical_inputs_used']
    assert not result['universal_phase_selection_impossibility_claimed']
    json.dumps(result)


def test_actual_all_lag_memory_has_identity_graded_potential_and_diffusion(result):
    row=result['actual_all_lag_memory']
    assert row['group_inverse_exact'] and row['memory_dimension']==12
    assert row['principal_minors_checked']==15
    assert row['diffusion_positive_semidefinite'] and row['diffusion_nonzero']
    assert row['potential_placement']=='I2 tensor A'
    assert s.simplify(s.sympify(row['potential_family_matrix'][0][1])-603*s.I*(3+s.sqrt(5))/6666724)==0
    assert not row['fixed_gain_fixed_t_ordering_survives_h_to_zero']
    assert not row['retained_potential_without_diffusion_is_same_source']


def test_fourier_label_shift_is_exact_root_permutation():
    omega=-s.Rational(1,2)+s.I*s.sqrt(3)/2
    p=a.coherent.P;z=(s.Integer(99)+20*s.I)/101
    q=s.diag(1,omega,omega**2)
    c=s.sqrt(2)*s.eye(3)+z*p+s.conjugate(z)*p*p
    for m in range(3):
        shifted=z*omega**(-m)
        cm=s.sqrt(2)*s.eye(3)+shifted*p+s.conjugate(shifted)*p*p
        assert a._zero((q.H**m)*c*(q**m)-cm)
        assert a._zero((q.H**m)*c*c*(q**m)-cm*cm)


def test_geometry_result_mapping_mutation_is_detached():
    value=a.native_geometry();value['group_order']=99
    assert a.native_geometry()['group_order']==6


def test_mcpr_fixed_point_and_positive_deformation_are_not_balanced_default(result):
    row=result['conditional_MCPR_phase_law']
    assert row['exact_affine_fixed_points_checked']
    assert row['positive_root_quotient_additional_premise']=='all three response roots r_j are nonnegative'
    assert row['phase_at_least_two_ninths_under_sign_premises']
    assert row['positive_root_quotient_strictly_below_two_thirds_for_positive_u']
    assert row['configured_coordinates_and_architecture_determine_phase']
    assert row['architecture_historically_target_informed']
    assert not row['equals_exact_balanced_Koide_action']
    assert not row['retained_coordinate_receipts_consumed']
    assert not row['parent_memory_or_pixel_closure_replayed']


def test_mcpr_declared_recurrence_direct_rational_controls():
    u=s.Rational(1,25);alpha_p=s.Rational(7,1000)
    chi=(u*u/512)/(1-u/77)
    zeta=(alpha_p*alpha_p/21)/(1-u/27-u*u/135)
    assert chi>0 and zeta>0
    assert u*u/512+u*chi/77==chi
    assert alpha_p*alpha_p/21+(u/27+u*u/135)*zeta==zeta


def test_actual_all_lag_potential_is_positive_cyclic_orientation_result(result):
    row=result['actual_all_lag_memory']
    A=s.Matrix([[s.sympify(x) for x in line] for line in row['potential_family_matrix']])
    R=s.Matrix([[s.sympify(x) for x in line] for line in result['source_halfturns']['native_R']])
    kappa=s.Rational(603,6666724)*(1+s.sqrt(5))
    assert a._zero(A*R-R*A)
    assert a._zero(A-s.I*kappa*(R-R.T))
    assert row['commutes_with_actual_C3']
    assert not row['orientation_generator_supplies_scalar_phase_without_control']
    assert not result['complete_base_memory_duties']['duties'][0]['commutes_with_actual_C3']
