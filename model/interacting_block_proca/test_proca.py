"""Sign-sensitive field filters, actual source admission and complete records."""
from fractions import Fraction as F
import sympy as s
import pytest
from . import model


def test_exact_polarization_constraints_and_time_component():
    result=model.algebra.symbolic_identities()
    assert result['Lorenz_polarization_constraint']
    assert model.literal_fourier_pair()['wrong_field_only_A0_fails_Hermiticity']


def test_induced_rotation_and_nontrivial_boost_action():
    result=model.algebra.induced_rotation_and_boost()
    assert result['exact_noncollinear_boost_cocycle_checked']
    assert result['exact_proper_cube_polarization_intertwiners']==24
    reference=model.algebra.retained_reference_rotation()
    assert reference['exact_right_reference_blocks_checked']==24**2
    assert reference['native_reference_dimension_remains']==3072


def test_same_native_source_and_complete_reference():
    parent,admission=model.source()
    assert admission['same_H_selected_actual_ground_reference_and_clock']
    assert admission['actual_native_regular_reference_cover_order']==3072
    assert parent['construction']['measured_inputs']==[]


def test_all_four_instrument_records_keep_positive_leakage():
    result=model.complete_instrument_control()
    for row in result['literal_complete_contraction_control']:
        p=list(map(F,row['four_complete_record_probabilities']))
        assert sum(p)==1 and p[2]+p[3]>0


def test_finite_Weyl_entry_compiler_has_correct_vacuum_and_adjoint():
    eta=(1+s.I)/5
    assert s.simplify(model.finite_displacement_element(0,0,eta)-s.exp(-s.Rational(1,25)))==0
    assert s.simplify(model.finite_displacement_element(1,0,eta)-eta*s.exp(-s.Rational(1,25)))==0
    for n in range(4):
        for m in range(4):
            assert s.simplify(s.conjugate(model.finite_displacement_element(n,m,eta))-model.finite_displacement_element(m,n,-eta))==0


def test_filtered_readout_moments_and_actual_ground_error_are_controlled():
    for a,side in ((F(1,4),64),(F(1,8),128)):
        row=model.resources.finite_source_resources(a,side)
        assert F(row['complete_source_vs_finite_free_record_TV_upper_exact'])<=F(19,20000)
        assert row['actual_replicas_per_macrocell_N']>row['comparison_total_number_cutoff_K']
        assert row['same_H_state_reference_clock_as_scalar_readouts']


def test_covariance_error_and_complete_record_transfer_shrink_quadratically():
    first=model.resources.continuum_bounds(F(1,4))
    second=model.resources.continuum_bounds(F(1,8))
    assert F(first['complete_Weyl_record_spacing_error_upper'])==4*F(second['complete_Weyl_record_spacing_error_upper'])
    assert not first['exact_finite_Lorentz_or_boost_symmetry_claim']


def test_complete_actual_ground_temporal_and_spatial_vector_records():
    for a,side in ((F(1,4),64),(F(1,8),128)):
        row=model.readout.controlled_experiment(a,side)
        source=model.resources.finite_source_resources(a,side,accuracy=F(1,10**9))
        assert F(source['complete_source_vs_finite_free_record_TV_upper_exact'])<=F(row['complete_record_TV_error_upper'])
        assert len(row['complete_sixteen_records'])==16
        assert row['exact_supported_Fourier_modes']==19
        for marginal in row['strictly_positive_actual_sign_marginals']:
            assert F(marginal['actual_minus_sign_probability_interval'][0])>0


def test_full_log_dimension_includes_source_reference_and_four_bits():
    from .capacity import logical_capacity,endpoint
    import mpmath as mp
    row=logical_capacity(1,1,2)
    assert row['fermion_modes']==240 and row['five_state_scalar_registers']==8
    assert row['119_state_link_registers']==24
    assert row['reference_dimension']*row['retained_pointer_and_control_dimension']==49152
    direct=mp.iv.log(2**240*5**8*119**24*49152)
    lo,hi=map(F,row['logHilbertDimension_nats_interval_exact'])
    assert lo<=endpoint(direct,True) and endpoint(direct)<=hi
    display_lo,display_hi=map(F,row['logHilbertDimension_nats_display_interval'])
    assert display_lo<=lo and hi<=display_hi


def test_parent_custody_rejects_mismatched_digest():
    from .certificate import validate_pins
    with pytest.raises(ValueError,match='source custody mismatch'):
        validate_pins({model.PARENT:'0'*64})


def test_transitive_source_closure_is_portable_and_contains_all_new_modules():
    from .certificate import input_hashes
    pins=input_hashes()
    assert len(pins)>351
    assert all(not key.startswith('/') and '..' not in key.split('/') for key in pins)
    for name in ('capacity.py','readout.py','algebra.py','PROOF.md','test_proca.py'):
        assert 'model/interacting_block_proca/'+name in pins
