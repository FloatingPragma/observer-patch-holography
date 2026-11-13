import numpy as np
from .response import full_native_control,ward_control,finite_response,continuum_reference

def test_complete180_current_matches_native_factor():
    for r in full_native_control():
        assert r['all_negative_positive_transition_columns']==8100
        assert r['absolute_difference']<1e-12

def test_gauge_contact_is_required_and_not_an_inserted_cancellation():
    r=ward_control()
    assert abs(r['total_pure_gauge_energy_response'])<1e-10
    assert abs(r['paramagnetic_second_derivative'])>10

def test_complete_records_and_response_scale_are_not_shape_only():
    r=finite_response(4,.5,[.1,1.,5.],records=True)
    vals=[x['native_subtracted_current_response']for x in r['rows']]
    assert 0<vals[0]<vals[1]<vals[2]
    inst=r['complete_native_current_instrument'];scale=float(inst['current_norm_bound_exact'])
    total=sum(sum(f['all_momentum_spin_summed_pair_weights'])for f in inst['families'])*(8/3)/scale**2
    assert abs(total+inst['failure_probability_diagnostic']-1)<1e-15
    c=continuum_reference(2.)
    assert abs(c['radial_current_integral']-c['independently_integrated_parameter_identity'])<1e-12
