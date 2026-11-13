from fractions import Fraction as F
from copy import deepcopy
import json
import pytest
from . import source as s

def test_conservation_derivative_and_uniform_past_bound():
    result=s.run()
    assert result['uniform_past_certificate']['global_past_upper_exact']=='-1/5'
    assert all(F(row['adiabatic_sound_speed_squared_interval'][1])<F(-1,5)
               for row in result['retained_background_samples'])
    assert result['required_completion']['background_conservation_does_not_determine_rest_sound_speed']

def test_invalid_background_bound_cannot_be_accepted():
    parent=json.loads((s.ROOT/'model/openwave_radiation_capacity/results.json').read_text())
    changed=deepcopy(parent);changed['inputs']['Omega_dark0']='9/10'
    with pytest.raises(AssertionError):s.past_bound(changed)
    changed=deepcopy(parent)
    next(x for x in changed['future_event_branch']if x['redshift']==1)['Omega_dark_interval'][1]='1/2'
    with pytest.raises(AssertionError):s.past_bound(changed)

def test_complete_sound_spectrum_and_stable_countercontrol():
    modes=s.periodic_sound_modes()
    assert sum(x['multiplicity']for x in modes)==8
    assert sum(x['multiplicity']for x in modes if x['Laplacian_eigenvalue_exact']=='0')==1
    for r in (F(0),F(1,3),F(1)):
        hubble=4*r/(3*(3+r));assert 0<=hubble<=F(1,3)

def test_canonical_reconstruction_has_exact_positive_potential_curvature():
    result=s.canonical_completion()
    assert result['V_second_derivative_over_H_squared_bounds_exact']==['2','6']
    assert result['rest_frame_sound_speed_squared_exact']=='1'
    assert not result['source_scalar_transfer_proved']
