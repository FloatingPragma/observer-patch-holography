from fractions import Fraction as F
import pytest
from .source import shift,character_control,full_uniform_work,nuclear,noncommuting_native_control,oscillator_controls

def test_full_character_trace_and_charge_work():
    assert character_control()['native_suffix_charge_square_mean_exact']=='68/63'
    assert full_uniform_work(3)['bare_Gauss_electric_mean_source_exact']=='8955/9152'
    assert full_uniform_work(4)['bare_Gauss_electric_mean_source_exact']=='2949/1360'

def test_complete_oscillator_square_and_zero_coupling_control():
    assert shift()==F(1,8) and shift(F(0))==0
    assert oscillator_controls()[-1]['distance_to_exact_full_oscillator_ground_diagnostic']<2e-14
    with pytest.raises(ValueError):shift(F(1),F(0))
    row=nuclear(4,2,'alpha',F(0))
    assert row['full_ground_GeV_interval']==row['parent_full_ground_GeV_interval']

def test_triton_protection_does_not_fake_helium_matches_or_static_dynamics():
    assert noncommuting_native_control()['static_elimination_exact_for_full_dynamics'] is False
    t=nuclear(3,1,'triton');h=nuclear(3,2,'helion');a=nuclear(4,2,'alpha')
    assert t['full_ground_GeV_interval']==t['parent_full_ground_GeV_interval']
    assert h['comparison']=='excluded_above' and a['comparison']=='excluded_below'
