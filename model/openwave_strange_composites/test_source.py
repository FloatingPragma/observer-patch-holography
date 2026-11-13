from fractions import Fraction as F
import numpy as np
from .certificate import calculate
from .symmetry import charge_conjugation
from .source import source
from .nonleptonic import word


def test_full_C_map_and_light_operator_transfer():
    d=charge_conjugation()
    assert d['complete_meson_dimension']==256 and d['complete_meson_permutation_error']==0
    assert source('baryon')['light_transfer_error']<1e-10
    assert source('meson')['light_transfer_error']<1e-10


def test_complete_native_current_and_held_out_failures():
    d=calculate();rows={x['name']:x for x in d['mass_rows']}
    assert abs(rows['K+']['computed_GeV']-.493677)<1e-10
    assert rows['K+']['computed_GeV']==rows['K-']['computed_GeV']
    assert abs(rows['Sigma0']['difference_GeV'])>.07
    assert d['zero_Cabibbo_control']['exclusive_proton_probability']==0
    r=d['complete_native_current_instrument']['complete_original_native_occupation_records']
    total=r['preparation_failure']+sum(r['success_and_current_success_occupation'])+sum(r['success_and_current_failure_occupation'])
    assert abs(total-1)<1e-8


def test_nonleptonic_vertex_does_not_override_full_source_threshold():
    d=calculate()['native_nonleptonic_vertex_and_obstruction']
    assert d['native_N5_vertex_norm_squared']>0
    assert F(d['strong_Lambda_to_channel_gap_GeV'])>F('.092')
    assert F(d['gap_after_two_full_weak_norm_shifts_GeV'])>0
    assert not d['strong_channel_is_kinematically_open']
    assert word(0,[(True,0),(True,0)])==(None,0)
