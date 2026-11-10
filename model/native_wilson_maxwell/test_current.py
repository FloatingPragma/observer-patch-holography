import sympy as s
from model.native_wilson_maxwell import current as m


def test_transverse_current_does_not_follow_from_scalar_energy():
    x=m.transverse_countercontrol()
    assert x['scalar_square_root_minimal_coupling_vertex']=='0'
    assert x['nonzero_transverse_spin_flip_amplitude']=='3/5'
    assert x['Dirac_density_matrix']==[['4/5','0'],['0','4/5']]


def test_full_Wilson_bond_and_pair_terms():
    x=m.wilson_controls()
    assert x['complete_one_bond_one_body_spectrum']=={'-1/a':2,'0':4,'+1/a':2}
    assert s.Matrix(x['pair_density_nonzero_control'])!=s.zeros(2)
    assert x['canonical_local_charge_spectrum']==[-2,-1,0,1,2]
    assert not x['local_density_equals_original_native_onsite_N_plus_minus_N_minus']
