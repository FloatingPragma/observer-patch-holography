"""Full-column source Ward identities and exact momentum bounds."""
from fractions import Fraction as F
import math
import pytest
from . import boost_current as b

@pytest.mark.parametrize('axis',range(3))
def test_actual_three_axis_full_fock_identity(axis):
    row=b.source_current(axis)
    assert row['full_Fock_dimension']==64
    assert row['composite_current_Hilbert_Schmidt_norm_squared']=='1/16'
    assert row['all_occupation_sectors_checked']
    assert row['one_particle_centering_control_checked']
    assert not row['physical_boost_or_Lorentz_covariance_derived']

@pytest.mark.parametrize('axis',range(3))
def test_zero_contact_removes_remainder(axis):
    row=b.source_current(axis,F(0))
    assert row['composite_current_Hilbert_Schmidt_norm_squared']=='0'
    assert not row['interaction_and_centering_remainder_nonzero']

@pytest.mark.parametrize('axis',[-1,3,True,1.0])
def test_axis_contract(axis):
    with pytest.raises(ValueError):b.source_current(axis)

@pytest.mark.parametrize('contact',[F(-1),F(1,3),.5])
def test_literal_source_coefficients(contact):
    with pytest.raises(ValueError):b.source_current(0,contact)

def test_actual_sine_translation_error_scales_quadratically():
    rows=[b.alias_bound(F(1,n),F(2),3)for n in(8,16,32)]
    assert rows[0]['free_current_over_speed_squared_momentum_error']==4*rows[1]['free_current_over_speed_squared_momentum_error']
    for n,row in zip((8,16,32),rows):
        for j in range(-20,21):
            k=j/10
            error=3*abs(math.sin(2*k/n)*n/2-k)
            assert error<=float(row['free_current_over_speed_squared_momentum_error'])
        assert not row['interacting_history_momentum_error_bounded']

@pytest.mark.parametrize('a,p,q',[(0,1,1),(1,-1,1),(1,1,-1),(1,1,True),(.5,1,1)])
def test_momentum_scope_contract(a,p,q):
    with pytest.raises(ValueError):b.alias_bound(a,p,q)
