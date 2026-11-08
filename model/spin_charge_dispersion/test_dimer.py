from fractions import Fraction as F
import numpy as np
import pytest

from . import dimer


def test_exact_dimer_source_algebra_and_periodic_counts():
    row=dimer.pair_certificate()
    assert row['trial_spin_dot_energy_per_site']=='-31/40'
    assert row['trial_shifted_H_over_J_per_site']=='-61/40'
    for side in (4,6,8):
        torus=dimer.torus_certificate(side)
        assert torus['sites']==side**3
        assert torus['distinct_source_bonds']==3*side**3
        assert torus['inter_dimer_bonds']==5*side**3//2


def test_pair_and_independent_dimer_product_bond_expectations():
    pauli=[np.array([[0,1],[1,0]])/2,
           np.array([[0,-1j],[1j,0]])/2,
           np.diag([1,-1])/2]
    a=(np.sqrt(3)+np.sqrt(2))/np.sqrt(10)
    b=(np.sqrt(3)-np.sqrt(2))/np.sqrt(10)
    pair=np.array([0,a,-b,0])
    h=sum(np.kron(x,x) for x in pauli)
    assert abs(np.vdot(pair,h@pair)+.35)<1e-14
    product=np.kron(pair,pair)
    # The middle sites belong to different dimers and have opposite means.
    inter=sum(np.kron(np.kron(np.eye(2),np.kron(x,x)),np.eye(2)) for x in pauli)
    assert abs(np.vdot(product,inter@product)+.24)<1e-14


def test_exact_current_and_improved_order_margin():
    row=dimer.certificate()
    e=F(row['unshifted_ground_spin_energy_e_lower'])
    assert e/3-F(1,4)==F(row['source_current_coefficient_kappa_lower'])==F(1,120)
    assert F(row['complete_inverse_current_measure_over_J_lower'])==F(1,240)
    lo,hi=map(F,row['sqrt_e_over_18_outward_interval'])
    assert lo*lo<=e/18<=hi*hi
    b=F(row['staggered_Bragg_weight_lower'])
    assert F(row['selected_magnetization_lower'])==18*b/(1+81*b)>F(11,100)
    assert row['source_current_coefficient_is_a_pole_speed'] is False


@pytest.mark.parametrize('side',[2,3,5,0])
def test_periodic_source_outside_valid_counting_domain_rejected(side):
    with pytest.raises(ValueError):
        dimer.torus_certificate(side)
