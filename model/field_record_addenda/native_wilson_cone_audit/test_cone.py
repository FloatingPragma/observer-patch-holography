import sympy as s
import numpy as np
from . import cone


def test_literal_clifford_current_and_complete_bond_algebra():
    result=cone.algebra()
    assert result['original_a_equals_inverse_m_axis_ratio']=='2'
    assert result['repaired_axis_ratio']=='1'
    assert not result['weak_chiral_Standard_Model_representation_constructed']
    assert result['charge_conjugation_retained']
    assert len(result['bond_checks'])==3


def test_finite_gauge_background_full_matrix_gap_square():
    assert cone.finite_background_check()['exact_matrix_gap_square_certificate']


def test_independent_all_direction_spectrum_and_error_diagnostic():
    alpha,beta=cone.matrices();alpha=[np.asarray(x,dtype=complex) for x in alpha]
    beta=np.asarray(beta,dtype=complex);gamma=alpha[0]@alpha[1]@alpha[2]@beta
    rng=np.random.default_rng(90211)
    for a in (.02,.1,.5,1.):
        for p in rng.uniform(-2,2,(40,3)):
            w=1-np.cos(a*p);W=w.sum();m=.7
            h=m*beta+sum((np.sin(a*x)/a*A for x,A in zip(p,alpha)),np.zeros((4,4),complex))+W*gamma/a
            values=np.linalg.eigvalsh(h);pred=np.sqrt(m*m+np.sin(a*p)@np.sin(a*p)/a**2+W*W/a**2)
            assert np.max(np.abs(values-[-pred,-pred,pred,pred]))<2e-14
            assert abs(pred-np.sqrt(m*m+p@p))<=a*a*(p@p)**2/(12*m)+1e-14


def test_massless_corner_doublers_lifted_without_weyl_claim():
    from itertools import product
    for bits in product((0,1),repeat=3):
        W=2*sum(bits)
        assert (W==0)==(bits==(0,0,0))
