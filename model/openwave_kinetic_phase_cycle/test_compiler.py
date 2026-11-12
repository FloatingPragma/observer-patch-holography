"""Source sign, complete-carrier and finite-error controls."""
from fractions import Fraction as F
import numpy as np
from scipy.linalg import expm
import pytest
from .compiler import (exact_sign_controls,spatial_inversion_certificate,
    trace_bound,rounds_for_error,compiled_block,exterior_permutation)


def test_exact_complete_sign_control_and_native_real_space_kinetic():
    result=exact_sign_controls()
    assert not result['leg_parity_alone_reverses_total_H']
    assert result['actual_CAR_permutation_checks']==8
    assert spatial_inversion_certificate(3)['spin_spatial_rows']==54


def test_fractional_time_resource_ceiling_and_minimal_rounds():
    for t in (F(0),F(1,3),F(1),F(13,7)):
        n=rounds_for_error(F(1,10**8),t,radius=F(7,3))
        assert trace_bound(n,t,radius=F(7,3))<=F(1,10**8)
        if n>1:
            try:previous=trace_bound(n-1,t,radius=F(7,3))
            except ValueError:continue
            assert previous>F(1,10**8)
    with pytest.raises(ValueError):trace_bound(1,F(1),radius=2)


@pytest.mark.parametrize('kappa',(F(-1),F(1,2),F(1)))
def test_noncommuting_complete_source_evolution(kappa):
    row=compiled_block(F(1,5),kappa,rounds=64,precision=48)
    assert abs(sum(float(x)for x in row['all_final_ancillary_probabilities'])-1)<1e-14
    assert float(row['physical_embedding_Frobenius_error'])<float(F(row['complete_unhalved_trace_bound']))/2


def test_spatial_exterior_sign_matters_for_odd_correlations():
    swap=np.asarray(exterior_permutation((2,3,0,1)),complex)
    assert np.linalg.norm(swap@swap-np.eye(16))==0
    assert swap[5,5]==-1
    with pytest.raises(ValueError):exterior_permutation((0,0))


def test_generic_noncommuting_symmetric_fourth_order_converges():
    # Independent two-row generators test the formula beyond the rank-one
    # native active block and detect mistaken first-order pulse ordering.
    sx=np.array([[0,1],[1,0]],complex);sz=np.diag([1.,-1.]);sy=np.array([[0,-1j],[1j,0]])
    generators=[sx+.2*sz,sz-.3*sx,sy+.1*sx,.4*sx-.7*sz]
    p=1/(4-4**(1/3));target=sum(generators)/4
    def run(n):
        def two(h):
            u=np.eye(2,dtype=complex)
            for j,w in ((0,1/8),(1,1/8),(2,1/8),(3,1/4),(2,1/8),(1,1/8),(0,1/8)):u=expm(-1j*h*w*generators[j])@u
            return u
        a=two(p/n);b=two((1-4*p)/n)
        return np.linalg.norm(np.linalg.matrix_power(a@a@b@a@a,n)-expm(-1j*target))
    assert run(8)/run(16)>14
