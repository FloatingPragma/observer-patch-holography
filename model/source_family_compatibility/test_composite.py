"""Independent full-ground, retarded-generator and remainder checks."""
from fractions import Fraction as F
import json

import numpy as np
import pytest
from scipy.linalg import expm

from . import composite as c


def _numeric(sites):
    row=c.finite_case(sites)
    data=c.source((sites,1,1),row['c'],row['U'])
    num,den,_=c.ground(data,row['ground_energy']);rho=num/den
    root=np.sqrt(row['ground_degeneracy'])*rho
    k=data['H']/data['scale']-float(row['ground_energy'])*np.eye(len(rho))
    def embed(operators):
        return np.stack([np.concatenate([(a.T@root).ravel(),(a@root).conj().ravel()])
                         for a in operators],axis=1)
    v,w=embed(data['a']),embed(data['b']);size=len(rho)
    def action(columns,kind,value=0):
        modes=columns.shape[1]
        blocks=columns.reshape(2,size,size*modes)
        if kind=='generator':return np.concatenate([k@blocks[0],-k.conj()@blocks[1]]).reshape(-1,modes)
        if kind=='time':
            return np.concatenate([expm(-1j*value*k)@blocks[0],
                                   expm(1j*value*k.conj())@blocks[1]]).reshape(-1,modes)
        if kind=='resolvent':
            return np.concatenate([np.linalg.solve(value*np.eye(size)-k,blocks[0]),
                                   np.linalg.solve(value*np.eye(size)+k.conj(),blocks[1])]).reshape(-1,modes)
        raise ValueError('test action')
    one=data['one']/data['scale'];C=np.array(row['C'],float);alpha=float(row['U']/2)
    block=np.block([[one,alpha*np.eye(len(one))],[alpha*np.eye(len(one)),C]])
    return row,data,rho,v,w,action,block


@pytest.mark.parametrize('sites',[2,3])
def test_each_dressed_family_car_and_mixed_operators_not_independent(sites):
    data=c.source((sites,1,1),F(3,2),F(4))
    identity=np.eye(len(data['H']),dtype=int)
    for i,b in enumerate(data['b']):
        for j,d in enumerate(data['b']):
            assert not np.any(c._anti(b,d))
            np.testing.assert_array_equal(c._anti(b,d.T),identity if i==j else np.zeros_like(identity))
    assert np.any(c._anti(data['a'][0],data['b'][0].T))
    assert np.any(c._anti(data['a'][0],data['b'][1].T))


@pytest.mark.parametrize('sites',[2,3])
def test_full_ground_certificate_matches_independent_complete_diagonalization(sites):
    row,data,rho,*_=_numeric(sites)
    eigenvectors_values,vectors=np.linalg.eigh(data['H']/data['scale'])
    selected=np.abs(eigenvectors_values-float(row['ground_energy']))<1e-12
    assert selected.sum()==row['ground_degeneracy']
    assert abs(eigenvectors_values[0]-float(row['ground_energy']))<2e-14
    projector=vectors[:,selected]@vectors[:,selected].T/selected.sum()
    np.testing.assert_allclose(projector,rho,atol=5e-15,rtol=0)
    for a,b in zip(data['a'],data['b'],strict=True):
        assert abs(np.trace(rho@(a.T@a))-.5)<1e-15
        assert abs(np.trace(rho@(b.T@b))-.5)<1e-15


def test_ground_energy_and_positive_pivot_mutations_fail():
    data=c.source((3,1,1),F(2),F(3,7))
    with pytest.raises(ArithmeticError,match='negative'):
        c.ground(data,F(-44,14))
    with pytest.raises(ArithmeticError,match='nonzero row'):
        c.psd_rank(np.array([[0,1],[1,1]],dtype=int))
    with pytest.raises(ValueError,match='symmetric'):
        c.psd_rank(np.array([[1,1],[0,1]],dtype=int))


@pytest.mark.parametrize('sites',[2,3])
def test_exact_response_block_and_gram_match_actual_retarded_generator(sites):
    row,data,rho,v,w,action,block=_numeric(sites);whole=np.concatenate([v,w],axis=1)
    np.testing.assert_allclose(whole.T@whole,np.eye(whole.shape[1]),atol=2e-15,rtol=0)
    evolved=action(whole,'generator');compressed=whole.T@evolved
    np.testing.assert_allclose(compressed,block,atol=5e-15,rtol=0)
    residual=action(w,'generator')-whole@block[:,len(data['a']):]
    np.testing.assert_allclose(residual.T@residual,np.array(row['residual_Gram'],float),atol=2e-14,rtol=0)
    np.testing.assert_allclose(whole.T@residual,0,atol=5e-15,rtol=0)


def test_three_site_residual_is_exact_positive_and_not_dimer_closure():
    row=c.finite_case(3)
    assert row['residual_squared_trace']==F(542232,13225)
    expected=np.sort([18258/13225,79350/13225,173508/13225]*2)
    np.testing.assert_allclose(np.linalg.eigvalsh(np.array(row['residual_Gram'],float)),expected,atol=2e-15,rtol=0)
    assert row['four_spin_eta_bond_expectations']=={'0,1':F(-179,115),'1,2':F(-179,115)}
    assert not row['mixed_CAR_operators_are_independent']


def test_direct_original_moments_detect_discarded_composite_memory():
    row,data,_,v,w,action,block=_numeric(3)
    vectors=v.copy();moments=[]
    for _ in range(4):
        vectors=action(vectors,'generator');moments.append(v.T@vectors)
    one=data['one']/data['scale'];alpha=float(data['u']/2);C=block[len(one):,len(one):]
    np.testing.assert_allclose(moments[1],one@one+alpha**2*np.eye(len(one)),atol=1e-14,rtol=0)
    np.testing.assert_allclose(moments[2],one@one@one+2*alpha**2*one+alpha**2*C,atol=2e-14,rtol=0)
    # R couples out only after reaching the second block. Its positive norm
    # first changes the original fourth moment by alpha^2 R*R.
    fourth=np.linalg.matrix_power(block,4)[:len(one),:len(one)]
    np.testing.assert_allclose(moments[3]-fourth,alpha**2*np.array(row['residual_Gram'],float),
                               atol=2e-14,rtol=0)
    for _ in range(2):vectors=action(vectors,'generator')
    sixth=v.T@vectors;truncated=np.linalg.matrix_power(block,6)[:len(one),:len(one)]
    assert np.linalg.norm(sixth-truncated)>1e-3


@pytest.mark.parametrize('time',[.03125,.125,.5])
def test_actual_time_error_obeys_original_response_quadratic_bound(time):
    row,data,_,v,w,action,block=_numeric(3);whole=np.concatenate([v,w],axis=1)
    exact=action(v,'time',time)
    approximate=whole@expm(-1j*time*block)[:,:v.shape[1]]
    error=np.linalg.norm(exact-approximate,2)
    bound=np.sqrt(float(row['residual_squared_norm_upper']))*min(time,float(row['U']/2)*time*time/2)
    assert 0<error<bound


@pytest.mark.parametrize('z',[1+2j,-.5+8j,2-3j])
def test_actual_resolvent_error_obeys_volume_independent_second_order_bound(z):
    row,_,_,v,w,action,block=_numeric(3);whole=np.concatenate([v,w],axis=1)
    exact=whole.T@action(whole,'resolvent',z)
    approximate=np.linalg.inv(z*np.eye(len(block))-block)
    assert 0<np.linalg.norm(exact-approximate,2)<float(row['residual_squared_norm_upper'])/abs(z.imag)**3


def test_dimer_composite_time_response_is_exact_on_original_full_ground():
    row,_,_,v,w,action,block=_numeric(2);whole=np.concatenate([v,w],axis=1)
    assert row['residual_squared_norm_upper']==0
    np.testing.assert_allclose(action(whole,'time',.73),whole@expm(-.73j*block),atol=3e-15,rtol=0)


def test_wrong_spin_eta_bond_readout_is_rejected(monkeypatch):
    c.finite_case.cache_clear()
    original=c.bond_observable
    monkeypatch.setattr(c,'bond_observable',lambda *args:np.zeros_like(original(*args)))
    with pytest.raises(ArithmeticError,match='spin/eta'):
        c.finite_case(3)
    c.finite_case.cache_clear()


def test_native_twist_keeps_complete_source_and_export_columns():
    row=c.native_twist()
    assert row['native_pulses']==10 and row['full_input_columns']==4 and row['failure_slots']==36
    assert row['numerical_admission_passed']
    assert F(row['complete_dilation_operator_error_interval'][1])<F(1,10**8)
    assert not row['odd_transfer_or_reservoir_service_derived']


def test_omitted_native_twist_complement_is_rejected(monkeypatch):
    from model.field_thermal import thermal
    original=thermal.execute_with_exports
    def corrupt(*args,**kwargs):
        result,events=original(*args,**kwargs)
        return result,(events[0][:-1],*events[1:])
    monkeypatch.setattr(thermal,'execute_with_exports',corrupt)
    with pytest.raises(ArithmeticError,match='complement'):
        c.native_twist()


def test_exact_certificate_retains_coefficients_and_no_physical_assignment():
    row=c.certificate()
    assert row['independent_CAR_flavours_added']==row['physical_species_added']==0
    assert not row['all_volume_residual_discarded']
    assert row['finite_full_ground_controls'][1]['U']==F(3,7)
    json.dumps(row,default=str,allow_nan=False)
