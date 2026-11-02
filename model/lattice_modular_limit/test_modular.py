from copy import deepcopy
from fractions import Fraction as F
from math import comb
import numpy as np
import pytest
from . import algebra,bounds,controls,packet


def functional(matrix,fun):
    eig,vectors=np.linalg.eigh(matrix)
    return (vectors*fun(eig))@vectors.conj().T


def test_exact_endpoint_and_full_CAR_identities():
    assert algebra.certificate()['lattice_discriminant']=='4096*r**2*(r + 1)**2'
    assert controls.certificate()['relative_entropy_at_theta_pi_over_two']=='3*log(4)/5'


def test_regular_log_series_controls_entire_closed_spectral_interval():
    x=np.linspace(0,1,1001)
    for d in (F(1,2),F(1,8)):
        for N in (4,16):
            r=1/(1+2*d);q=1-2*x
            approx=sum(float(2*r**(2*j+1)/F(2*j+1))*q**(2*j+1) for j in range(N))
            exact=np.log((1-x+float(d))/(x+float(d)))
            tail=2*r**(2*N+1)/(F(2*N+1)*(1-r*r))
            assert np.max(abs(approx-exact))<=float(tail)+1e-14


def test_regularized_noncommuting_matrix_intertwiner_budget():
    C=np.diag([.2,.5,.8]);u=np.array([[.8,-.6,0],[.6,.8,0],[0,0,1.]])
    A=u@np.diag([.2001,.5,.7999])@u.T;f=np.array([1.,0.,0.])
    d=F(1,4);N=8;Y=np.eye(3)-2*C
    residual=max(np.linalg.norm((A-C)@np.linalg.matrix_power(Y,j)@f) for j in range(2*N-1))
    certified=F(int(np.ceil(residual*10**12)),10**12)
    fun=lambda x:np.log((1-x+float(d))/(x+float(d)))
    error=np.linalg.norm(functional(A,fun)@f-functional(C,fun)@f)
    assert error<=float(bounds.regular_log_budget(d,N,certified,map_norm=1))


def test_weighted_phase_Bernstein_bound_including_endpoints():
    x=np.linspace(0,1,2001);N=64
    def phase(y,t):
        out=np.zeros_like(y,dtype=complex);mask=(y>0)&(y<1);z=y[mask]
        out[mask]=z*(1-z)*np.exp(-1j*t*np.log((1-z)/z));return out
    for t in (-2.,0.,.7,2.):
        values=phase(np.arange(N+1)/N,t)
        poly=sum(values[j]*comb(N,j)*x**j*(1-x)**(N-j) for j in range(N+1))
        assert np.max(abs(poly-phase(x,t)))<=(1+abs(t))/(2*np.sqrt(N))


def test_endpoint_atoms_prevent_unweighted_flow_convergence():
    phases=[];weighted=[]
    for n in range(2,9):
        p=1/(1+np.exp(n*np.pi));k=np.log((1-p)/p)
        phases.append(np.exp(-1j*k));weighted.append(p*(1-p)*np.exp(-1j*k))
    assert min(abs(phases[j+1]-phases[j]) for j in range(len(phases)-1))>1.99
    assert abs(weighted[-1])<abs(weighted[0])/10**6


def test_embedded_complement_correction_is_required():
    d=.5;h=lambda x:np.log((1-x+d)/(x+d))
    W=np.array([[1.],[0.]])
    padded=W@np.array([[.2]])@W.T
    raw=functional(padded,h);compressed=W@np.array([[h(.2)]])@W.T
    complement=np.eye(2)-W@W.T
    assert np.allclose(raw,compressed+h(0)*complement)
    assert np.linalg.norm(raw-compressed)>1


def test_regularizer_removal_is_uniform_only_on_the_weighted_core():
    x=np.geomspace(1e-14,.5,1000);x=np.unique(np.r_[x,1-x])
    k=np.log((1-x)/x)
    for d in (.1,.01,.0001):
        kd=np.log((1-x+d)/(x+d));w=x*(1-x)
        assert np.max(w*abs(kd-k))<=d+1e-15
        assert np.max(abs(kd-k))>1
        for t in (.2,2.):
            assert np.max(w*abs(np.exp(-1j*t*kd)-np.exp(-1j*t*k)))<=t*d+1e-15


def test_matching_core_vectors_are_needed_for_generator_limit():
    # Endpoint-free fixed two-mode limit, with a rotating eigenbasis.
    C=np.diag([.1,.7]);g=np.array([1.,2.]);f=C@(np.eye(2)-C)@g
    target=functional(C,lambda x:np.log((1-x)/x))@f
    errors=[]
    for e in (.1,.01,.001):
        V=np.array([[np.cos(e),-np.sin(e)],[np.sin(e),np.cos(e)]])
        A=V@np.diag([.1+e/10,.7-e/10])@V.T
        fa=A@(np.eye(2)-A)@g
        errors.append(np.linalg.norm(functional(A,lambda x:np.log((1-x)/x))@fa-target))
    assert errors[-1]<errors[0]/50


def test_exact_resource_budget_and_simplex_norm():
    assert bounds.core_flow_budget(144,1,F(1,4*144*3**144))==F(9,16)
    # Two maps can be individually isometric without their combined map being so.
    T=np.array([[1.,1.]])
    assert not np.allclose(T.T@T,np.eye(2))
    assert np.linalg.norm(T,2)==pytest.approx(np.sqrt(2))


@pytest.mark.parametrize('call',[
    lambda:bounds.regular_log_budget(0,2,0),
    lambda:bounds.regular_log_budget(.1,2,0),
    lambda:bounds.regular_log_budget(F(1,2),True,0),
    lambda:bounds.core_flow_budget(0,1,0),
    lambda:bounds.core_flow_budget(2,-1,0),
    lambda:bounds.core_flow_budget(2,1,True),
    lambda:bounds.core_regularization_budget(-1,1),
])
def test_invalid_bound_domain(call):
    with pytest.raises(ValueError):call()


@pytest.fixture(scope='module')
def receipt():return packet.build()


@pytest.mark.parametrize('key',['uniform_covariance_gap_assumed','finite_valley_maps_are_exact_isometries',
    'bare_logs_converge_on_all_unmodified_probes','operator_norm_modular_flow_convergence_claimed',
    'absolute_continuum_entropy_computed','shrinking_region_area_or_geometry_limit_proved',
    'nonzero_U_modular_limit_proved','native_modular_gate_execution_supplied',
    'Einstein_or_particle_identification_endpoint_closed'])
def test_unsupported_scope_cannot_be_promoted(receipt,key):
    value=deepcopy(receipt);value['theorem'][key]=True
    with pytest.raises(ArithmeticError):packet.validate(value)


def test_nested_boolean_aliases_constants_and_historical_custody_are_strict(receipt):
    cases=[]
    value=deepcopy(receipt);value['theorem']['source_U']=False;cases.append(value)
    value=deepcopy(receipt);value['algebra']['uniform_spectral_gap_assumed']=0;cases.append(value)
    value=deepcopy(receipt);value['bounds']['uniform_eight_valley_map_norm_upper']='1';cases.append(value)
    value=deepcopy(receipt);value['historical_context']['historical_full_source_replay_performed']=True;cases.append(value)
    value=deepcopy(receipt);value['dependency_custody_pass']=True;cases.append(value)
    for value in cases:
        with pytest.raises(ArithmeticError):packet.validate(value)


def test_historical_snapshot_is_immutable_and_independent_of_git_history(tmp_path,monkeypatch):
    expected=packet.historical_context()
    monkeypatch.chdir(tmp_path)
    monkeypatch.setenv('PATH','')
    assert packet.historical_context()==expected
    forged=tmp_path/'forged.json'
    forged.write_bytes(packet.HISTORICAL.read_bytes()+b' ')
    monkeypatch.setattr(packet,'HISTORICAL',forged)
    with pytest.raises(ArithmeticError,match='immutable historical'):
        packet.historical_context()
