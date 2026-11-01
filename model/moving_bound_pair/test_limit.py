"""Independent source-symbol and comparison controls for moving bound pairs."""
from fractions import Fraction as F
from itertools import product
import math
import numpy as np
import pytest
from ..band_edge_contact import limit as edge
from ..bound_pair.fiber import SINGLET, SWAP
from ..contact_scattering.scattering import free_blocks
from . import limit as m


def test_exact_green_majorant():
    c=m.green_comparison_certificate()
    assert c['green_difference_coefficient']==F(1153,540)<F(9,4)
    assert c['rounded_green_difference_coefficient']==3
    assert c['all_edge_valleys']==8


@pytest.mark.parametrize('K', [[0.,0.,0.],[.2,-.3,.4],[.5,.5,.5]])
def test_exact_edge_all_eight_extrema_and_global_upper_bound(K):
    K=np.array(K); speed=.7; B=m.free_edge(K,speed)
    for signs in product([-1,1],repeat=3):
        q=np.array(signs)*math.pi/2
        k=q+K/2
        block=free_blocks(k[None,:],K,speed)[0]
        assert np.linalg.eigvalsh(block)[-1]==pytest.approx(B,abs=2e-15)
    rng=np.random.default_rng(991)
    q=rng.uniform(-math.pi,math.pi,(200,3))
    eigen=np.linalg.eigvalsh(free_blocks(q+K/2,K,speed))
    assert np.max(eigen)<=B+1e-14


def test_pauli_green_formula_independent_four_spin_resolvent():
    rng=np.random.default_rng(881)
    for _ in range(12):
        K=rng.uniform(-.5,.5,3); q=rng.uniform(-math.pi,math.pi,(10,3))
        speed=.9; E=m.free_edge(K,speed)+.3
        actual=np.linalg.solve(E*np.eye(4)-free_blocks(q+K/2,K,speed),
                              np.broadcast_to(SINGLET,(len(q),4))[...,None])[...,0]
        scalar=np.einsum('i,ki->k',SINGLET.conj(),actual)
        assert np.max(abs(scalar-m.singlet_kernel(q,K,E,speed)))<3e-13


def test_pyramid_difference_against_complete_cartesian_torus():
    a=.2; P=np.array([1.,.5,-.25]); K=a*P; z=2.; speed=math.sqrt(3)/a**2
    side=64
    q1=-math.pi/2+(np.arange(side)+.5)*math.pi/side
    q=np.stack(np.meshgrid(q1,q1,q1,indexing='ij'),axis=-1).reshape(-1,3)
    moving=m.singlet_kernel(q,K,m.free_edge(K,speed)+z,speed).mean()
    fixed=m.singlet_kernel(q,np.zeros(3),6/a**2+z,speed).mean()
    direct=(moving-fixed)/a**3
    pyramid=m.green_correction_diagnostic(a,P,z,order=48)
    assert abs(direct-pyramid)<2e-8


@pytest.mark.parametrize('a,P,z', [(.1,[1,.5,0],1.),(.05,[0,1,-1],.4),
                                 (.02,[2,-1,.5],0.),(.001,[1,1,1],2.)])
def test_uniform_exterior_difference_bound(a,P,z):
    correction=m.green_correction_diagnostic(a,P,z,order=24)
    assert abs(correction)<=3*a*max(abs(x) for x in P)**2
    if z>0:
        full=m.scaled_green_diagnostic(a,P,z,order=24)
        assert abs(full+math.sqrt(z)/math.pi)<=a*(11*z+3*max(abs(x) for x in P)**2)


def test_zero_momentum_comparison_exactly_reduces():
    assert m.green_correction_diagnostic(.001,[0,0,0],1.,order=12)==0
    assert m.scaled_green_diagnostic(.001,[0,0,0],1.,order=12)==edge.scaled_green_diagnostic(.001,1.).real


@pytest.mark.parametrize('gamma,cap,time,eps',[(F(1),F(1),F(1),F(1,100)),
 (F(1,10),F(10),F(3),F(1,1000)),(F(3),F(0),F(0),F(1,10)),
 (F(3,2),F(2,7),F(2,3),F(2))])
def test_exact_wavepacket_spacing_schedule(gamma,cap,time,eps):
    cert=m.spacing_certificate(gamma,cap,time,eps)
    assert cert['operator_error_upper']<=eps
    a=cert['spacing_upper']
    assert 256*a*(gamma**2+cap**2)<=gamma
    assert a*gamma<=F(1,256) and a*cap<=F(1,512)


@pytest.mark.parametrize('a,g,p',[(F(1,511),1,1),(F(0),1,1),(.001,1,1),
                                 (F(1,1024),0,1),(F(1,1024),1,-1),(True,1,1)])
def test_unproved_pole_domains_rejected(a,g,p):
    with pytest.raises(ValueError):m.moving_pole_certificate(a,g,p)


@pytest.mark.parametrize('a',[1/512,1/1024])
def test_moving_pole_and_propagation_diagnostics(a):
    P=[1.,.5,-.25]
    coarse=m.pole_diagnostic(a,1.,P,order=20)
    fine=m.pole_diagnostic(a,1.,P,order=40)
    assert abs(coarse['reflected_band_energy']-fine['reflected_band_energy'])<2e-6
    cert=m.propagation_certificate(F(a),1,1,F(1,2))
    assert fine['absolute_energy_error']<float(cert['reflected_band_energy_error_upper'])
    source_phase=np.exp(-.5j*fine['reflected_band_energy'])
    limit_phase=np.exp(-.5j*fine['target'])
    assert abs(source_phase-limit_phase)<float(cert['operator_error_upper'])
    assert .25<=fine['offset']<=4


def test_momentum_dependent_subtraction_would_delete_kinetics():
    a=1/1024; P=[1.,.5,-.25]
    result=m.pole_diagnostic(a,1.,P,order=24)
    wrong_energy=-result['offset']
    assert abs(wrong_energy-result['target'])>.3
    assert result['absolute_energy_error']<.001
    assert result['reflected_band_energy']-wrong_energy==pytest.approx(m.edge_deficit(a,P))


def test_reversing_time_is_not_an_ignorable_scalar():
    result=m.pole_diagnostic(1/4096,1.,[0,0,0],order=8)
    wrong=np.exp(1j*result['reflected_band_energy'])
    correct=np.exp(-1j*result['target'])
    assert abs(wrong-correct)>1.6
    assert abs(np.exp(-1j*result['reflected_band_energy'])-correct)<.001


def test_fixed_lattice_momentum_has_no_uniform_finite_kinetic_limit():
    K=np.array([.25,0,0])
    first=m.edge_deficit(.01,K/.01)
    second=m.edge_deficit(.005,K/.005)
    assert second==pytest.approx(4*first)
    with pytest.raises(ValueError):m.moving_pole_certificate(F(1,100),1,25)


def test_source_bound_fiber_vector_and_exchange():
    side=4; K=np.array([math.pi/2,0,0]); speed=.125; E=2.
    addresses=np.array(list(product(range(side),repeat=3)))
    k=addresses*2*math.pi/side
    blocks=free_blocks(k,K,speed)
    d=np.tile(SINGLET,len(k))/math.sqrt(len(k))
    resolved=np.linalg.solve(E*np.eye(4)-blocks,
             np.broadcast_to(SINGLET/math.sqrt(len(k)),(len(k),4))[...,None])[...,0].reshape(-1)
    U=1/np.vdot(d,resolved).real
    psi=resolved/np.linalg.norm(resolved)
    applied=np.einsum('kij,kj->ki',blocks,psi.reshape(-1,4)).reshape(-1)+U*d*np.vdot(d,psi)
    assert np.linalg.norm(applied-E*psi)<1e-14
    partner=(np.array([1,0,0])-addresses)%side
    index=(partner[:,0]*side+partner[:,1])*side+partner[:,2]
    assert np.linalg.norm(psi.reshape(-1,4)+psi.reshape(-1,4)[index]@SWAP.T)<1e-14
    assert np.vdot(d,psi).real>0


def test_exterior_transition_bound_has_explicit_pole_exclusion():
    cert=m.exterior_transition_certificate(F(1,4096),1,1,4,F(1,2))
    a=1/4096; P=[1.,.5,0]; z=4.
    actual=1/(-1/math.pi-m.scaled_green_diagnostic(a,P,z,order=24))
    target=math.pi/(math.sqrt(z)-1)
    assert abs(actual-target)<float(cert['renormalized_transition_error_upper'])
    with pytest.raises(ValueError):m.exterior_transition_certificate(F(1,100),1,1,4,F(1,2))
    with pytest.raises(ValueError):m.exterior_transition_certificate(F(1,4096),1,1,4,0)


@pytest.mark.parametrize('bad',[[1,2],[True,False,True],[1,np.nan,2],['a','b','c']])
def test_bad_momentum_rejected(bad):
    with pytest.raises(ValueError):m.edge_deficit(.1,bad)
