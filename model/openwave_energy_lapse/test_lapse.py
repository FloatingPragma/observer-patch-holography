from fractions import Fraction as F
import numpy as np
import pytest
from scipy.sparse import diags,kron,eye
from scipy.sparse.linalg import spsolve
from model.native_equilibration.vacuum_repair import hamiltonian
from .source import native,prepare,KAPPA,BITS,energy_interval,source_matrices

def lap(n):
    a=diags([-np.ones(n-1),2*np.ones(n),-np.ones(n-1)],[-1,0,1]);i=eye(n)
    return(kron(kron(a,i),i)+kron(kron(i,a),i)+kron(kron(i,i),a)).tocsr()

def minimize(density,n=5):
    L=lap(n);q=np.zeros(n**3)
    for _ in range(6):
        r=float(KAPPA)*(L@q)+density*np.exp(q)
        q-=spsolve(float(KAPPA)*L+diags(density*np.exp(q)),r)
    return q,float(float(KAPPA)*(q@(L@q))/2+density@np.expm1(q))

def test_native_whole_matrix_ground_and_operator_charge():
    m,levels,gaps=native();a=np.asarray(m,float)
    np.testing.assert_allclose(a,hamiltonian(2,2,F(1))/np.sqrt(14),atol=1e-15,rtol=0)
    e,u=np.linalg.eigh(a);K=a-e[0]*np.eye(8)
    assert abs(e[0]-float((levels[0].lo+levels[0].hi)/2))<1e-14
    for j in range(8):
        assert abs(u[:,j]@K@u[:,j]-float((gaps[j].lo+gaps[j].hi)/2))<1e-14
    rng=np.random.default_rng(612)
    for _ in range(20):
        v=rng.normal(size=8)+1j*rng.normal(size=8);v/=np.linalg.norm(v)
        assert np.vdot(v,K@v).real>=-1e-14
    q=-.0001;d=(np.exp(q+1e-6)*K-np.exp(q-1e-6)*K)/(2e-6)
    np.testing.assert_allclose(d,np.exp(q)*K,atol=5e-10,rtol=0)
    _,_,_,oldJ=source_matrices()
    assert np.linalg.norm(a@np.asarray(oldJ,float)-np.asarray(oldJ,float)@a)>.1

def test_nonlinear_finite_minimum_contains_independent_newton_solution():
    r=prepare(5,(((2,2,2),1),((3,2,2),7)));density=np.zeros(125)
    for site,k in(((2,2,2),1),((3,2,2),7)):
        e=energy_interval(k);density[np.ravel_multi_index(site,(5,5,5))]=float((e.lo+e.hi)/2)
    q,minimum=minimize(density)
    lo,hi=map(F,r['minimum_energy_interval'])
    assert float(lo)-1e-18<=minimum<=float(hi)+1e-18
    assert np.max(abs(q-r['integer_geometry'].ravel()/2**BITS))<1e-11
    assert F(r['fixed_point_contraction_rational_upper'])<F(r['window'])

def test_mixed_minimum_derivative_against_independent_variations():
    n=5;a=np.ravel_multi_index((2,2,2),(n,n,n));b=np.ravel_multi_index((3,2,2),(n,n,n))
    density=np.zeros(n**3);density[a]=.7;density[b]=1.3;q,_=minimize(density)
    h=float(KAPPA)*lap(n)+diags(density*np.exp(q));v=np.zeros(n**3);v[b]=1
    target=-np.exp(q[a]+q[b])*spsolve(h,v)[a]
    step=.001;values={}
    for sa in(-1,1):
        for sb in(-1,1):
            rho=density.copy();rho[a]+=sa*step;rho[b]+=sb*step;values[sa,sb]=minimize(rho)[1]
    finite=(values[1,1]-values[1,-1]-values[-1,1]+values[-1,-1])/(4*step**2)
    assert target<0 and abs(finite-target)<1e-12

def test_stiffness_and_geometry_mass_are_real_controls():
    L=lap(5);rho=np.zeros(125);rho[62]=1
    q=-spsolve(L,rho);half=-spsolve(2*L,rho);screened=-spsolve(L+10*eye(125),rho)
    np.testing.assert_allclose(half,q/2,atol=1e-15,rtol=0)
    assert abs(screened[0])<abs(q[0])/20

def test_invalid_geometry_and_source_inputs_fail():
    for kwargs in({'n':4},{'kappa':0},{'sources':(((2,2,2),0),)},{'sources':(((99,2,2),1),)},{'window':F(1,10)}):
        with pytest.raises(ValueError):prepare(**kwargs)
