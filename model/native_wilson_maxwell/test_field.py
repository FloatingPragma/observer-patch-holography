"""Independent exact cochain, finite photon and full Gauss controls."""
import itertools,math
import numpy as np
import sympy as s
from scipy.linalg import expm
from . import field,gauss

def test_actual_cube_and_open_box_chain_complex():
    data=field.cube();vertices=data['vertices'];D=data['D'];C=data['C'];G=data['G']
    # Independent scalar graph construction by Hamming adjacency.
    L=s.zeros(8)
    for i,x in enumerate(vertices):
        for j,y in enumerate(vertices):
            if sum(a!=b for a,b in zip(x,y,strict=True))==1:L[i,j]=-1;L[i,i]+=1
    assert D.T*D==L and L*G==s.eye(8)-s.ones(8)/8
    assert G==s.Rational(11,12)*s.eye(8)-s.Rational(1,4)*L+s.Rational(1,48)*L**2-s.Rational(11,96)*s.ones(8)
    assert C.rows==6 and C.rank()==5
    for n in (1,2):
        v,e,f,d,c=field.cochains(n)
        assert d.rank()==len(v)-1 and c.rank()==len(e)-len(v)+1
        assert c*d==s.zeros(c.rows,d.cols)
        d2=field.volume_boundary(n);assert d2*c==s.zeros(n**3,len(e))
    assert list(data['D2'])==[-1,1,-1,1,-1,1]

def test_full_finite_photon_boundary_and_actual_link_norm():
    data=field.cube();cut=field.cutoff(3)
    assert len(cut['states'])==56 and len(cut['annihilators'])==5
    assert s.simplify(cut['Q_edge_norm_squared_at_cap3']-(13+5*s.sqrt(6))/48)==0
    assert s.simplify(s.Rational(9,16)-cut['Q_edge_norm_squared_at_cap3'])>0
    counter=field.boundary_counterexample()
    assert counter['exact_Q_commutator_matrix_element']==s.Rational(3,32)
    Q=field.link_coordinates()[0];vac=s.eye(56)[:,0];vector=vac;c=data['Q_edge_vacuum_variance']
    for degree in range(1,9):
        vector=(Q*vector).applyfunc(field.simple)
        if degree in (2,4,6):
            exact=s.Rational(math.prod(range(1,degree,2)))*c**(degree//2)
            assert field.simple(vector[0]-exact)==0
        if degree==8:
            assert field.simple(vector[0]-81*c**4)==0
            assert field.simple(vector[0]-105*c**4)!=0

def test_literal_cutoff_exponential_and_its_remainder():
    Q=np.array(field.link_coordinates()[0]).astype(complex)
    for coupling in (-.5,.01,.2,1.):
        U=expm(1j*coupling*Q)
        np.testing.assert_allclose(U.conj().T@U,np.eye(56),atol=2e-14)
        assert np.linalg.norm(U-np.eye(56),2)<=abs(coupling)*.75+1e-14
        assert np.linalg.norm(U-np.eye(56)-1j*coupling*Q,2)<=(coupling*.75)**2/2+1e-14
    P1=np.diag([int(sum(n)==1) for n in field.cutoff(3)['states']]);vac=np.eye(56)[:,0]
    # The one-photon projection removes all even Taylor powers from vacuum.
    coupling=.1
    error=P1@(expm(1j*coupling*Q)-1j*coupling*Q)@vac
    assert np.linalg.norm(error)<=abs(coupling*.75)**3/6+1e-14

def test_full_kinematic_partial_shifts_and_Gauss_instrument():
    value=gauss.calculate();assert value['complete_two_endpoint_boundary_cases']==157464
    states=list(itertools.product(range(4),range(-4,5),range(-32,33)))
    index={state:i for i,state in enumerate(states)}
    def constraints(state):
        bits,ell,R=state;return -ell-(bits&1)-R,ell-((bits>>1)&1)
    hop_entries=[];pump_entries=[]
    for col,(bits,ell,R) in enumerate(states):
        if bits==1 and ell<4:hop_entries.append((index[(2,ell+1,R)],col))
        if bits&1==0 and R>-32:pump_entries.append((index[(bits|1,ell,R-1)],col))
    for row,col in hop_entries+pump_entries:assert constraints(states[row])==constraints(states[col])
    assert len({col for row,col in pump_entries})==len(pump_entries)
    # K^dag K and the retained failure effect sum to I on all2340 inputs,
    # including nonphysical and invalid boundary states.
    success=np.zeros(len(states),dtype=int)
    for row,col in pump_entries:success[col]+=1
    assert set(success)=={0,1} and np.all(success+(1-success)==1)
    physical=[state for state in states if constraints(state)==(0,0)]
    assert len(physical)==4
    assert all(ell==((bits>>1)&1) and R==-bits.bit_count() for bits,ell,R in physical)

def test_exact_finite_current_flow_sign():
    K=s.Matrix([[0,0],[1,0]]);H=K+K.H;rho_head=s.diag(0,1)
    J_phase=s.I*(K-K.H);I_flow=-J_phase
    assert s.I*(H*rho_head-rho_head*H)==I_flow
    assert J_phase==-I_flow
    theta=s.Symbol('theta',real=True)
    family=s.exp(s.I*theta)*K+s.exp(-s.I*theta)*K.H
    assert family.diff(theta).subs(theta,0)==J_phase
