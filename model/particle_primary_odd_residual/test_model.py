"""Independent CAR/Laurent actions, full occupation forms and source joins."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as s

from . import model as m


def _move(mask,i,j):
    if not mask&(1<<j):return None
    sign=(-1)**((mask&((1<<j)-1)).bit_count());mask^=1<<j
    if mask&(1<<i):return None
    sign*=(-1)**((mask&((1<<i)-1)).bit_count())
    return mask|(1<<i),sign


def _spin():
    pairs=m.response.beta_pairs()
    chart=s.Matrix.hstack(pairs[0][0],pairs[1][0],pairs[0][1],pairs[1][1])
    gamma=m.response.density.gamma_data();beta=gamma[3]
    return [np.array(s.simplify(chart.H*((-beta+beta*g)/2)*chart),dtype=complex) for g in gamma[:3]]


def test_literal_32_mode_CAR_with_tuple_Haar_exponents_matches_second_action():
    spin=_spin();sites=list(product(range(2),repeat=3));where={x:j for j,x in enumerate(sites)}
    zero=(0,)*12;edges=[];terms=[]
    for x in sites:
        for axis in range(3):
            if x[axis]:continue
            y=list(x);y[axis]=1;y=tuple(y);edge=len(edges);edges.append((x,y))
            for i,j in product(range(4),repeat=2):
                c=spin[axis][i,j]
                if c:
                    terms.extend(((4*where[x]+i,4*where[y]+j,edge,1,c),
                                  (4*where[y]+j,4*where[x]+i,edge,-1,c.conjugate())))
    omega=sum(1<<(4*x+j) for x in range(8) for j in(2,3));phi={}
    for x,j in product(range(8),range(2)):
        mask,sign=_move(omega,4*x+j,4*x+2+j);phi[(mask,zero)]=sign/4
    def act(vector):
        out={}
        for (mask,exponent),a in vector.items():
            for i,j,edge,step,c in terms:
                moved=_move(mask,i,j)
                if moved:
                    e=list(exponent);e[edge]+=step;key=(moved[0],tuple(e))
                    out[key]=out.get(key,0)+a*c*moved[1]
        return {key:a for key,a in out.items() if a!=0}
    first=act(phi);parts={2:{},4:{}}
    for key,a in first.items():
        mask=key[0]
        D=sum(bool(mask&(1<<(4*x+j))) for x in range(8) for j in(0,1))
        D+=sum(not mask&(1<<(4*x+j)) for x in range(8) for j in(2,3))
        parts[D][key]=a
    assert [sum(abs(a)**2 for a in parts[k].values()) for k in(2,4)]==[3,10.5]
    second={k:act(v) for k,v in parts.items()}
    actual=[[sum(a.conjugate()*second[j].get(key,0) for key,a in second[k].items())
             for j in(2,4)] for k in(2,4)]
    assert actual==[[46.5,61.5],[61.5,348]]
    assert [len(second[k]) for k in(2,4)]==[1928,16456]
    # The first norms are constant even at nontrivial link coordinates.
    phases=np.exp(1j*np.arange(1,13)*.37)
    for k in(2,4):
        vector={}
        for (mask,e),a in parts[k].items():
            vector[mask]=vector.get(mask,0)+a*np.prod(phases**np.array(e))
        assert abs(sum(abs(a)**2 for a in vector.values())-({2:3,4:10.5}[k]))<1e-12
    control=m.q6_second_action_control()
    assert control['parent_second_action_Gram']==[[F(93,2),F(123,2)],[F(123,2),F(348)]]


def _action(h,vector):
    terms=[(i,j,h[i,j]) for i,j in product(range(h.shape[0]),repeat=2) if h[i,j]]
    out={}
    for mask,a in vector.items():
        for i,j,c in terms:
            moved=_move(mask,i,j)
            if moved:out[moved[0]]=out.get(moved[0],0)+a*c*moved[1]
    return {mask:a for mask,a in out.items() if a!=0}


def _sum_vectors(*vectors):
    out={}
    for vector in vectors:
        for key,a in vector.items():out[key]=out.get(key,0)+a
    return {key:a for key,a in out.items() if a!=0}


def _norm(vector):return sum(abs(a)**2 for a in vector.values())


def test_independent_16_mode_two_native_factorization_including_magnetic_cross_terms():
    h=np.zeros((8,8),complex);h[:4,4:]=1j*_spin()[1];h[4:,:4]=h[:4,4:].conj().T
    X=np.zeros((16,16),complex);X[:8,:8]=h
    Y=np.zeros_like(X);Y[8:,8:]=h
    omega=sum(1<<j for j in(2,3,6,7,10,11,14,15));phi={}
    for site,j in product((8,12),range(2)):
        mask,sign=_move(omega,site+j,site+j+2);phi[mask]=sign/2
    yf=_action(Y,phi);parts={2:{},4:{}}
    for mask,a in yf.items():
        D=sum(bool(mask&(1<<j)) for j in(8,9,12,13))
        D+=sum(not mask&(1<<j) for j in(10,11,14,15))
        parts[D][mask]=a
    tx,t2,t4=.25,.5,.75
    W=_sum_vectors({k:t2*v for k,v in parts[2].items()}, {k:t4*v for k,v in parts[4].items()})
    Xt={k:tx*v for k,v in _action(X,phi).items()}
    zX=1;sX=tx*tx;qX=tx;MX=2.5*tx*tx
    z6=_norm(yf);s6=_norm(W);q6=sum(a.conjugate()*W.get(k,0) for k,a in yf.items()).real
    M6=_norm(_action(Y,W))
    expected=MX+M6+zX*s6+sX*z6+4*qX*q6
    assert abs(_norm(_action(X+Y,_sum_vectors(Xt,W)))-expected)<1e-12
    # A scalar magnetic coefficient commutes with the full matter algebra.
    mag,mag_t=.5,.25
    vt=_sum_vectors(Xt,W,{k:mag_t*v for k,v in phi.items()})
    actual=_norm(_sum_vectors(_action(X+Y,vt),{k:mag*v for k,v in vt.items()}))
    predicted=expected+(zX+z6)*mag_t**2+(sX+s6)*mag**2+4*(qX+q6)*mag*mag_t+mag**2*mag_t**2
    assert abs(actual-predicted)<1e-12


def _creator(modes,j):
    out=np.zeros((1<<modes,1<<modes))
    for mask in range(1<<modes):
        if not mask&(1<<j):out[mask|(1<<j),mask]=(-1)**((mask&((1<<j)-1)).bit_count())
    return out


def test_pair_creation_and_trace_zero_number_forms_on_every_small_Fock_sector():
    creators=[_creator(6,j) for j in range(6)];numbers=[x@x.T for x in creators]
    D=sum(numbers);L=sum(lam*creators[j]@creators[j+3] for j,lam in enumerate((1,2,3)))
    assert np.linalg.eigvalsh(14*(np.eye(64)+D/2)-L.T@L)[0]>-1e-11
    assert (L.T@L)[0,0]==14  # Removing the vacuum term breaks the form bound.
    N=numbers[0]+numbers[1]+numbers[2]-numbers[3]-numbers[4]-numbers[5]
    assert np.linalg.eigvalsh(3*D-N.T@N)[0]>-1e-12
    assert max(np.linalg.eigvalsh(N))==3
    # The hole estimate requires zero trace of the one-body coefficient.
    assert D[-1,-1]==6 and (D.T@D)[-1,-1]>3*D[-1,-1]


def test_actual_full_cube_one_body_commutator_and_coherent_normalization():
    spin=_spin();sites=list(product(range(2),repeat=3));where={x:j for j,x in enumerate(sites)}
    h=np.zeros((32,32),complex);edge=0
    for x in sites:
        for axis in range(3):
            if x[axis]:continue
            y=list(x);y[axis]=1;y=tuple(y);i=4*where[x];j=4*where[y]
            block=np.exp(.23j*(edge+1))*spin[axis];edge+=1
            h[i:i+4,j:j+4]=block;h[j:j+4,i:i+4]=block.conj().T
    p=[4*x+j for x in range(8) for j in(0,1)];o=[4*x+j for x in range(8) for j in(2,3)]
    h=h[np.ix_(p+o,p+o)];S=np.zeros_like(h);S[:16,16:]=np.eye(16)
    r=h@S-S@h;Bstar=h[16:,:16]
    assert np.allclose(r[:16,:16],-Bstar) and np.allclose(r[16:,16:],Bstar)
    assert np.allclose(r[16:,:16],0)
    assert abs(np.linalg.norm(r[:16,16:],'fro')**2-48)<1e-12
    assert abs(np.trace(Bstar))<1e-14 and np.linalg.norm(Bstar,2)<=1.5+1e-12
    assert max((k+1)*(16-k) for k in range(16))==72
    assert F(72,512)==F(9,64)
    assert m.coherent_commutator_control()['probe_normalization_squared']==512


@pytest.fixture(scope='module')
def source_bound():
    member=m.contract.member(0)
    return member,m.bind_member(member)


def test_actual_source_odd_eigenvalue_window_excludes_its_ground(source_bound):
    result=source_bound[1]['result']['odd_trial_residual']
    assert F('5.365')<F(result['source_residual_upper'])<F('5.366')
    assert F(result['inherited_source_residual_upper'])>100
    assert F(result['eigenvalue_interval_lower_minus_ground_upper'])>F('0.5479')
    assert F(result['some_A1u_excitation_interval'][1])<F('15.009')
    assert not result['source_residual_squared_is_exact_after_cutoff']
    assert not result['trial_replaces_selected_polynomial_state']
    assert not result['selected_probe_weight_in_this_eigenvalue_interval_proved']
    assert not result['isolated_eigenvalue_or_particle_mass_proved']


def test_actual_selected_coherent_weight_and_full_failure_budgets(source_bound):
    r=source_bound[1]['result']['selected_coherent_response']
    assert F(r['q6_defect_expectation_upper'])<F(46,25)
    assert F(r['selected_centered_second_moment_upper'])<F(13767,22400)
    assert F(r['original_transfer_weight_lower'])>F(7057,1317120)>F(5,1024)
    assert F(r['original_final_excitation_weight_lower'])>F(10151,10536960)>F(1,2048)
    assert F(r['cooling_original_weight_loss_upper'])==F(9,2048)
    assert r['original_transfer_interval']==['3/4','45/4']
    assert r['original_final_excitation_interval']==['3/4','53/4']
    assert r['finite_reported_transfer_interval']==['23/32','361/32']
    assert F(r['finite_reported_original_weight_lower'])>F(39,8192)
    assert r['all_failure_records_and_reference_retained']
    assert not r['probabilities_conditioned_on_success']
    assert not r['final_state_cubic_irrep_is_fixed_by_operator_irrep']
    assert not r['response_attached_to_the_trial_eigenvalue']
    enclosing=source_bound[1]['result']['odd_trial_residual']['some_A1u_excitation_interval']
    assert F(enclosing[0])<F(3,4)<F(53,4)<F(enclosing[1])


def test_complete_coherent_probe_retains_entangled_reference_and_failure():
    T=m.branch.channels.car_bilinear(4,0,2)+m.branch.channels.car_bilinear(4,1,3)
    K=T/(16*s.sqrt(2));effect=K.H*K
    values=np.linalg.eigvalsh(np.array(effect,dtype=complex))
    assert values.max()<1 and values.min()>-1e-14
    rootvalues,rootvectors=np.linalg.eigh(np.eye(16)-np.array(effect,dtype=complex))
    failure=(rootvectors*np.sqrt(rootvalues))@rootvectors.conj().T
    K=np.array(K,dtype=complex);psi=np.eye(16,dtype=complex).reshape(-1)/4
    success=np.kron(K,np.eye(16))@psi;failed=np.kron(failure,np.eye(16))@psi
    assert abs(np.vdot(success,success)+np.vdot(failed,failed)-1)<1e-12
    assert np.vdot(failed,failed).real>.99
    assert abs(np.vdot(success,success)-np.trace(K.conj().T@K)/16)<1e-12


def test_residual_does_not_force_a_single_line_or_a_selected_probe_overlap():
    H=np.diag([-10.,-4.,-2.]);trial=np.array([0.,1.,1.])/np.sqrt(2)
    center=trial@H@trial;radius=np.linalg.norm((H-center*np.eye(3))@trial)
    assert (center,radius)==pytest.approx((-3,1))
    assert sum(center-radius-1e-12<=e<=center+radius+1e-12 for e in np.diag(H))==2
    # A stationary faithful preparation does not identify a probe residue:
    # this complete probe acts only on the disjoint ground coordinate.
    rho=np.diag([.8,.1,.1]);K=np.diag([.5,0,0]);P=np.diag([0.,1.,1.])
    assert np.trace(rho@K.T@P@K)==0 and trial@P@trial>0.999


def test_source_and_claim_mutations_are_rejected(source_bound):
    member,receipt=source_bound;m.verify_bound(member,receipt)
    forged=deepcopy(receipt)
    forged['result']['selected_coherent_response']['response_attached_to_the_trial_eigenvalue']=True
    with pytest.raises(ValueError):m.verify_bound(member,forged)
    forged=deepcopy(member);forged['state']['target_ground_excess']='0'
    with pytest.raises(ValueError):m.bind_member(forged)
    with pytest.raises(ValueError):m.bind_member(m.contract.member(1))
    with pytest.raises(ValueError):m.odd_parent_moment([1.]*6,[1]*5)
    with pytest.raises(ValueError):m.odd_parent_moment([1]*5,[1]*5)
    assert all((m.contract.ROOT/path).is_file() for path in m.SOURCE_INPUTS)
