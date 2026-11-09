"""Exact tent geometry, full-CAR quantum-link controls and source bindings."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
from scipy import sparse
import sympy as s

from . import model as m
from model.particle_local_wilson_density import density


@pytest.mark.parametrize('radius',[1,2,4,8])
def test_literal_periodic_tent_sums_match_closed_form(radius):
    R=radius;L=2*R*R
    values=[max(F(),1-F(min(k,L-k),R)) for k in range(L)]
    q=m.tent_sums(R)
    assert sum(values)==q['one_dimensional_sum']
    assert sum(x*x for x in values)==q['one_dimensional_square_sum']
    assert sum((values[(k+1)%L]-values[k])**2 for k in range(L))==q[
        'one_dimensional_gradient_square_sum']
    assert sum(abs(values[(k+1)%L]-values[k]) for k in range(L))==2
    assert sum(x>0 for x in values)==2*R-1


@pytest.mark.parametrize('radius',[1,2,4])
def test_address_sum_is_exactly_uniform_without_periodic_source_bonds(radius):
    R=radius;L=2*R*R;a=1/R
    t=np.array([max(0,1-min(k,L-k)/R) for k in range(L)])
    profile=t[:,None,None]*t[None,:,None]*t[None,None,:]
    f0=profile;f1=np.roll(profile,1,axis=0)
    k=np.sum(1-np.cos(np.pi*a*(f1-f0)))
    D=float(m.tent_sums(R)['one_axis_three_dimensional_gradient_square_sum'])
    assert 2*a*a*D-1e-12<=k<=5*a*a*D+1e-12
    for shift in ((0,0,0),(L-2,L-1,1)):
        left=np.roll(profile,shift,axis=(0,1,2))
        right=np.roll(left,1,axis=0)
        assert abs(np.sum(1-np.cos(np.pi*a*(right-left)))-k)<1e-12
    if R==1:assert k==pytest.approx(4)
    # The actual source has open interior bonds; none of these operations
    # inserts an edge connecting the two physical faces.
    assert 3*L*L*(L-1)<3*L**3


def _bilinear(modes,i,j):
    rows=[];cols=[];values=[]
    for mask in range(1<<modes):
        if not mask&(1<<j):continue
        sign=(-1)**((mask&((1<<j)-1)).bit_count());state=mask^(1<<j)
        if state&(1<<i):continue
        sign*=(-1)**((state&((1<<i)-1)).bit_count())
        rows.append(state|(1<<i));cols.append(mask);values.append(sign)
    return sparse.csr_matrix((values,(rows,cols)),shape=(1<<modes,1<<modes),dtype=complex)


def _lift(matrix):
    result=sparse.csr_matrix((1<<len(matrix),1<<len(matrix)),dtype=complex)
    for i,j in zip(*np.nonzero(matrix)):result+=matrix[i,j]*_bilinear(len(matrix),int(i),int(j))
    return result


def test_full_car_quantum_link_twists_work_reference_and_cutoff():
    g=density.gamma_data();beta=np.array(g[3],complex)
    S=np.array((-g[3]+g[3]*g[0])/2,complex)
    forward=np.zeros((8,8),complex);forward[:4,4:]=S
    onsite=np.zeros((8,8),complex);onsite[:4,:4]=3*beta;onsite[4:,4:]=3*beta
    hop=_lift(forward);h0=_lift(onsite)
    shift=sparse.csr_matrix([[0,0,0],[1,0,0],[0,1,0]],dtype=complex)
    J=sparse.kron(shift,hop)+sparse.kron(shift.conj().T,hop.conj().T)
    electric=sparse.diags([15.75,0,15.75],dtype=complex)
    H=(J+sparse.kron(sparse.eye(3),h0)+sparse.kron(electric,sparse.eye(256))).toarray()
    N0=np.array([(x&15).bit_count()-2 for x in range(256)])
    N1=np.array([(x>>4).bit_count()-2 for x in range(256)])
    phase=np.tile(np.pi*(N0/3-N1/7),3)
    Q=64*np.eye(768)-H;rho=Q@Q;rho/=np.trace(rho)
    energy=np.trace(rho@H).real;outputs=[];work=0
    rng=np.random.default_rng(710);psi=rng.normal(size=(768,2))+1j*rng.normal(size=(768,2))
    psi/=np.linalg.norm(psi)
    for sign in (-1,1):
        u=np.exp(sign*1j*phase)
        shifted=u.conjugate()[:,None]*H*u[None,:]
        work+=np.trace(rho@shifted).real/2-energy/2
        outputs.append(u[:,None]*psi/np.sqrt(2))
    expected=-(1-np.cos(np.pi*(F(1,3)+F(1,7))))*np.trace(rho@J.toarray()).real
    assert abs(work-expected)<1e-11 and work>0
    reference=sum((x.conj().T@x for x in outputs),np.zeros((2,2),complex))
    assert np.max(abs(reference-psi.conj().T@psi))<1e-13
    # U commutes with an electric-coordinate projection, so performing
    # the source compression does not replace the twist by its mean.
    link_projection=np.repeat([0,1,1],256)
    assert np.max(abs(link_projection*np.exp(1j*phase)-np.exp(1j*phase)*link_projection))==0


def test_corner_star_gram_and_full_fock_singular_value_sum():
    g=density.gamma_data();beta=g[3]
    legs=[(-beta+beta*g[j])/2 for j in range(3)]
    gram=sum((h.H*h for h in legs),s.zeros(4))
    assert gram.eigenvals()=={(3-s.sqrt(3))/2:2,(3+s.sqrt(3))/2:2}
    # A bipartite one-body star has eigenvalues +/-singular values.
    # The full-Fock maximum is their positive sum, over all15 rows.
    sigma_sum=30*(s.sqrt((3-s.sqrt(3))/2)+s.sqrt((3+s.sqrt(3))/2))
    squared=s.simplify(4*sigma_sum**2)
    assert s.simplify(squared-3600*(3+s.sqrt(6)))==0
    assert F(49,20)**2>6
    assert 3600*(3+F(49,20))==19620


def test_tightness_and_positive_mean_allow_an_elastic_limit():
    for n in (1,3,7,12):
        a=F(1,2**n)
        U=np.array([[np.sqrt(float(1-a)),-np.sqrt(float(a))],
                    [np.sqrt(float(a)),np.sqrt(float(1-a))]])
        assert np.max(abs(U.T@U-np.eye(2)))<1e-14
        weights=[1-a,a];energies=[F(),1/a]
        assert sum(weights)==1 and sum(w*E for w,E in zip(weights,energies))==1
        assert sum(w*E*E for w,E in zip(weights,energies))==1/a
        assert sum(w for w,E in zip(weights,energies) if E>10)<=F(1,10)
        if n>=7:assert sum(w for w,E in zip(weights,energies) if 1<=E<=10)==0


def test_complete_origin_restriction_has_original_reference_distance():
    # Orthogonal retained labels make the complete distance exactly the
    # omitted origin fraction, independently of entangled source states.
    rng=np.random.default_rng(17);psi=rng.normal(size=(2,2))+1j*rng.normal(size=(2,2))
    psi/=np.linalg.norm(psi)
    probabilities=np.array([1/4]*4);retained=np.array([0,1/2,1/2,0])
    norms=[]
    for j in range(4):
        U=np.diag([1,np.exp(1j*j)]);out=(U@psi).ravel()
        delta=(probabilities[j]-retained[j])*np.outer(out,out.conjugate())
        norms.append(np.linalg.svd(delta,compute_uv=False).sum())
    assert abs(sum(norms)/2-.5)<1e-14


def test_exact_finite_probability_window_and_phase_error_domains():
    bound=m.positive_window(F(1),F(4))
    assert bound=={'original_transfer_interval':[F(1,4),F(16)],
                   'original_weight_lower':F(1,16)}
    assert m.phase_implementation_bound(2,'1/100000')==F(81,10000)
    for wrong in (True,.1,-1):
        with pytest.raises(ValueError):m.phase_implementation_bound(2,wrong)
    with pytest.raises(ValueError):m.positive_window(2,1)


@pytest.fixture(scope='module')
def receipts():
    return {n:m.bind_member(m.contract.member(n)) for n in (0,1,2,32)}


def test_actual_source_uniform_work_tightness_and_record_budgets(receipts):
    for n,receipt in receipts.items():
        v=receipt['result']
        assert F(v['uniform_original_signed_first_moment_lower'])>F(3,50)
        assert 0<F(v['original_signed_first_moment_lower'])<=F(
            v['original_signed_first_moment_upper'])<=F(14175,32)
        assert v['original_total_spectral_weight']==1
        assert v['uniform_original_probability_tightness_proved']
        assert v['limit_may_be_purely_elastic']
        assert not v['uniform_fixed_positive_band_noncollapse_proved']
        assert not v['source_H_or_boundary_conditions_changed']
        assert not v['full_primary_histogram_executed']
        w=F(v['finite_original_window']['original_weight_lower'])
        assert F(v['finite_reported_original_weight_lower'])>=3*w/4>0
        assert v['physical_mass_table_GeV']==[]
    bulk=receipts[32]['result']['retained_origin']
    assert bulk['profile_unwrapped_on_retained_origins']
    assert F(bulk['signed_first_moment_lower'])>F(1,20)
    assert not bulk['regular_momentum_dictionary_transferred_to_this_probe']


def test_actual_source_and_normalization_mutations_are_rejected(receipts):
    member=m.contract.member(0);m.verify_bound(member,receipts[0])
    wrong=deepcopy(receipts[0]);wrong['result']['original_total_spectral_weight']='1/2'
    with pytest.raises(ValueError):m.verify_bound(member,wrong)
    wrong=deepcopy(receipts[0]);wrong['result']['limit_may_be_purely_elastic']=False
    with pytest.raises(ValueError):m.verify_bound(member,wrong)
    changed=deepcopy(member);changed['source']['side']=4
    with pytest.raises(ValueError):m.bind_member(changed)
    with pytest.raises(ValueError):m.verify_bound(m.contract.member(1),receipts[0])


def test_source_paths_and_proof_control_bytes():
    assert all((m.contract.ROOT/p).is_file() for p in m.SOURCE_INPUTS)
    text=(m.contract.ROOT/'model/particle_primary_local_fluctuation/PROOF.md').read_bytes()
    assert all(c>=32 or c==10 for c in text)
