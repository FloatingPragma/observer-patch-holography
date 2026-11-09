"""Independent full-Fock, quantum-link and complete-record controls."""
from copy import deepcopy
from fractions import Fraction as F

import numpy as np
import pytest
from scipy import sparse
import sympy as s

from . import model as r


def _car(modes,i,j):
    rows=[];cols=[];values=[]
    for mask in range(1<<modes):
        if not mask&(1<<j):continue
        sign=(-1)**((mask&((1<<j)-1)).bit_count());reduced=mask^(1<<j)
        if reduced&(1<<i):continue
        sign*=(-1)**((reduced&((1<<i)-1)).bit_count())
        rows.append(reduced|(1<<i));cols.append(mask);values.append(sign)
    return sparse.csr_matrix((values,(rows,cols)),shape=(1<<modes,1<<modes),dtype=complex)


@pytest.mark.parametrize('overlap',[s.Rational(0),s.Rational(3,5),s.Rational(1)])
def test_all_occupation_spectrum_including_two_particle_maximum(overlap):
    car=lambda i,j:s.Matrix(_car(4,i,j).real.astype(int).toarray())
    root=s.sqrt(1-overlap**2)
    operator=car(2,1)-overlap*car(0,2)-root*car(0,3)
    square=operator.H*operator
    z=s.Symbol('lambda')
    expected=z**8*(z-1)**4*(z*z-2*z+overlap**2)**2
    assert s.expand(square.charpoly(z).as_expr()-expected)==0
    for number in range(5):
        indices=[j for j in range(16) if j.bit_count()==number]
        eigen=np.linalg.eigvalsh(np.array(square.extract(indices,indices),dtype=float))
        target=0 if number in (0,4) else 1 if number in (1,3) else float(1+root)
        assert abs(eigen[-1]-target)<1e-13
    if overlap==0:
        # A one-body bound misses the true full-Fock maximum.
        assert max(square.eigenvals())==2


def test_actual_source_corner_modes_are_not_orthogonal():
    data=r.source_probe_control()
    assert len(data['actual_corner_Gram_controls'])==16
    assert {row['neighbor_overlap_squared'] for row in data['actual_corner_Gram_controls']}=={'1/4'}
    assert data['native_dimension']==1
    beta=r.density.gamma_data()[3]
    for p,m in r.beta_pairs():
        b=p*m.H
        assert beta*b-b*beta==2*b
        assert b.H*b*b.H*b==b.H*b


def test_actual_wilson_pair_commutator_with_quantum_link_compression():
    modes=8;dim=1<<modes
    beta=np.array(r.density.gamma_data()[3],dtype=complex)
    gamma=np.array(r.density.gamma_data()[0],dtype=complex)
    spin=(-beta+beta@gamma)/2
    onsite=sum((3*beta[i,j]*(_car(modes,i,j)+_car(modes,4+i,4+j))
                for i in range(4) for j in range(4)),sparse.csr_matrix((dim,dim)))
    hop=sum((spin[i,j]*_car(modes,i,4+j) for i in range(4) for j in range(4)),
            sparse.csr_matrix((dim,dim)))
    p,m=(np.array(v,dtype=complex).ravel() for v in r.beta_pairs()[0])
    b=sum((p[i]*m[j].conjugate()*_car(modes,i,j) for i in range(4) for j in range(4)),
          sparse.csr_matrix((dim,dim)))
    assert sparse.linalg.norm(onsite@b-b@onsite-6*b)<1e-12
    link=sparse.csr_matrix([[0,0,0],[1,0,0],[0,1,0]],dtype=complex)
    electric=sparse.diags([15.75,0,15.75],format='csr')
    H=sparse.kron(sparse.eye(3),onsite)+sparse.kron(electric,sparse.eye(dim))
    H+=sparse.kron(link,hop)+sparse.kron(link.H if hasattr(link,'H') else link.conj().T,hop.conj().T)
    B=sparse.kron(sparse.eye(3),b)
    residual=H@B-B@H-6*B
    bound=sparse.linalg.eigsh(residual.conj().T@residual,k=1,which='LA',return_eigenvectors=False)[0]
    assert .1<bound<=1+1e-11
    # The link is a finite quantum operator; it is not replaced by a phase.
    assert sparse.linalg.norm(link.conj().T@link-sparse.eye(3))>0


def _pair_controls():
    return tuple(s.Matrix(_car(4,i,i+2).real.astype(int).toarray()) for i in range(2))


def test_complete_probe_retains_failure_on_arbitrary_reference():
    pairs=_pair_controls();rows=r.averaged_instrument(pairs)
    assert len(rows)==3
    assert sum((k.H*k for k in rows),s.zeros(16))==s.eye(16)
    psi=s.zeros(32,1);psi[0]=s.sqrt(s.Rational(1,3));psi[21]=s.sqrt(s.Rational(2,3))
    outputs=[s.kronecker_product(k,s.eye(2))*psi for k in rows]
    assert sum((v.H*v)[0] for v in outputs)==1
    expected=sum((psi.H*s.kronecker_product(b.H*b,s.eye(2))*psi)[0] for b in pairs)/2
    assert s.simplify(sum((v.H*v)[0] for v in outputs[:-1])-expected)==0


def test_aggregate_failure_is_not_separate_failure_channel():
    pairs=_pair_controls();rows=r.averaged_instrument(pairs)
    psi=s.zeros(16,1);psi[0]=psi[4]=1/s.sqrt(2)
    aggregate=rows[-1]*psi*psi.H*rows[-1].H
    separate=sum(((s.eye(16)-b.H*b)*psi*psi.H*(s.eye(16)-b.H*b)/2 for b in pairs),s.zeros(16))
    assert s.simplify(s.trace(aggregate-separate))==0
    assert r.apparatus.clean(aggregate-separate)!=s.zeros(16)


def test_spectral_second_moment_and_transfer_law_in_stationary_mixed_state():
    energies=np.array([0.,3.,6.,9.]);prob=np.array([.7,.1,.15,.05])
    H=np.diag(energies);rho=np.diag(prob)
    B=np.zeros((4,4),dtype=complex);B[2,0]=1;B[3,1]=.5;B[0,3]=.25
    residual=H@B-B@H-6*B
    spectral=sum(prob[i]*abs(B[j,i])**2*(energies[j]-energies[i]-6)**2 for i in range(4) for j in range(4))
    assert abs(spectral-np.trace(rho@residual.conj().T@residual))<1e-14
    success=sum(prob[i]*abs(B[j,i])**2 for i in range(4) for j in range(4))
    assert success<1
    band=sum(prob[i]*abs(B[j,i])**2 for i in range(4) for j in range(4) if 2<=energies[j]-energies[i]<=10)
    assert band>=success-spectral/16
    # Downward transfers exist in a stationary mixed state.
    assert prob[3]*abs(B[0,3])**2>0


def test_complete_two_clock_probability_law_with_idle_entangled_reference():
    H=s.diag(0,1);rho=s.diag(s.Rational(3,4),s.Rational(1,4))
    U=s.diag(1,-1);energy=r.spectroscopy.energy_instrument(U,2)
    probe=r.averaged_instrument((s.Matrix([[0,0],[1,0]]),))
    psi=s.Matrix([s.sqrt(3)/2,0,0,s.Rational(1,2)])
    total=0
    for a,A in enumerate(energy):
        for b,K in enumerate(probe):
            for c,C in enumerate(energy):
                output=s.kronecker_product(C*K*A,s.eye(2))*psi
                weight=(output.H*output)[0]
                law=sum(rho[i,i]*abs(K[j,i])**2*abs(A[i,i])**2*abs(C[j,j])**2
                        for i in range(2) for j in range(2))
                assert s.simplify(weight-law)==0
                total+=weight
    assert total==1
    reference_after=sum((s.kronecker_product(A,s.eye(2))*psi*psi.H*s.kronecker_product(A.H,s.eye(2))
                         for A in energy),s.zeros(4))
    assert reference_after!=psi*psi.H


@pytest.fixture(scope='module')
def member():return r.contract.member(0)


def test_actual_selected_original_band_and_finite_complete_clock_bound(member):
    receipt=r.bind_member(member);r.verify_bound(member,receipt)
    value=receipt['result']
    assert F(value['original_success_weight_lower'])>F(11,16)
    assert F(value['original_transfer_weight_lower'])>F(1,2)
    assert F(value['original_transfer_weight_simple_lower'])==F(1,2)
    assert value['original_transfer_interval']==[2,10]
    assert value['original_final_excitation']['inferred_final_excitation_interval']==['2','21/2']
    assert F(value['original_final_excitation_weight_simple_lower'])==F(3,8)
    assert value['finite_reported_transfer_interval']==['15/8','81/8']
    assert F(value['finite_reported_original_event_weight_lower'])>F(15,32)
    assert F(value['complete_clock_and_tail_loss_upper'])<=F(1,32)
    J=value['finite_clock_plan']['clock_plan']['clock_size']
    assert value['complete_fine_history_outcomes']==17*J*J
    assert value['all_failures_and_reference_retained']
    assert not value['original_weight_conditioned_on_success']
    assert not value['blank_or_trial_substituted_for_selected_state']
    assert not value['numerical_ground_energy_subtraction_required']
    assert not value['actual_primary_histogram_executed']
    assert value['physical_mass_table_GeV']==[]


def test_bound_and_source_mutations_rejected(member):
    receipt=r.bind_member(member)
    changed=deepcopy(receipt);changed['result']['original_transfer_weight_simple_lower']='1'
    with pytest.raises(ValueError):r.verify_bound(member,changed)
    changed=deepcopy(member);changed['clock']['GeV_per_tick']='1'
    with pytest.raises(ValueError):r.bind_member(changed)
    with pytest.raises(ValueError):r.bind_member(r.contract.member(1))


def test_evidence_paths_and_invalid_instruments():
    assert all((r.contract.ROOT/name).is_file() for name in r.SOURCE_INPUTS)
    with pytest.raises(ValueError):r.averaged_instrument(())
    with pytest.raises(ValueError):r.averaged_instrument((2*s.eye(2),))
