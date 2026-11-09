"""Literal Clifford, full-Fock and complete-reference symmetry controls."""
from copy import deepcopy

import numpy as np
import pytest
import sympy as s

from . import model as m


def car(count,index):
    a=np.zeros((1<<count,1<<count))
    for state in range(1<<count):
        if state&(1<<index):
            a[state^(1<<index),state]=(-1)**((state&((1<<index)-1)).bit_count())
    return a


def conjugator(spin):
    """Signed-permutation CAR transform, with vacuum mapped to full."""
    count=len(spin);out=np.zeros((1<<count,1<<count))
    assert np.all(np.sum(abs(spin),axis=1)==1)
    for state in range(1<<count):
        occupied=[i for i in range(count) if state&(1<<i)]
        mask=(1<<count)-1;amplitude=1
        for i in reversed(occupied):
            j=int(np.flatnonzero(spin[i])[0])
            amplitude*=spin[i,j]*(-1)**((mask&((1<<j)-1)).bit_count())
            mask^=1<<j
        out[mask,state]=amplitude
    return out


def fock_lift(onebody):
    n=len(onebody);cs=[car(n,i) for i in range(n)]
    result=np.zeros((1<<n,1<<n),complex)
    for i,j in zip(*np.nonzero(onebody)):
        result+=onebody[i,j]*(cs[i].T@cs[j])
    return result


def test_actual_spin_native_and_binary_cube_controls():
    data=m.exact_controls()
    assert data['oriented_Wilson_axes_checked']==3
    assert data['binary_cube_lifts_checked']==48
    assert data['native_integer_charge_trace']==0
    assert data['all_native_blocks_required_for_Gauss_phase_cancellation']


def test_complete_CAR_unitary_includes_both_parities():
    M=np.array(m.density.gamma_data()[1],dtype=float);C=conjugator(M)
    cs=[car(4,j) for j in range(4)]
    assert np.array_equal(C.T@C,np.eye(16)) and np.array_equal(C@C,np.eye(16))
    for i in range(4):
        wanted=sum(M[i,j]*cs[j].T for j in range(4))
        assert np.array_equal(C@cs[i]@C.T,wanted)
    N=sum(c.T@c for c in cs)
    assert np.array_equal(C@N@C.T,4*np.eye(16)-N)
    parity=np.diag([(-1)**i.bit_count() for i in range(16)])
    assert np.array_equal(C@parity,parity@C)


def test_actual_two_site_Wilson_Fock_operator_needs_both_transformations():
    gamma=[np.array(g,dtype=complex) for g in m.density.gamma_data()]
    beta=gamma[3];spin=np.kron(np.eye(2),gamma[1].real)
    C=conjugator(spin);phase=(3+4j)/5
    S=(-beta+beta@gamma[0])/2
    def H(z):
        h=np.kron(np.eye(2),3*beta)
        h[:4,4:]=z*S;h[4:,:4]=z.conjugate()*S.conj().T
        return fock_lift(h)
    h=H(phase);hc=H(phase.conjugate())
    assert np.linalg.norm(C@hc@C.T-h)<1e-12
    assert np.linalg.norm(C@h@C.T-h)>1
    assert np.linalg.norm(hc-h)>1


def test_missing_native_blocks_leave_a_gauge_determinant_phase():
    M=np.array(m.density.gamma_data()[1],dtype=float);C=conjugator(M)
    z=(3+4j)/5
    lifted=np.diag([z**i.bit_count() for i in range(16)])
    inverse=lifted.conj()
    assert np.linalg.norm(C@lifted@C.T-z**4*inverse)<1e-13
    assert np.linalg.norm(C@lifted@C.T-inverse)>1
    blocks=m.channels.native_number_control()['blocks']
    assert sum(b['Q6']*b['internal_dimension'] for b in blocks)==0
    assert all(sum(x['Q6']*x['internal_dimension'] for x in blocks if x is not b)!=0 for b in blocks)


def test_inversion_commutes_only_on_even_CAR_sectors():
    beta=m.density.gamma_data()[3]
    C=conjugator(np.array(m.density.gamma_data()[1],dtype=float))
    lifted=np.zeros((16,16))
    for i in range(16):
        left=[j for j in range(4) if i&(1<<j)]
        for k in range(16):
            right=[j for j in range(4) if k&(1<<j)]
            if len(left)==len(right):lifted[i,k]=beta.extract(left,right).det()
    parity=np.diag([(-1)**i.bit_count() for i in range(16)])
    assert np.array_equal(C@lifted@C.T,parity@lifted)
    assert np.linalg.norm(C@lifted-lifted@C)>1


def test_complete_source_reference_records_preserve_probabilities():
    gamma=[np.array(g,dtype=complex) for g in m.density.gamma_data()]
    C=conjugator(gamma[1].real);H=fock_lift(3*gamma[3]);Q=20*np.eye(16)-H
    root=Q/np.linalg.norm(Q);rho=root@root.conj().T
    assert np.linalg.norm(C@rho@C.T-rho)<1e-14
    b=car(4,0)/2
    failure=np.diag(np.sqrt(1-np.diag(b.T@b)))
    assert np.linalg.norm(b.T@b+failure.T@failure-np.eye(16))<1e-14
    for K in (b,failure):
        conjugate=C@K@C.T
        out=K@root;paired=conjugate@root
        assert np.linalg.norm(paired-C@out@C.T)<1e-14
        assert abs(np.linalg.norm(paired)**2-np.linalg.norm(out)**2)<1e-14
    assert np.linalg.norm(C@root-root)>1  # Source-only action does not fix the purification.


@pytest.fixture(scope='module')
def bound():
    member=m.contract.member(0)
    return member,m.bind_member(member)


def test_actual_members_and_opposite_occupations(bound):
    member,value=bound;r=value['result'];m.verify_bound(member,value)
    assert r['opposite_occupation_sector_spectra_equal']
    assert r['paired_complete_energy_record_probabilities_equal']
    assert r['Hamiltonian_time_not_reversed']
    assert r['spatial_inversion_relation']=='C I C^-1 = CAR_parity I'
    assert not r['conjugation_is_a_physical_antiparticle_identification']
    assert r['physical_mass_table_GeV']==[]
    for n in (1,2,8):
        mem=m.contract.member(n);v=m.bind_member(mem);m.verify_bound(mem,v)
        assert all(x['full_occupation']==2*x['anchor_occupation'] for x in v['result']['native_occupation_conjugation'])


def test_rehashed_forgery_and_wrong_source_are_rejected(bound):
    member,value=bound;bad=deepcopy(value)
    bad['result']['Lorentz_spin_or_CPT_proved']=True
    bad['result_sha256']=m.contract.digest(bad['result'])
    with pytest.raises(ValueError):m.verify_bound(member,bad)
    wrong=deepcopy(member);wrong['source']['spacing']='2'
    with pytest.raises(ValueError):m.bind_member(wrong)
