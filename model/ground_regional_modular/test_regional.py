from fractions import Fraction as F

import numpy as np
import pytest

from . import regional as r
from ..interacting_vacuum import ground
from ..spin_resolved_spectrum import spin


def full_fock_oracle(sites):
    kept=[2*x+s for x in sites for s in (0,1)]
    order=kept+[x for x in range(16) if x not in kept]
    rank={x:j for j,x in enumerate(order)}
    out=np.zeros((1<<len(kept),)*2,complex)
    for norm,state in zip(spin.NORMS,spin.ground_family(),strict=True):
        vector=np.zeros(1<<16,complex)
        for mask,(a,b) in state.items():
            occupied=[j for j in range(16) if mask&(1<<j)]
            inversions=sum(rank[x]>rank[y] for i,x in enumerate(occupied) for y in occupied[i+1:])
            vector[mask]=(-1)**inversions*(a+np.sqrt(3)*b)/(36*np.sqrt(norm))
        tensor=vector.reshape((2,)*16,order='F').transpose(order)
        split=tensor.reshape(1<<len(kept),-1,order='F')
        out+=split@split.conj().T/5
    return out


def test_exact_source_census_and_complete_support():
    value=r.certificate()
    assert len(value['single_sites'])==8
    assert len(value['all_site_pairs'])==28
    assert sum(bool(x['bond_sign']) for x in value['all_site_pairs'])==12
    assert all(x['support_rank']+x['kernel_rank']==16 for x in value['all_site_pairs'])


@pytest.mark.parametrize('sites',[(0,),(0,1),(1,6),(0,1,4),(0,2,7)])
def test_signed_trace_against_independent_full_fock_tensor(sites):
    assert np.max(np.abs(r.floating(r.density(sites))-full_fock_oracle(sites)))<3e-15


@pytest.mark.parametrize('sites',[(0,1,2),(0,2,7),(1,3,5)])
def test_general_natural_orbital_spectrum_on_three_sites(sites):
    a=ground.band_data()[0][np.ix_(sites,sites)]
    c=np.eye(3)/2-a/(2*np.sqrt(3))
    predicted=r.generic_spectrum(np.linalg.eigvalsh(c))
    actual=np.linalg.eigvalsh(r.floating(r.density(sites)))
    assert np.max(np.abs(actual-predicted))<5e-15
    assert np.count_nonzero(actual>1e-12)==20


def test_full_cube_and_empty_orbitals_are_not_discarded():
    spectrum=r.generic_spectrum([0]*4+[1]*4)
    assert len(spectrum)==65536
    assert np.count_nonzero(spectrum)==5
    assert np.array_equal(spectrum[-5:],np.full(5,0.2))


def test_spin_integral_normalizes_each_occupied_subset():
    for k in range(9):
        assert r.spin_twirl_weights(k)==(F(1,k+1),)*(k+1)


def supported_data():
    rho=r.floating(r.density((0,1)))
    eig,vec=np.linalg.eigh(rho)
    q=vec[:,eig>1e-12];weights=eig[eig>1e-12]
    h=r.pair_hopping(1)/16
    n=np.diag([x.bit_count() for x in range(16)])
    central=np.diag(np.log1p(np.diag(n)))
    k=np.log(6)*np.eye(16)+16*np.log(2+np.sqrt(3))*h+central
    return rho,q,weights,h,central,k


def test_supported_modular_identity_and_actual_central_domain():
    rho,q,weights,h,central,k=supported_data()
    assert np.max(np.abs(q.T@k@q-np.diag(-np.log(weights))))<4e-14
    assert np.max(np.abs(central@h-h@central))==0
    assert np.linalg.norm(central@r.creator(4,0)-r.creator(4,0)@central)>0.5
    assert weights.min()>1/48
    # A bare energy plus scalar cannot give different vacuum and triplet K.
    assert h[0,0]==h[5,5]==0
    assert abs(k[5,5]-k[0,0]-np.log(3))<1e-15


def entropy_positive(x):
    values=np.linalg.eigvalsh(x)
    assert np.min(values)>0
    return -np.dot(values,np.log(values))


def test_noncommuting_first_law_and_uniform_remainder():
    rho,q,weights,h,central,k=supported_data()
    psi=np.arange(1,9)+1j*np.arange(8,0,-1)
    psi=psi/np.linalg.norm(psi)
    tau=np.diag(weights);sigma=np.outer(psi,psi.conj());delta=sigma-tau
    assert np.linalg.norm(tau@delta-delta@tau)>0.01
    derivative=np.trace(delta@(q.T@k@q)).real
    pairing=16*np.log(2+np.sqrt(3))*np.trace(delta@(q.T@h@q)).real+np.trace(delta@(q.T@central@q)).real
    assert abs(derivative-pairing)<5e-15
    for epsilon in (1/128,1/64,1/32):
        changed=tau+epsilon*delta
        relative=epsilon*derivative-entropy_positive(changed)+entropy_positive(tau)
        bound=24*epsilon**2*np.linalg.norm(delta)**2/(1-epsilon)
        assert 0<relative<bound
    step=1e-6
    numerical=(entropy_positive(tau+step*delta)-entropy_positive(tau-step*delta))/(2*step)
    assert abs(numerical-derivative)<2e-8


def test_outside_support_has_singular_entropy_derivative():
    _,_,weights,_,_,_=supported_data()
    base=-np.dot(weights,np.log(weights))
    slopes=[]
    for epsilon in (1e-3,1e-5,1e-7):
        spectrum=np.r_[(1-epsilon)*weights,epsilon]
        difference=-np.dot(spectrum,np.log(spectrum))-base
        h2=-epsilon*np.log(epsilon)-(1-epsilon)*np.log1p(-epsilon)
        assert abs(difference-(-epsilon*base+h2))<8e-16
        slopes.append(difference/epsilon)
    assert slopes[2]>slopes[1]+4 and slopes[1]>slopes[0]+4


@pytest.mark.parametrize('bad',[(),(1,0),(0,0),(-1,),(8,),(True,),[0]])
def test_invalid_regions_are_rejected(bad):
    with pytest.raises(ValueError):r.density(bad)


def test_invalid_occupation_and_unretained_dense_domains():
    for bad in ([-0.1],[1.1],[float('nan')],[[0.5]]):
        with pytest.raises(ValueError):r.generic_spectrum(bad)
    with pytest.raises(ValueError):r.density((0,1,2,3))
