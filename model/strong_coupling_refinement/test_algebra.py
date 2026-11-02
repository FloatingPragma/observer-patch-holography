import numpy as np
import pytest
from scipy.linalg import expm, eigh
from . import algebra, bounds


@pytest.mark.parametrize('sides', [(2,1,1),(3,1,1),(2,2,1),(2,2,2)])
def test_complete_half_filled_car_denominator_and_exchange(sides):
    row = algebra.exact_audit(sides)
    assert row['PTP_zero'] and row['one_hop_denominator_exactly_U'] and row['second_order_exchange_exact']


def test_cube_contains_signed_flux_but_exchange_is_sign_independent():
    n, edges = algebra.geometry((2,2,2))
    assert any(e[2] == -1 for e in edges)
    masks,T,D = algebra.operators(n, edges)
    P = algebra.spin_embedding(n,masks)
    positive = tuple((x,y,1) for x,y,s in edges)
    assert (algebra.exchange_four(n,edges)-algebra.exchange_four(n,positive)).nnz == 0
    # Dropping negative source links changes the actual carrier hopping.
    _,Tpositive,_ = algebra.operators(n,positive)
    assert (T-Tpositive).nnz > 0
    # Truncating to singly occupied states before squaring loses all exchange.
    assert (P.T@T@P).nnz == 0 and (P.T@T@T@P).nnz > 0


@pytest.mark.parametrize('sides,c,U,t', [((2,1,1),2.,16.,.7),((3,1,1),1.,16.,.5),((2,2,1),1.,32.,.4)])
def test_full_car_low_band_energies_and_complete_dynamics(sides,c,U,t):
    n, edges = algebra.geometry(sides)
    masks,T,D = algebra.operators(n,edges)
    P = algebra.spin_embedding(n,masks).toarray()
    H = c*T.toarray()/2+U*np.diag(D)
    Heff = c*c/(4*U)*algebra.exchange_four(n,edges).toarray()
    energies,V = eigh(H)
    m = 1<<n
    tau = c*len(edges)
    assert energies[m] >= U-tau-1e-10
    assert energies[m-1] < 1e-10
    assert max(abs(energies[:m]-eigh(Heff,eigvals_only=True))) <= tau**3/(U-tau)**2+1e-10
    actual = np.linalg.norm(expm(-1j*t*H)@P-P@expm(-1j*t*Heff),2)
    assert actual <= 4*tau/U+3*t*tau**3/U**2+1e-10
    # Canonical polar isometry is compared on the same spin carrier.
    spectral = V[:,:m]@V[:,:m].T
    gram = P.T@spectral@P
    lam,Q = eigh(gram)
    W = spectral@P@(Q@np.diag(lam**(-.5))@Q.T)
    assert np.linalg.norm(W-P,2) <= 2*tau/U+1e-10
    assert np.linalg.norm(W.T@H@W-Heff,2) <= 3*tau**3/U**2+1e-10


def test_all_charge_sectors_separated_from_half_filled_ground():
    n,edges=algebra.geometry((3,1,1));c=1;U=16
    sector_min=[]
    for N in range(2*n+1):
        masks,T,D=algebra.operators(n,edges,N)
        K=c*T.toarray()/2+np.diag(U*D-U*N/2+U*n/2)
        sector_min.append(eigh(K,eigvals_only=True)[0])
    assert sector_min[n] <= 0
    assert min(e for N,e in enumerate(sector_min) if N!=n)-sector_min[n] >= U/2-c*len(edges)


def test_actual_dimer_source_entangles_and_wrong_exchange_sign_fails():
    n,edges=algebra.geometry((2,1,1));masks,T,D=algebra.operators(n,edges)
    P=algebra.spin_embedding(n,masks).toarray(); c=128;U=c*c
    H=c*T.toarray()/2+U*np.diag(D)
    K=algebra.exchange_four(n,edges).toarray()/4
    initial=np.eye(4)[:,1];t=np.pi/2
    ideal=expm(-1j*t*K)@initial
    reduced=ideal.reshape(2,2)@ideal.reshape(2,2).conj().T
    assert np.allclose(reduced,np.eye(2)/2,atol=1e-14)
    psi=expm(-1j*t*H)@P@initial
    delta=np.linalg.norm(psi-P@ideal)
    assert delta < 5/64
    assert abs(np.vdot(P@ideal,psi))**2 > .5
    wrong=P@expm(+1j*t*K)@initial
    assert abs(np.vdot(wrong,P@ideal)) < 1e-13
    assert algebra.interaction_audit()['adjacent_bond_commutator_nonzero_entries'] > 0


@pytest.mark.parametrize('edges', [((0,0,1),),((0,1,0),),((0,1,1),(0,1,-1))])
def test_invalid_or_duplicate_source_edges_rejected(edges):
    with pytest.raises(ValueError):algebra.operators(2,edges)
