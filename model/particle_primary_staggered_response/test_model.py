"""Actual phased cube blocks, CAR sectors and complete original records."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

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


def test_actual_cube_with_nontrivial_link_phases_has_exact_cancellation():
    pairs=m.odd.response.beta_pairs();chart=s.Matrix.hstack(pairs[0][0],pairs[1][0],pairs[0][1],pairs[1][1])
    gamma=m.odd.response.density.gamma_data();beta=gamma[3]
    spin=[np.array(chart.H*((-beta+beta*g)/2)*chart,dtype=complex) for g in gamma[:3]]
    sites=list(product(range(2),repeat=3));ix={x:i for i,x in enumerate(sites)}
    h=np.zeros((32,32),complex);R=np.zeros_like(h);D=np.zeros_like(h)
    for x in sites:
        n=4*ix[x];w=(-1)**sum(x)
        R[n:n+2,n+2:n+4]=w*np.eye(2)
        D[n:n+2,n:n+2]=np.eye(2);D[n+2:n+4,n+2:n+4]=-np.eye(2)
        for axis in range(3):
            if x[axis]:continue
            y=list(x);y[axis]=1;y=tuple(y);q=4*ix[y]
            phase=np.exp(1j*(ix[x]+2*axis+1)/3)
            h[n:n+4,q:q+4]=phase*spin[axis]
            h[q:q+4,n:n+4]=phase.conjugate()*spin[axis].conj().T
    positive=[4*x+j for x in range(8) for j in (0,1)]
    negative=[4*x+j for x in range(8) for j in (2,3)]
    comm=h@R-R@h
    assert np.linalg.norm(comm[np.ix_(positive,negative)])<1e-14
    assert np.linalg.norm(comm[np.ix_(negative,positive)])<1e-14
    J=1j*(R-R.conj().T)
    assert np.linalg.norm(h@J-J@h)<1e-14
    assert np.linalg.norm(h@(R+R.conj().T)-(R+R.conj().T)@h)>1
    number=np.block([[comm[np.ix_(positive,positive)],np.zeros((16,16))],
                     [np.zeros((16,16)),-comm[np.ix_(negative,negative)].T]])
    assert abs(np.trace(number))<1e-14
    assert np.linalg.norm(number,2)<=1.5+1e-14
    assert np.linalg.norm(3*D@R-R@(3*D)-6*R)<1e-14


def test_full_occupation_pair_commutator_and_Hermitian_first_moment():
    cs=[car(4,i) for i in range(4)];n=[c.T@c for c in cs]
    T=cs[0].T@cs[1]-cs[2].T@cs[3]
    D=n[0]+(np.eye(16)-n[1])+n[2]+(np.eye(16)-n[3])
    assert np.linalg.norm(T.T@T-T@T.T-(2*np.eye(16)-D))<1e-14
    H=3*(n[0]-n[1]+n[2]-n[3]);A=1j*(T-T.T)
    bracket=A@(H@A-A@H)-(H@A-A@H)@A
    assert np.linalg.norm(bracket-12*(2*np.eye(16)-D))<1e-14
    rho=np.diag(np.exp(-np.diag(H)));rho/=np.trace(rho)
    energies=np.diag(H);law=(rho.diagonal()[None,:]*abs(A)**2)
    mean=np.sum(law*(energies[:,None]-energies[None,:]))
    assert abs(mean-np.trace(rho@bracket)/2)<1e-13


def test_number_sector_bound_and_actual_quasispin_norm():
    for q in range(33):assert min(q,32-q)**2<=16*q
    # Sixteen equal pair ladders have maximum spin8. Retaining all
    # magnetizations verifies the full-Fock cap, not just its vacuum value.
    ladder=[(8-z)*(8+z+1) for z in range(-8,9)]
    assert max(ladder)==72 and F(72,512)==F(9,64)
    for spin2 in range(17):
        j=F(spin2,2)
        for k in range(spin2+1):
            z=-j+k
            assert (j-z)*(j+z+1)<=72


def test_cube_rotation_and_inversion_are_operator_labels():
    c=m.clifford_control()
    assert c['proper_cube_rotation_count']==24
    assert c['combined_operator_type']=='A2g'
    assert c['mask_inversion_parity']*c['local_scalar_inversion_parity']==1
    assert c['CAR_parity']=='even'
    assert not c['mask_is_conserved_physical_momentum']


@pytest.fixture(scope='module')
def bound():
    member=m.contract.member(0)
    return member,m.bind_member(member)


def test_unchanged_selected_state_has_sharper_original_channel_band(bound):
    member,value=bound;r=value['result']
    assert F(r['selected_centered_second_moment_upper'])<F(207,1600)
    assert F(r['original_transfer_weight_lower'])>F(17,1280)
    assert r['original_transfer_band']==[3,9]
    assert F(r['final_excitation']['original_final_excitation_weight_lower'])>F(23,5120)
    assert list(map(F,r['final_excitation']['inferred_final_excitation_interval']))==[3,10]
    assert F(r['Hermitian_quadrature_signed_energy_mean_lower'])>F(531,1600)
    assert F(r['complete_clock_event_loss_upper'])<=F(1,8192)
    assert F(r['reported_transfer_original_weight_lower'])>F(539,40960)
    assert r['reported_transfer_interval']==['95/32','289/32']
    assert r['all_failures_and_reference_retained'] and not r['selected_state_changed']
    assert not r['final_state_irrep_fixed_by_operator']
    assert not r['dispersion_or_rest_mass_identified'] and r['physical_mass_table_GeV']==[]
    m.verify_bound(member,value)


def test_source_and_rehashed_response_forgery_are_rejected(bound):
    member,value=bound;forged=deepcopy(value)
    forged['result']['original_transfer_weight_lower']='1'
    forged['result_sha256']=m.contract.digest(forged['result'])
    with pytest.raises(ValueError):m.verify_bound(member,forged)
    with pytest.raises(ValueError):m.bind_member(m.contract.member(1))
