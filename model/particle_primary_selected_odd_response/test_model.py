"""Full-CAR, gauge, original-state and same-clock odd-response checks."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as s
from scipy.linalg import expm

from . import model as m


def test_exact_six_mode_residual_on_every_occupation():
    control=m.exact_controls()
    R=s.zeros(64)
    for i in range(3):
        row=[0,1,2];row[i]=3+i
        R+=m.creation_word(row,6)
    G=R.H*R
    assert G.eigenvals()=={s.Integer(0):50,s.Integer(1):6,s.Integer(2):6,s.Integer(3):2}
    assert control['universal_residual_squared_norm']==3
    assert control['native_to_coordinate_basis_permutation']==[0,1,4,2,5,7,10,3,6,8,11,9,12,13,14]


def test_all_four_spin_monomials_are_actual_gauge_singlets():
    native=m.sectors.native
    blocks=dict(m.sectors.native_blocks());weak=blocks[-3];e=blocks[6][0]
    for alpha,b in product(range(2),repeat=2):
        modes=[15*alpha+i for i in weak]+[15*b+e]
        mask=sum(1<<i for i in modes)
        for generator in native.native_currents():
            terms=[(15*r+i,15*r+j,c) for r in range(4)
                   for (i,j),c in generator.todok().items() if c]
            assert m.sectors._act({mask:s.S.One},terms)=={}


def test_literal_Wilson_hopping_images_are_orthogonal_with_arbitrary_links():
    response=m.sectors.branch.response;density=response.density
    pairs=response.beta_pairs();plus=s.Matrix.hstack(pairs[0][0],pairs[1][0])
    g=density.gamma_data();beta=g[3]
    permutation=m.exact_controls()['native_to_coordinate_basis_permutation']
    blocks=dict(m.sectors.native_blocks());indices=list(blocks[-3])+list(blocks[6])
    columns=s.zeros(60,3)
    for j,index in enumerate(indices):
        for spin in range(4):columns[15*spin+index,j]=plus[spin,0 if j<2 else 1]
    images=[]
    for axis,sign in product(range(3),(1,-1)):
        c=s.Rational(3,5);d=s.Rational(4,5)
        link=density.native_weak_link(axis,c,d).extract(permutation,permutation)
        hop=(-beta+beta*g[axis])/2
        h=s.kronecker_product(hop if sign==1 else hop.H,link if sign==1 else link.H)
        images.append(h*columns)
    full=s.Matrix.vstack(*images)
    assert density.clean(full.H*full)==3*s.eye(3)


def test_success_is_unconditional_and_vacancy_bound_is_operator_positive():
    probes=[]
    for alpha,b in product(range(2),repeat=2):
        probes.append(np.array(m.creation_word((2*alpha,2*alpha+1,4+b),6),complex))
    success=sum(A.conj().T@A for A in probes)/4
    occupation=np.diag([mask.bit_count()/2 for mask in range(64)])
    assert np.linalg.eigvalsh(success-np.eye(64)+occupation).min()>-1e-13
    assert np.linalg.norm(success,2)==1
    rho=np.diag(np.arange(1,65)/2080)
    assert sum(np.trace(rho@(A.conj().T@A+np.eye(64)-A.conj().T@A)).real for A in probes)/4==pytest.approx(1)


def test_even_exchange_retains_failure_memory_and_reference():
    A=np.array(m.creation_word((0,1,2),6),complex)
    B=np.array(m.creation_word((3,4,5),6),complex)
    X=A@B.conj().T;U=expm(np.pi/2*(X-X.conj().T))
    assert np.max(abs(U.conj().T@U-np.eye(64)))<1e-14
    numbers=np.diag([mask.bit_count() for mask in range(64)])
    assert np.max(abs(U@numbers-numbers@U))<1e-14
    # Start the apparatus triple full. Its empty/full outputs form the
    # complete creation instrument, without resetting the source.
    success=U[np.ix_(range(8),range(56,64))]
    failure=U[np.ix_(range(56,64),range(56,64))]
    local=np.array(m.creation_word((0,1,2),3),complex)
    assert np.max(abs(success-local))<1e-14
    assert np.max(abs(failure-(np.eye(8)-local.conj().T@local)))<1e-14
    rho=np.diag(np.arange(1,9)/36);purifier=np.sqrt(rho)
    branches=[K@purifier for K in (success,failure)]
    assert sum(np.vdot(v,v).real for v in branches)==pytest.approx(1)
    assert np.max(abs(sum(v.T@v.conj() for v in branches)-rho))<1e-14


def test_uniform_selected_density_uses_full_Fock_pivots():
    value=m.weak_e_allocation()
    assert value['retained_fractions']==(F(),F(21,50),F(),F(),F(21,50))
    assert value['nonconstant_link_margin_lower']>0
    assert value['uniform_weak_plus_e_defect_density_upper']<F(341,100)
    q6=value['native_rows'][-1]
    assert q6['all_positive_pivots'] and q6['full_Fock_dimension']==256
    assert m.refinement.q6_positive_shift(q6['onsite_per_endpoint'],q6['inverse_electric_margin'],q6['occupation_loss_upper'])


@pytest.fixture(scope='module')
def receipts():
    return {n:m.bind_member(m.contract.member(n)) for n in (0,1,2,8)}


def test_actual_original_weight_and_capacity_scale(receipts):
    for n,bound in receipts.items():
        r=bound['result'];a=F(1,2**n)
        assert F(r['original_creation_success_lower'])>=F(59,400)
        assert F(r['original_transfer_weight_lower'])>=F(1,50)
        assert F(r['finite_record_original_weight_lower'])>=F(1,64)
        assert list(map(F,r['original_transfer_interval']))==[1/(2*a),35/(2*a)]
        assert F(r['centered_second_moment_upper'])==(F(9,2) if n==0 else 9/a**2)
        assert r['source_native_increment']==[0,2,0,0,1]
        assert r['charge_conjugate_source_increment']==[0,-2,0,0,-1]
        assert not r['source_marginal_replaced_by_trial']
        assert not r['division_by_success_probability']
        assert r['all_failures_and_quantum_reference_retained']
        assert not r['all_fixed_energy_weight_vanishes_proved']
        assert not r['critical_density_Weyl_survival_decided']
        assert r['physical_mass_table_GeV']==[]
    zero=receipts[0]['result']['member_zero_sharper_band']
    assert F(zero['original_transfer_weight_lower'])>F(23,32)
    assert F(zero['final_excitation_transport']['original_final_excitation_weight_lower'])>F(21,32)


def test_original_chebyshev_bound_on_an_independent_stationary_model():
    rng=np.random.default_rng(114)
    H=np.diag(np.arange(8,dtype=float));rho=np.diag(np.exp(-np.arange(8)))
    rho/=rho.trace()
    A=rng.normal(size=(8,8))+1j*rng.normal(size=(8,8));A/=np.linalg.norm(A,2)
    center=2.;radius=1.5
    weights=abs(A)**2*np.diag(rho)[None,:]
    transfer=np.diag(H)[:,None]-np.diag(H)[None,:]
    R=H@A-A@H-center*A
    second=np.trace(rho@R.conj().T@R).real
    assert abs(second-np.sum((transfer-center)**2*weights))<1e-12
    assert weights[abs(transfer-center)<=radius].sum()>=weights.sum()-second/radius**2-1e-12


def test_actual_source_and_forged_mass_claim_are_rejected(receipts):
    member=m.contract.member(0);m.verify_bound(member,receipts[0])
    changed=deepcopy(member);changed['clock']['Schrodinger']='exp(-2itH)'
    with pytest.raises(ValueError):m.bind_member(changed)
    wrong=deepcopy(receipts[0]);wrong['result']['all_fixed_energy_weight_vanishes_proved']=True
    wrong['result_sha256']=m.contract.digest(wrong['result'])
    with pytest.raises(ValueError):m.verify_bound(member,wrong)
    wrong=deepcopy(receipts[0]);wrong['result']['physical_mass_table_GeV']=[{'electron':'0.000511'}]
    wrong['result_sha256']=m.contract.digest(wrong['result'])
    with pytest.raises(ValueError):m.verify_bound(member,wrong)
