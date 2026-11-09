"""Full-residual forms, complete preparation and retained quantum records."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as s

from . import model as m


@pytest.fixture(scope='module')
def inputs():
    member=m.contract.member(0)
    inherited=m.odd.branch.bind_member(member)
    return member,inherited['result']['odd_endpoint_trial']


def test_arbitrary_weight_formula_reproduces_inherited_full_residual(inputs):
    member,trial=inputs
    actual=m.residual_form(trial,trial['rational_hopping_coefficients'],trial['rational_loop_coefficients'])
    old=m.odd.refined_bound(member)
    assert F(actual['source_energy_center'])==F(old['source_trial_energy'])
    assert abs(float(actual['full_residual_squared_upper']-F(old['source_residual_squared_upper'])))<1e-25
    improved=m.residual_form(trial,m.HOPPING_WEIGHTS,m.LOOP_WEIGHTS)
    assert improved['source_energy_center']>actual['source_energy_center']
    assert improved['full_residual_squared_upper']<actual['full_residual_squared_upper']
    assert improved['full_residual_squared_upper']<F(49,2)


@pytest.fixture(scope='module')
def source_moment_matrices(inputs):
    # Independent assembly of first/second action matrices on normalized
    # packets. The second moment includes an orthogonal complementary space.
    _,trial=inputs
    n=np.array([float(F(x['norm_squared'])) for x in trial['hopping_packets']]+[float(F(trial['loop_packet_norm_squared']))]*5)
    d=np.array([float(F(x['relative_H0_energy'])) for x in trial['hopping_packets']])
    A=np.array([[float(F(x)) for x in row] for row in trial['loop_packet_matrix']])
    C=np.zeros((11,11));C[6:,6:]=np.array(m.odd.variational.character_control()['cubic_loop_singlet_matrix'],float)/(-420)
    D=np.zeros((11,11));D[:6,:6]=np.diag(d);D[6:,6:]=A
    D0=D-C;gram=np.zeros((11,11))
    basis=[tuple(F(i==j) for i in range(11)) for j in range(11)]
    def square(t):return float(m.odd.odd_parent_moment(t[:6],t[6:])['exact_parent_second_action_norm_squared'])
    for i in range(11):gram[i,i]=square(basis[i])
    for i,j in product(range(11),repeat=2):
        if i<j:gram[i,j]=gram[j,i]=(square(tuple(x+y for x,y in zip(basis[i],basis[j])))-gram[i,i]-gram[j,j])/2
    first=np.zeros((12,12));first[0,1:]=n;first[1:,0]=n;first[1:,1:]=np.diag(n)@D
    second=np.zeros((12,12));second[0,0]=sum(n)
    second[0,1:]=n@D;second[1:,0]=second[0,1:]
    second[1:,1:]=np.diag(n)@(D0@D0+D0@C+C@D0)+gram
    norms=np.sqrt(np.r_[1,n])
    first=first/norms[:,None]/norms[None,:];second=second/norms[:,None]/norms[None,:]
    return n,first,second


def test_full_complement_extension_checks_residual_and_probability(inputs,source_moment_matrices):
    _,trial=inputs;n,first,second=source_moment_matrices
    outside=second-first@first
    eigen,rotation=np.linalg.eigh(outside)
    assert eigen.min()>-1e-10 and eigen.max()>200
    coupling=(rotation*np.sqrt(np.maximum(eigen,0)))@rotation.T
    H=np.block([[first,coupling],[coupling,np.zeros((12,12))]])
    t=np.array([float(x) for x in m.HOPPING_WEIGHTS+m.LOOP_WEIGHTS])
    vector=np.r_[1,-t*np.sqrt(n),np.zeros(12)];vector/=np.linalg.norm(vector)
    e=vector@H@vector;residual=H@vector-e*vector
    certified=m.residual_form(trial,m.HOPPING_WEIGHTS,m.LOOP_WEIGHTS)
    assert abs(e-float(certified['relative_Rayleigh_energy']))<1e-12
    assert abs(residual@residual-float(certified['full_residual_squared_upper']))<1e-10
    spectrum,modes=np.linalg.eigh(H);weight=np.abs(modes.T@vector)**2
    inside=np.abs(spectrum-e)<=23/4
    assert weight[inside].sum()>=1-float(certified['full_residual_squared_upper'])/(23/4)**2-1e-12
    # Keeping only the12-dimensional compression would omit actual norm.
    compressed=vector[:12];compressed_residual=first@compressed-e*compressed
    assert residual@residual>compressed_residual@compressed_residual+20


def test_residual_form_on_independent_rational_weight_choices(inputs,source_moment_matrices):
    _,trial=inputs;n,first,second=source_moment_matrices
    for scale in (F(0),F(1,2),F(1),F(3,2),F(-1)):
        th=tuple(scale*x for x in m.HOPPING_WEIGHTS);tl=tuple(scale*x for x in m.LOOP_WEIGHTS)
        exact=m.residual_form(trial,th,tl)
        v=np.r_[1,-np.array([float(x) for x in th+tl])*np.sqrt(n)];v/=np.linalg.norm(v)
        e=v@first@v;variance=v@second@v-e*e
        assert abs(e-float(exact['relative_Rayleigh_energy']))<1e-11
        assert abs(variance-float(exact['full_residual_squared_upper']))<1e-9


def test_actual_packet_householder_preserves_norm_and_parity(inputs):
    _,trial=inputs
    data=m.residual_form(trial,m.HOPPING_WEIGHTS,m.LOOP_WEIGHTS)
    n=data['packet_norms_squared'];t=m.HOPPING_WEIGHTS+m.LOOP_WEIGHTS
    psi=s.Matrix([0,1]+[-s.Rational(x.numerator,x.denominator)*s.sqrt(s.Rational(y.numerator,y.denominator)) for x,y in zip(t,n)])
    psi=psi/s.sqrt(s.Rational(data['trial_norm_squared'].numerator,data['trial_norm_squared'].denominator))
    blank=s.eye(13)[:,0];v=blank-psi
    assert s.simplify((psi.H*psi)[0])==1 and (blank.H*psi)[0]==0
    assert s.simplify((v.H*v)[0])==2
    assert s.simplify(blank-v*(v.H*blank)[0])==psi
    spatial=s.diag(1,*([-1]*12))
    assert spatial*psi==-psi and spatial*blank==blank
    # Both actual packet and blank states are in the even CAR block. Spatial
    # inversion is deliberately changed by this intervention.
    assert s.eye(13)*psi==psi


def test_complete_graded_exchange_retains_arbitrary_reference():
    control=m.swap_reference_control()
    assert control['exact_recovery_of_entangled_source_reference']
    d=4;parity=np.array([0,0,1,1]);swap=np.zeros((16,16))
    for i,j in product(range(d),repeat=2):swap[j*d+i,i*d+j]=(-1.)**(parity[i]*parity[j])
    assert np.array_equal(swap.T@swap,np.eye(16))
    target=np.array([3.,4.,0.,0.])/5
    source_ref=np.arange(1,13)+1j*np.arange(12,0,-1);source_ref=source_ref/np.linalg.norm(source_ref)
    initial=np.einsum('ir,j->ijr',source_ref.reshape(4,3),target).reshape(-1)
    final=np.kron(swap,np.eye(3))@initial
    expected=np.einsum('i,jr->ijr',target,source_ref.reshape(4,3)).reshape(-1)
    assert np.linalg.norm(final-expected)<1e-14
    tensor=final.reshape(4,12)
    assert np.linalg.norm(tensor.conj().T@tensor-np.outer(source_ref.conj(),source_ref))<1e-14
    # Ordinary swap on odd-odd amplitudes has a physically different sign.
    assert swap[2*4+2,2*4+2]==-1


def test_complete_clock_on_prepared_nonstationary_source_keeps_memory():
    source=s.Matrix([s.Rational(3,5),s.Rational(4,5)])
    memory_ref=s.Matrix([1,0,0,s.I])/s.sqrt(2)
    U=s.diag(1,s.I);rows=m.spectroscopy.energy_instrument(U,4)
    joint=s.kronecker_product(source,memory_ref)
    probabilities=[]
    for row in rows:
        after=s.kronecker_product(row,s.eye(4))*joint
        probabilities.append(s.simplify((after.H*after)[0]))
    assert probabilities==[s.Rational(9,25),0,0,s.Rational(16,25)]
    assert sum(probabilities)==1
    assert U*(source*source.H)!=(source*source.H)*U


@pytest.fixture(scope='module')
def bound(inputs):return inputs[0],m.bind_member(inputs[0])


def test_selected_source_preparation_intervention_and_clock_bounds(bound):
    member,receipt=bound;r=receipt['result'];t=r['physical_trial']
    assert F(t['full_residual_upper'])<F('4.947619')
    assert F(t['some_physical_A1u_excitation_interval'][0])>1
    assert F(t['some_physical_A1u_excitation_interval'][1])<15
    assert F(r['original_prepared_spectral_weight_lower'])>F(137,529)
    assert F(r['reported_original_probability_lower'])>F(1,4)
    assert F(r['reported_interval_lower_minus_ground_upper'])>F(1,5)
    assert F(r['source_work_mean_interval'][0])>5
    assert r['additional_memory_qubit_bound']==696
    assert r['initial_state_is_selected_primary_preparation']
    assert not r['source_marginal_preserved_by_preparation_intervention']
    assert r['memory_retains_displaced_source_and_its_original_reference_correlations']
    assert not r['postselection_or_success_normalization']
    assert not r['this_is_the_selected_K00_response']
    assert not r['physical_species_or_relativistic_rest_mass_proved']
    assert r['physical_mass_table_GeV']==[]
    m.verify_bound(member,receipt)


def test_wrong_source_or_forged_prepared_weight_is_rejected(bound):
    member,receipt=bound
    with pytest.raises(ValueError):m.bind_member(m.contract.member(1))
    forged=deepcopy(receipt);forged['result']['reported_original_probability_lower']='1'
    forged['result_sha256']=m.contract.digest(forged['result'])
    with pytest.raises(ValueError):m.verify_bound(member,forged)
    altered=deepcopy(member);altered['clock']['Schrodinger']='exp(-2itH)'
    with pytest.raises(ValueError):m.verify_bound(altered,receipt)
