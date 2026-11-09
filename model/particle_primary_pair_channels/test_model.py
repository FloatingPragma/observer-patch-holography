"""Actual Clifford/CAR channels and original spectral-law controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product

import numpy as np
import pytest
import sympy as s

from . import model as m


def test_exact_full_fock_pauli_effects_and_all_defect_sectors():
    probes,D=m.local_controls()
    for probe in probes:
        effect=probe.H*probe
        assert set(effect.eigenvals())=={0,s.Rational(1,2),1}
        for defects in range(5):
            indices=[i for i in range(16) if D[i,i]==defects]
            block=effect.extract(indices,indices)
            assert min(block.eigenvals())>=max(0,1-s.Rational(defects,2))
        assert probe*probe!=s.zeros(16)  # A Pauli combination is not a rank-one CAR pair.
        assert probe**3==s.zeros(16)


def test_physical_rotation_representation_and_inversion_are_exact():
    rotations,inversion=m.coefficient_representation()
    matrices=[np.array(r,dtype=np.int64) for _,r in rotations]
    keys={tuple(r.flat) for r in matrices}
    for a in matrices:
        assert np.array_equal(a.T@a,np.eye(32,dtype=np.int64))
        assert tuple((np.array(inversion,dtype=np.int64)@a).flat) not in keys
        for b in matrices:assert tuple((a@b).flat) in keys
    assert inversion*inversion==s.eye(32)
    for _,r in rotations:assert r*inversion==inversion*r


def test_cubic_parity_ranks_and_independent_adapted_basis_car_contraction():
    p=m.cubic_projectors();W,labels=m.coefficient_transform('cubic')
    expected={name+parity:rank for parity in ('g','u')
              for name,rank in (('A1',1),('A2',1),('E',2),('T1',6),('T2',6))}
    assert {name:int(v.trace()) for name,v in p.items()}==expected
    probes,D=m.local_controls()
    for name,rank in expected.items():
        selected=[j for j,label in enumerate(labels) if label[0]==name]
        assert W.extract(selected,range(32)).H*W.extract(selected,range(32))==p[name]
        local=[sum((W[j,a]*probes[a] for a in range(4)),s.zeros(16)) for j in selected]
        comm=sum((b.H*b-b*b.H for b in local),s.zeros(16))
        assert m.response.apparatus.clean(comm)==s.Rational(rank,16)*(2*s.eye(16)-D)


def test_walsh_zero_mode_and_parity_are_operator_channel_labels():
    W,labels=m.coefficient_transform('walsh');_,inversion=m.coefficient_representation()
    for j,(k,a) in enumerate(labels):
        assert W[j,:]*inversion==(-1)**(sum(k)+1)*W[j,:]
    scalar=W[0,:].T
    assert scalar*scalar.T==m.cubic_projectors()['A1u']
    vectors=sum((W[j,:].T*W[j,:] for j in (1,2,3)),s.zeros(32))
    assert m.cubic_projectors()['T1u']*vectors==vectors
    assert vectors!=m.cubic_projectors()['T1u']  # Multiplicity two is retained.


def test_complete_coherent_channel_control_and_entangled_reference():
    probes,_=m.local_controls();N=4;dimension=16
    # One actual local quartet executes both branches of a coefficient transform.
    W=s.Matrix([[(-1)**((i&j).bit_count())/2 for j in range(4)] for i in range(4)])
    W=W.applyfunc(lambda x:s.Rational(str(x)))
    complements=[]
    for b in probes:
        effect=b.H*b
        root=s.eye(dimension)+(2*s.sqrt(2)-3)*effect+(2-2*s.sqrt(2))*effect**2
        assert m.response.apparatus.clean(root*root)==s.eye(dimension)-effect
        complements.append(root)
    kraus=[m.response.apparatus.clean(sum((W[j,i]*row[i] for i in range(N)),s.zeros(dimension))/s.sqrt(N))
           for row in (probes,complements) for j in range(N)]
    assert m.response.apparatus.clean(sum((k.H*k for k in kraus),s.zeros(dimension)))==s.eye(dimension)
    psi=s.zeros(32,1);psi[24]=s.sqrt(s.Rational(2,3));psi[7]=s.sqrt(s.Rational(1,3))
    records=[s.kronecker_product(k,s.eye(2))*psi for k in kraus]
    assert s.simplify(sum((v.H*v)[0] for v in records))==1
    success=sum((v.H*v)[0] for v in records[:N])
    incoherent=sum((psi.H*s.kronecker_product(b.H*b,s.eye(2))*psi)[0] for b in probes)/N
    assert s.simplify(success-incoherent)==0
    assert any(m.response.apparatus.clean(kraus[j].H*kraus[j]-probes[j].H*probes[j]/N)!=s.zeros(dimension)
               for j in range(N))


def test_polynomial_log_slope_and_signed_commutator_window_on_actual_measures():
    energies=np.array([0.,.02,.3,2.,6.,11.,45.]);beta=2.;order=5;scale=128
    x=beta*(energies.max()-energies)/(2*scale)
    t=np.array([sum(y**j/float(s.factorial(j)) for j in range(order+1)) for y in x])
    tm=np.array([sum(y**j/float(s.factorial(j)) for j in range(order)) for y in x])
    assert np.all(beta*tm/t<=beta+1e-15)
    p=t**(2*scale);p=p/p.sum()
    B=np.zeros((7,7));B[4,0]=.8;B[1,0]=.3;B[0,4]=.1;B[6,0]=.05
    mu=[(energies[j]-energies[i],p[i]*B[j,i]**2) for i,j in product(range(7),repeat=2)]
    chi=np.trace(np.diag(p)@(B.T@B-B@B.T))
    direct=sum((p[i]-p[j])*B[j,i]**2 for i,j in product(range(7),repeat=2))
    assert abs(chi-direct)<1e-15 and chi>0
    moment=sum(weight*(omega-6)**2 for omega,weight in mu)
    result=m.passive_window(str(chi),2,6,str(moment),'1/100',40)
    actual=sum(weight for omega,weight in mu if .01<=omega<=40)
    assert actual+1e-14>=float(result['original_transfer_weight_lower'])
    assert sum(weight for omega,weight in mu if omega<0)>0
    # A population inversion invalidates the passive sign, even if stationary.
    inverted=p[::-1]
    assert any(energies[j]<energies[i] and inverted[i]-inverted[j]>0
               for i,j in product(range(7),repeat=2))


@pytest.fixture(scope='module')
def bound():
    member=m.contract.member(0)
    return member,m.bind_member(member)


def test_actual_selected_channel_weights_clocks_and_mass_boundary(bound):
    member,receipt=bound;m.verify_bound(member,receipt);row=receipt['result']
    assert F(row['spatially_averaged_each_Pauli_transfer_weight_simple_lower'])==F(5,16)
    assert F(row['each_Walsh_pair_record_original_success_lower'])>F(11,512)
    assert len(row['cubic_parity_channels'])==10
    for channel in row['cubic_parity_channels']:
        rank=channel['isotypic_rank']
        assert F(channel['original_transfer_weight_lower'])>=F(11*rank,1024)
        assert F(channel['finite_reported_original_weight_lower'])>=F(11*rank-1,1024)
        assert F(channel['original_final_excitation_weight_simple_lower'])==F(11*rank-8,1024)
        assert channel['weight_is_isotypic_aggregate']
    assert row['absolute_seconds_per_source_tick'] is None
    assert row['physical_mass_table_GeV']==[]
    assert not row['Fourier_label_is_conserved_particle_momentum']
    assert not row['full_primary_spectral_matrix_or_histogram_executed']
    assert row['reference_idle_and_retained']
    assert not row['uniform_refinement_weight_or_fixed_lower_endpoint_proved']


def test_source_level_and_result_mutations_are_rejected(bound):
    member,receipt=bound
    altered=deepcopy(receipt);altered['result']['cubic_parity_channels'][0]['original_transfer_weight_lower']='1'
    with pytest.raises(ValueError):m.verify_bound(member,altered)
    altered=deepcopy(member);altered['clock']['GeV_scale']='6'
    with pytest.raises(ValueError):m.bind_member(altered)
    with pytest.raises(ValueError):m.bind_member(m.contract.member(1))
    with pytest.raises(ValueError):m.coefficient_transform('particle momentum')
    with pytest.raises(ValueError):m.passive_window(1,2,6,3,1,5)
    with pytest.raises(ValueError):m.passive_window(.5,2,6,3,1,40)


def test_all_source_evidence_paths_exist():
    assert all((m.contract.ROOT/name).is_file() for name in m.SOURCE_INPUTS)
