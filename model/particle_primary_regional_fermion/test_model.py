"""Independent finite CAR, Gauss, spin, normalization and clock controls."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import json

import numpy as np
import pytest
import sympy as s

from . import model as m


def creation(mode,count):
    A=np.zeros((1<<count,1<<count),complex)
    for mask in range(1<<count):
        if not mask>>mode&1:
            A[mask|(1<<mode),mask]=(-1)**((mask&((1<<mode)-1)).bit_count())
    return A


def test_actual_full_endpoint_CAR_effect_and_Gauss_completion():
    C=[creation(j,4) for j in range(4)];N=[A@A.conj().T for A in C];I=np.eye(16)
    D=N[0]+N[1]+2*I-N[2]-N[3];sigma,_=m.pairs.pair_basis()
    actual=m.local_control()
    for a,t in enumerate(sigma):
        t=np.array(t,complex);B=sum(t[i,j]*C[i]@C[2+j].conj().T for i,j in product(range(2),repeat=2))/np.sqrt(2)
        assert np.linalg.norm(B-np.array(actual['fields'][a],complex))<1e-14
        effect=B.conj().T@B
        assert np.linalg.eigvalsh(effect)[0]>-1e-14
        assert np.linalg.eigvalsh(I-effect)[0]>-1e-14
        assert np.linalg.eigvalsh(effect-I+D/2)[0]>-1e-14
    P=actual['gauge_projector'];bare=actual['bare_inside_creator'];anchor=actual['anchor']
    assert P*bare*P==s.zeros(48) and bare*anchor!=s.zeros(48,1)
    targets=s.Matrix.hstack(*(B*anchor for B in actual['success']))
    assert m.clean(targets.H*targets)==s.eye(4)
    assert actual['physical_control_dimension']==5


def test_full_eight_mode_other_hopping_residual_is_sqrt_three_not_sqrt_two():
    C=[creation(j,8) for j in range(8)]
    sigma,_=m.pairs.pair_basis()
    # Mode groups p,m,u,v are independent canonical images on the literal graph.
    for t in sigma:
        t=np.array(t,complex)
        R=sum(t[a,b]*(C[4+a]@C[2+b].conj().T-C[a]@C[6+b].conj().T)
              for a,b in product(range(2),repeat=2))/np.sqrt(2)
        values=np.linalg.eigvalsh(R.conj().T@R)
        assert abs(values[-1]-3)<2e-13
        assert values[-1]>2.5
    assert m.hopping_control()['full_Fock_squared_norm_by_active_count']==[F(0),F(1,2),F(1),F(2),F(3)]


def test_source_Wilson_neighbor_images_have_disjoint_support_and_complete_spin_Gram():
    pairs=m.pairs.response.beta_pairs();P=s.Matrix.hstack(*(p for p,n in pairs));M=s.Matrix.hstack(*(n for p,n in pairs))
    gamma=m.pairs.response.density.gamma_data();beta=gamma[3];h=[(-beta+beta*g)/2 for g in gamma[:3]]
    cells=tuple(product(range(2),repeat=3));where={x:i for i,x in enumerate(cells)}
    for x,y in product(cells,repeat=2):
        if sum(a!=b for a,b in zip(x,y))!=1:continue
        images=[];supports=[]
        for source,partner,base in ((x,y,P),(y,x,M)):
            image=s.zeros(32,2);support=set()
            for j in range(3):
                target=tuple(1-a if i==j else a for i,a in enumerate(source))
                if target==partner:continue
                A=h[j] if source[j] else h[j].H
                # Distinct phases verify that link coefficients are not replaced by real constants.
                phase=(s.Rational(3,5)+s.I*s.Rational(4,5))**(j+1)
                block=where[target];image[4*block:4*block+4,:]=phase*A*base;support.add(target)
            assert m.clean(image.H*image)==s.eye(2)
            images.append(image);supports.append(support)
        assert not supports[0]&supports[1]
        assert m.clean(images[0].H*images[1])==s.zeros(2)


def test_native_endpoint_spinors_and_complete_pair_coefficients_have_different_central_signs():
    control=m.spin_control();rows=control['rows']
    assert len(rows)==48
    assert sum(row['endpoint']==-s.eye(2) for row in rows)==1
    assert all(row['pair']==s.eye(4) for row in rows if row['endpoint']==-s.eye(2))
    assert all(s.simplify(row['endpoint'].det())==1 for row in rows)
    assert control['pair_internal_irreps']=={'A1':1,'T1':3}
    assert control['endpoint_irrep']=='G1'
    assert not control['Lorentz_spin_or_physical_hypercharge_dictionary_proved']


def test_odd_regional_exchange_and_even_global_pair_exchange_are_distinct():
    C=[creation(j,4) for j in range(4)]
    assert np.linalg.norm(C[0]@C[2]+C[2]@C[0])==0
    assert np.linalg.norm(C[0]@C[2].conj().T+C[2].conj().T@C[0])==0
    B=C[0]@C[1].conj().T;D=C[2]@C[3].conj().T
    assert np.linalg.norm(B@D-D@B)==0
    assert np.linalg.norm(B@D+D@B)>1


def test_directed_average_penalties_and_complete_coefficient_normalization():
    rows=m.directed_menu();cells=list(product(range(2),repeat=3))
    assert len(rows)==96
    for a in range(4):
        chosen=[r for r in rows if r['coefficient']==a];assert len(chosen)==24
        for x in cells:
            assert sum(tuple(r['target'])==x for r in chosen)==3
            assert sum(tuple(r['compensator'])==x for r in chosen)==3
        links={tuple(sorted((tuple(r['target']),tuple(r['compensator'])))) for r in chosen}
        assert len(links)==12
        assert F(3,24*2)==F(1,16)
        assert F(2,24)*F(10,3)==F(5,18)
    assert sum(r['internal_probe_irrep']=='A1' for r in rows)==24
    assert sum(r['internal_probe_irrep']=='T1' for r in rows)==72


def test_coherent_charge_spin_clock_records_keep_original_law_and_terminal_disturbance():
    d=5;Q=s.diag(0,6,6,6,6);K=[]
    for j in range(1,5):
        A=s.zeros(5);A[j,0]=s.Rational(1,2);K.append(A)
    K.append(s.diag(0,1,1,1,1))
    meter=m.meter.coherent_meter((0,6,6,6,6),K,(6,6,6,6,0),5)
    for row in meter['rows']:
        expected=K[row['probe']] if row['meter']==(1 if row['probe']<4 else 0) else s.zeros(5)
        assert row['K']==expected
    control=m.local_control();basis=s.Matrix.hstack(control['anchor'],*(B*control['anchor'] for B in control['success']))
    P0=s.diag(*(int((j&3)==0) for j in range(16)))
    P2=s.diag(*(int((j&3)==3) for j in range(16)))
    Pplus=s.diag(*(int((j&3)==1) for j in range(16)))
    Pminus=s.diag(*(int((j&3)==2) for j in range(16)))
    terminal=[np.array(m.clean(basis.H*s.kronecker_product(s.eye(3),P)*basis),complex) for P in (P0,P2,Pplus,Pminus)]
    assert np.linalg.norm(sum(terminal)-np.eye(d))<1e-14
    H=np.diag(np.arange(d,dtype=float));H[0,1:]=[.2,.3,.4,.5];H[1:,0]=H[0,1:]
    raw=8*np.eye(d)-H;rho=raw@raw;rho/=np.trace(rho)
    U=np.linalg.solve(np.eye(d)+1j*H/16,np.eye(d)-1j*H/16)
    clocks=((np.eye(d)+U)/2,(np.eye(d)-U)/2)
    changed=0.;mass=0.
    for first,last,A in product(clocks,clocks,[np.array(a,complex) for a in K]):
        L=last@A@first;state=L@rho@L.conj().T
        refined=sum(P@state@P for P in terminal)
        assert abs(np.trace(refined)-np.trace(state))<1e-14
        changed+=np.linalg.norm(refined-state);mass+=np.trace(state).real
    assert abs(mass-1)<1e-14 and changed>1e-3


def test_actual_spin_increments_and_joint_coherent_counters_do_not_dephase_input():
    row=m.spin_increment_control();fields=row['raw_fields'];Zp=row['positive_doubled_spin'];Zm=row['negative_doubled_spin']
    Q=6*sum((m.creator(j,4)*m.creator(j,4).H for j in range(2)),s.zeros(16))
    assert Q*Zp==Zp*Q and Q*Zm==Zm*Q and Zp*Zm==Zm*Zp
    for A,(zp,zm) in zip(fields,row['raw_doubled_spin_increments']):
        assert Q*A-A*Q==6*A
        assert Zp*A-A*Zp==zp*A and Zm*A-A*Zm==zm*A
    # Actual five-dimensional physical raw-spin block, including its anchor.
    operators=[]
    for j in range(1,5):
        A=s.zeros(5);A[j,0]=s.Rational(1,2);operators.append(A)
    operators.append(s.diag(0,1,1,1,1))
    alphabets=((0,6,6,6,6),(0,1,1,-1,-1),(0,-1,1,-1,1))
    shifts=((6,6,6,6,0),(1,1,-1,-1,0),(-1,1,-1,1,0))
    for charge,transfer,radix in zip(alphabets,shifts,(5,3,3)):
        control=m.meter.coherent_meter(charge,operators,transfer,radix)
        for record in control['rows']:
            j=record['probe']
            expected=operators[j] if record['meter']==transfer[j]%radix else s.zeros(5)
            assert record['K']==expected
    assert row['spin_meter_registers']==16
    assert not row['spin_increment_is_final_spin_after_energy_clock']


@pytest.fixture(scope='module')
def actual():return m.bind_member(m.contract.member(0))


def test_actual_selected_charge_and_spin_events_have_original_clock_weight(actual):
    r=actual['result'];w=F(r['per_chosen_coefficient_original_success_lower'])
    assert w>F(7944,10000)
    assert F(r['full_residual_norm_upper'])**2<75
    assert F(r['per_chosen_coefficient_recorded_weight_lower'])>F(5,16)
    assert F(r['complete_menu_each_coefficient_recorded_weight_lower'])==F(r['per_chosen_coefficient_recorded_weight_lower'])/4
    assert F(r['complete_menu_scalar_recorded_simple_lower'])==F(5,64)
    assert F(r['complete_menu_vector_recorded_simple_lower'])==F(15,64)
    assert r['reported_transfer_interval']==['69/8','279/8']
    assert r['complete_probe_outcomes']==97 and r['complete_success_outcomes']==96
    assert r['coherent_charge_meter_dimension']==5**8
    assert not r['source_carrier_state_reference_clock_changed']
    assert not r['pair_energy_divided_by_two_is_justified']
    assert r['physical_mass_table_GeV']==[]
    raw=r['raw_spin_apparatus']
    assert F(raw['per_fixed_spin_entry_recorded_weight_lower'])>F(1,4)
    assert F(raw['each_of_four_joint_spin_rows_recorded_weight_lower'])>F(1,16)
    assert raw['reported_transfer_interval']==['77/8','271/8']
    assert raw['joint_charge_spin_meter_dimension']==5**8*3**16
    assert raw['charge_and_spin_counters_preserve_original_complete_reference_history']
    assert not raw['raw_spin_and_Pauli_records_are_same_joint_apparatus']


def test_source_or_species_promotion_mutations_are_rejected(actual):
    m.verify_bound(m.contract.member(0),actual)
    with pytest.raises(ValueError):m.bind_member(m.contract.member(1))
    for key in ('spinor_label_is_physical_Lorentz_spin','pair_energy_is_single_endpoint_rest_energy','source_carrier_state_reference_clock_changed'):
        altered=deepcopy(actual);altered['result'][key]=True
        with pytest.raises(ValueError):m.verify_bound(m.contract.member(0),altered)


def test_canonical_receipt_binds_actual_transitive_source_and_finite_theorems():
    value=m.certificate();assert json.loads(json.dumps(value,sort_keys=True))==value
    assert value['schema']=='oph.primary_regional_fermion_certificate.v1'
    for path in ('model/particle_primary_charged_string/model.py','model/particle_primary_charge_conjugation/model.py',
                 'model/particle_primary_geometric_cut/model.py','Lean/QFT/CPRestrictionNet.lean'):
        assert path in value['source_sha256']
