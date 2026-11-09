"""Independent charge-counter, CAR/Gauss and complete-history controls."""
from copy import deepcopy
from fractions import Fraction as F

import pytest
import sympy as s

from . import model as m


def destroy(count,index):
    out=s.zeros(1<<count)
    for mask in range(1<<count):
        if mask&(1<<index):
            out[mask^(1<<index),mask]=(-1)**((mask&((1<<index)-1)).bit_count())
    return out


def test_actual_native_range_and_minimal_cyclic_transfer_menu():
    data=m.native_charge_menu()
    assert data['local_charge_range']==[-72,72]
    assert data['matter_modes_per_cell']==60 and data['difference_meter_radix']==5
    assert [m.transfer_record(k) for k in range(5)]==[0,6,None,None,-6]
    assert len({q%5 for q in (-6,0,6)})==3
    assert all(len({q%radix for q in (-6,0,6)})<3 for radix in range(1,5))
    assert 1%5==6%5  # The admitted counter is not a generic charge-difference PVM.
    with pytest.raises(ValueError):m.transfer_record(-1)
    with pytest.raises(ValueError):m.transfer_record(5)


def test_every_actual_bond_component_has_opposite_endpoint_record():
    rows=m.source_shift_table();success=rows[:-1]
    assert len(rows)==49 and rows[-1]['charge_shift']==[0]*8
    assert len({tuple(r['charge_shift']) for r in success})==24
    for row in success:
        assert sorted(row['charge_shift'])==[-6,0,0,0,0,0,0,6]
        assert sum(abs(a-b) for a,b in zip(row['target'],row['origin']))==1
    assert sum(r['matched_beta_component']==0 for r in success)==24


def test_three_state_control_has_a_literal_physical_CAR_link_embedding():
    c=[destroy(4,j) for j in range(4)];n=[x.H*x for x in c]
    E=s.diag(-1,0,1);p0=s.diag(0,1,0)
    plus=s.zeros(3);plus[2,1]=1
    minus=s.zeros(3);minus[0,1]=1
    # The minus sign of the positive-charge vector is its literal CAR sign.
    V=s.zeros(48,3);V[26,0]=1;V[35,1]=-1;V[12,2]=1
    qx=s.kronecker_product(s.eye(3),n[0]+n[1]-s.eye(16))
    qy=s.kronecker_product(s.eye(3),n[2]+n[3]-s.eye(16))
    flux=s.kronecker_product(E,s.eye(16))
    bp=s.kronecker_product(plus,c[0].H*c[3])
    bm=s.kronecker_product(minus,c[2].H*c[1])
    assert V.H*V==s.eye(3)
    assert (qx-flux)*V==(qy+flux)*V==s.zeros(48,3)
    assert V.H*qx*V==s.diag(0,1,-1)
    assert 6*V.H*qx*V==m.quantum_control()['Q']
    expected_plus=s.zeros(3);expected_plus[1,0]=1
    expected_minus=s.zeros(3);expected_minus[2,0]=1
    assert V.H*bp*V==expected_plus and V.H*bm*V==expected_minus
    assert V.H*s.kronecker_product(p0,s.eye(16))*V==s.diag(1,0,0)


def test_coherent_meter_factors_on_every_source_and_meter_basis_vector():
    q=m.quantum_control();W=q['meter']['coupling'];d=3;M=5
    for b,(K,shift) in enumerate(zip(q['operators'],(6,-6,0))):
        for source in range(d):
            for meter in range(M):
                v=s.eye(d*M)[:,source*M+meter]
                actual=m.clean(W*s.kronecker_product(K,s.eye(M))*W.H*v)
                expected=s.kronecker_product(K*s.eye(d)[:,source],s.eye(M)[:,(meter+shift)%M])
                assert actual==expected
        rows=[r for r in q['meter']['rows'] if r['probe']==b]
        assert len(rows)==M
        assert all(r['K']==(K if r['meter']==shift%M else s.zeros(d)) for r in rows)


def test_common_entangled_reference_and_energy_histories_are_preserved():
    q=m.quantum_control();root=q['root'];rho=q['rho'];operators=q['operators'];energies=q['energy_projectors']
    vec=s.Matrix([root[i,j] for i in range(3) for j in range(3)])
    assert s.simplify((vec.H*vec)[0])==1
    total=0;record=0
    for ei in energies:
        for b,K in enumerate(operators):
            for ef in energies:
                L=ef*K*ei;v=s.kronecker_product(L,s.eye(3))*vec
                probability=s.simplify((v.H*v)[0]);total+=probability
                assert probability==q['original_energy_probabilities'][record];record+=1
                for row in [r for r in q['meter']['rows'] if r['probe']==b]:
                    Lv=ef*row['K']*ei
                    actual=s.kronecker_product(Lv,s.eye(3))*vec
                    assert actual==(v if row['meter']==(6,-6,0)[b]%5 else s.zeros(9,1))
    assert total==1


def test_terminal_charge_is_complete_but_not_the_transfer_event_label():
    q=m.quantum_control()
    assert q['H']*q['Q']!=q['Q']*q['H']
    assert q['H']*q['rho']==q['rho']*q['H']
    assert q['post_energy_charge_not_equal_probe_endpoint_probability']>0
    # All energy probabilities survive the later charge PVM, even though its
    # source/reference postmeasurement states have a further disturbance.
    for ei in q['energy_projectors']:
        for K in q['operators']:
            for ef in q['energy_projectors']:
                L=ef*K*ei;state=L*q['rho']*L.H
                refinements=[s.diag(*(int(i==j) for i in range(3))) for j in range(3)]
                assert s.simplify(sum(s.trace(P*state*P) for P in refinements)-s.trace(state))==0


def test_premature_charge_measurement_changes_original_energy_probabilities():
    q=m.quantum_control()
    assert q['premature_charge_readout_record_distance']>s.Rational(1,10)
    with pytest.raises(ValueError,match='charge-shift'):
        m.coherent_meter((0,6,-6),q['operators'],(6,6,0),5)
    with pytest.raises(ValueError,match='complete'):
        m.coherent_meter((0,6,-6),q['operators'][:2],(6,-6),5)


def test_finite_coupling_errors_have_no_outcome_or_success_multiplier():
    assert m.encoder_error(8,F(1,4096))==F(1,256)
    assert m.encoder_error(8,1)==1
    assert m.encoder_error(8,0)==0
    with pytest.raises(ValueError):m.encoder_error(8,F(-1,100))


@pytest.fixture(scope='module')
def selected():
    member=m.contract.member(0)
    return member,m.bind_member(member)


def test_actual_original_charged_energy_bound_is_retained_on_one_history(selected):
    member,bound=selected;r=bound['result']
    assert r['difference_meter_registers']==8 and r['difference_meter_dimension']==5**8
    assert r['absolute_charge_instrument_outcomes']==145**8
    assert r['one_coupling_controlled_CAR_number_shifts']==480
    assert r['two_couplings_controlled_CAR_number_shifts']==960
    assert F(r['complete_encoder_half_trace_error_upper'])==0
    assert F(r['reported_joint_opposite_charge_transfer_original_weight_lower'])>F(1,4)
    assert r['reported_transfer_interval']==['77/8','271/8']
    assert not r['meter_read_before_cancellation'] and not r['initial_absolute_charge_record_retained']
    assert r['complete_original_source_reference_history_preserved_exactly_for_exact_gates']
    assert not r['source_clock_energy_zero_changed'] and not r['source_preparation_changed']
    assert not r['terminal_charge_inferred_equal_to_transfer_label']
    assert not r['large_source_or_meter_histogram_executed']
    assert r['positive_event_requires_probe_success'] and not r['positive_event_includes_failure_energy_records']
    assert r['meter_records_distinguish_proved_probe_transfer_alphabet']
    assert not r['arbitrary_local_charge_differences_resolved']
    assert r['physical_mass_table_GeV']==[]
    m.verify_bound(member,bound)


def test_changed_source_result_and_requested_precision_are_rejected(selected):
    member,bound=selected
    with pytest.raises(ValueError):m.bind_member(m.contract.member(1))
    wrong=deepcopy(bound);wrong['result']['reported_joint_opposite_charge_transfer_original_weight_lower']='1'
    wrong['result_sha256']=m.contract.digest(wrong['result'])
    with pytest.raises(ValueError):m.verify_bound(member,wrong)
    with pytest.raises(ValueError):m.verify_bound(member,bound,per_cell_encoder_error='1/4096')
    noisy=m.bind_member(member,per_cell_encoder_error='1/4096')['result']
    assert F(noisy['reported_joint_opposite_charge_transfer_original_weight_lower'])==F(bound['result']['reported_joint_opposite_charge_transfer_original_weight_lower'])-F(1,256)
