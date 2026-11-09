"""Independent algebra, original-error budgets and exact replay controls."""
import copy
from fractions import Fraction as F
import json
import subprocess
import sys

import pytest

from model.spin_twist_closure import field_gap
from . import algebra, bounds, packet


def test_exact_source_selection_double_commutators():
    data=algebra.local_certificate()
    assert all(row['all_symbolic_matrix_entries_zero'] for row in data['all_three_axes'])
    assert [row['selection_commutes_exactly'] for row in data['all_three_axes']]==[True,False,False]
    assert data['selection_overlap_coefficient']==6


def test_actual_finite_source_gap_and_centered_variance():
    assert all(row['centered_variance_gap_control_passed'] for row in algebra.finite_ground_controls())


@pytest.mark.parametrize('n',[2,3,8,32])
def test_same_family_concentration_window(n):
    row=bounds.family(n,F(7,3))
    assert row['field']/(row['coupling']*row['p']**2)==F(1,n)
    assert row['coupling']*row['p']**3/row['field']==F(1,n)
    order=bounds.order_concentration(n,1,1,F(7,3))
    assert order['variance_upper']==F(48*27,n)
    bond=bounds.bond_concentration(n,1,1,F(7,3))
    assert bond['variance_upper']==F(792,n)+F(24,n**6)
    assert bond['variance_upper']<=bond['simple_variance_upper']


def test_selection_cost_is_small_on_both_fields_and_all_energy_histories():
    q=bounds.field_energy(8,'Q',3,2)
    p=bounds.field_energy(8,'P',F(1,4),2)
    assert q['selection_scaled_cost']==F(1,2*8**5)
    assert p['selection_scaled_cost']==F(1,16)
    assert bounds.history_energy([1,4,9])==36


def test_scalar_coefficient_history_transport_decreases():
    values=[bounds.bond_transport(n,1,1,1,1)['distance_squared_upper'] for n in (32,64,128)]
    assert values[0]>values[1]>values[2]>0


@pytest.mark.parametrize('sites,edges',[(1,()),(2,((0,1),)),(4,((0,1),(1,2),(2,3),(0,3)))])
def test_count_budget_equals_literal_native_ramp_budget(sites,edges):
    original=field_gap.resource_plan(sites,edges,F(2),F(1,3),F(1,20))
    counted=bounds.ground_preparation_counts(sites,len(edges),len(original['physical_odd_sites']),
                                            F(2),F(1,3),F(1,20))
    assert counted['source_clock_duration']==original['source_clock_duration']
    assert counted['slices']==original['piecewise_constant_slices']
    assert counted['native_pulses_upper']==original['native_pulses_upper']
    assert counted['total_half_trace_upper']<=F(1,20)
    assert counted['all_native_failure_records_retained']
    assert not counted['ground_postselection']


def test_actual_cube_bipartition_and_counts():
    for radius in (1,2,3):
        coords=[(x,y,z) for x in range(-radius,radius+1)
                for y in range(-radius,radius+1) for z in range(-radius,radius+1)]
        row=bounds.cube_counts(radius)
        assert row['sites']==len(coords)
        assert row['odd_sites']==sum(sum(c)%2 for c in coords)
        assert row['edges']==sum(sum(c[axis]<radius for axis in range(3)) for c in coords)


def test_computable_local_preparation_uses_same_positive_field():
    row=bounds.local_preparation(2,1,1,F(1,10))
    assert row['prepared_marginal_half_trace_upper']<=F(1,10)
    assert row['finite_source']['radius']%2==0
    assert row['preparation']['gap_lower']==row['field']/8
    assert not row['ordered_zero_field_state_accuracy_claimed']
    assert not row['finite_source_matrix_allocated']


@pytest.mark.parametrize('operation',[
    lambda:bounds.family(1),lambda:bounds.family(True),lambda:bounds.family(2,0),
    lambda:bounds.family(2,1.0),lambda:bounds.order_concentration(2,F(1,2),1),
    lambda:bounds.field_energy(2,'bad',1,1),
    lambda:bounds.local_preparation(2,1,1,1),
    lambda:bounds.ground_preparation_counts(1,0,2,1,1,F(1,10))])
def test_invalid_or_unselected_inputs_rejected(operation):
    with pytest.raises(ValueError):operation()


def test_receipt_scope_and_mutation_rejection():
    data=packet.build()
    assert data['scope']['direct_single_pair_regular_positive_energy_history_subsequence']
    uniform=data['uniform_mode_original_weight_and_concentration']
    assert uniform['uniform_mode_centre_and_weighted_variance_vanish_on_joint_clock']
    assert not uniform['uniform_original_probability_floor_proved']
    spectral=data['complete_finite_uniform_spectral_readout']
    assert [row['n']for row in spectral['finite_cases']]==[2,3,4]
    assert not spectral['fixed_box_table_substituted_for_growing_volume_field']
    assert spectral['physical_mass_prediction_count']==0
    speed=data['uncut_speed_concentration_from_source_susceptibility']
    assert packet.canonical(speed)==packet.canonical(packet.serialize(packet.speed_concentration.certificate()))
    assert speed['finite_source_control']['complete_record_count']==6
    assert speed['exact_source_speed_tube']['original_inside_probability_lower']=='4/5'
    assert data['scope']['uncut_original_weight_speed_tail_from_one_source_susceptibility']
    assert not data['scope']['actual_source_speed_concentration_defect_vanishes']
    for name in (packet.uniform_spectral.WITNESS,packet.uniform_spectral.zero_field.WITNESS):
        assert data['source_bindings'][name]==packet.hashlib.sha256((packet.ROOT/name).read_bytes()).hexdigest()
    changed=copy.deepcopy(data)
    changed['scope']['equality_to_ordered_zero_field_state_clock_limit_proved']=True
    with pytest.raises(ValueError):packet.validate(changed)
    changed=copy.deepcopy(data)
    changed['uncut_speed_concentration_from_source_susceptibility']['actual_ordered_source_concentration_proved']=0
    with pytest.raises(ValueError):packet.validate(changed)


def test_runtime_imports_are_all_source_bound_in_fresh_process():
    code='''import json,sys
from pathlib import Path
from model.spin_joint_field_limit import packet
value=packet.build(); root=packet.ROOT.resolve(); names=value['source_bindings']; missing=[]
for module in tuple(sys.modules.values()):
    path=getattr(module,'__file__',None)
    if path:
        path=Path(path).resolve()
        if path.suffix=='.py' and path.is_relative_to(root):
            name=path.relative_to(root).as_posix()
            if name not in names:missing.append(name)
print(json.dumps(sorted(set(missing))))'''
    result=subprocess.run([sys.executable,'-W','error','-c',code],cwd=packet.ROOT,
                           check=True,capture_output=True,text=True)
    assert json.loads(result.stdout)==[]
