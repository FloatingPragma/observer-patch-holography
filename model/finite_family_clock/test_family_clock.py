"""Meaningful source, protocol and mutation checks for the finite package."""
import copy
import json
from pathlib import Path

import pytest
import sympy as s

from . import primitives as p
from .certificate import certificate,input_hashes,replay_differences,validate_pins
from .family import calculate as family
from .memory import calculate as memory


def test_native_event_and_actual_repeated_tick():
    value=family()
    assert value['native_event']['actual_Observer_repair_columns']==360
    assert value['exact_instrument']['two_tick_compression_for_all_phases']=='sqrt(3)/9 * I3'
    assert value['multi_tick'][1]['distance_from_M_power']>0.1
    assert value['actual_FR_Floquet']['FR_report_arccos_readout_is_eigenphase_of_this_U'] is False


def test_orientation_mutation_changes_the_instrument():
    prepared=p.state(s.sqrt(s.Rational(1,3)),0,s.sqrt(s.Rational(2,3)))
    bare=p.repair()*p.address_reflection()*prepared
    wrong=p.compression(p.analyzer(),bare)
    right=s.sqrt(2)/3*p.circulant()
    assert (wrong-right).applyfunc(s.simplify)!=s.zeros(3)


def test_complete_effect_and_collision_bound():
    value=family()
    effect=value['complete_binary_instrument']
    for a,b in zip(effect['probability'],effect['no_operation_probability'],strict=True):
        assert abs(a+b-1)<1e-14
    for row in value['collision_clock']['rows']:
        assert row['max_coherence_multiplier_error']<=row['diamond_error_upper_bound']
    assert value['collision_clock']['discarded_outcomes']==[]


def test_memory_complete_spectrum_and_same_H():
    value=memory()
    assert value['same_native_one_family']['complete_native_columns_compared']==320
    assert value['explicit_preparation']['exact_Householder_maps_blank_to_ground']
    for row in value['families']:
        assert row['composite_independently_compiled_from_native_CAR_and_scalar']
        assert row['complete_dark_spectral_weight_exact']=='0'
        assert row['diagnostic']['dark_multiplicity']==9
        assert len(row['diagnostic']['lines'])==2
        assert row['diagnostic']['complete_elementary_record_probability']==1
        assert row['diagnostic']['same_H_amplitude_identity_max_error']<1e-12
    assert value['physical_rest_mass_proved'] is False


def test_sign_obstruction_and_nonidentical_carriers():
    value=memory()
    lo=float(value['baseline_squared_circulant_ratios'][0])
    assert all(float(row['mu/e'])>lo for row in value['finite_memory_mass_rows'])
    dims=family()['register_and_Fock_are_distinct']
    assert dims['Fock_K_tensor_common_register_dimension']==192
    assert dims['Fock_K_tensor_three_family_one_particle_space_dimension']==262144


def test_source_closure_and_immutable_parents():
    pins=input_hashes()
    validate_pins(pins)
    for name in ('model/coherent_seam_repair/program.py',
                 'model/coherent_defect_wave/chiral_joint/model.py',
                 'model/coherent_defect_wave/chiral_higgs/model.py',
                 'model/finite_family_clock/primitives.py'):
        assert name in pins
    wrong=dict(pins)
    name='model/coherent_seam_repair/program.py'
    wrong[name]='0'*64
    with pytest.raises(ValueError,match='source custody mismatch'):
        validate_pins(wrong)


def test_receipt_replays_and_rejects_mutations():
    expected=json.loads(Path(__file__).with_name('results.json').read_text())
    actual=certificate()
    assert replay_differences(expected,actual)==[]
    for mutate in ('scope','exact','spectrum'):
        wrong=copy.deepcopy(actual)
        if mutate=='scope':
            wrong['scope']['primary_common_source_transfer']=True
        elif mutate=='exact':
            wrong['family_instrument_and_collision_clock']['exact_instrument']['trace_M_dagger_M']='1/9'
        else:
            wrong['same_H_quantum_memory']['families'][0]['diagnostic']['lines'][0]['energy_transfer']+=0.01
        assert replay_differences(wrong,actual)
