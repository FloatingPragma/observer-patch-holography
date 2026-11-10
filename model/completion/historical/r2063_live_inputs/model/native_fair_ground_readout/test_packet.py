"""Fail-closed source and complete readout custody controls."""
from copy import deepcopy
from fractions import Fraction as F
import json

import pytest

from . import packet, reference


@pytest.fixture(scope='module')
def fresh():
    return packet.produce()


def test_retained_receipt_matches_current_complete_owner(fresh):
    assert packet.check_payload(packet.strict_json.load(packet.OUTPUT), fresh)
    assert fresh['inherited_complete_energy_measurement'][
        'audited_full_66_word_receipt_sha256'] == packet.MEASUREMENT_RECEIPT_SHA256
    assert fresh['inherited_complete_energy_measurement']['full_66_word_execution_repeated_here'] is False


def test_actual_words_and_complete_score_use_one_native_chart(fresh):
    row = fresh['executed_fair_ground_preparation']
    assert row['ground_word_declared_resources']['coherent_work_qubits'] == 9
    assert row['ground_word_declared_resources']['native_pulses'] == row['actual_compiled_pulses']
    assert row['local_phase_word_sha256'] != row['matched_control_word_sha256']
    assert row['budget']['formal_raw_words'] == 25**11
    assert row['budget']['all_raw_records_and_first_processor_padding_failures_retained'] is True
    assert row['energy_readout_join']['complete_unnormalized_score_contrast_strict_lower'] == '3/50'
    assert row['energy_readout_join']['stable_particle_pole_or_mass_identified'] is False


@pytest.mark.parametrize('mutation', ('word', 'erased_outcomes', 'error', 'score',
    'source', 'thermal_as_ground', 'particle', 'typed_alias'))
def test_retained_scientific_mutations_fail(fresh, mutation):
    changed = deepcopy(fresh)
    row = changed['executed_fair_ground_preparation']
    if mutation == 'word': row['ground_word_sha256'] = '0'*64
    elif mutation == 'erased_outcomes': row['budget']['formal_raw_words'] = 66*2**11
    elif mutation == 'error': row['budget']['complete_halfdiamond_error_strict_upper'] = '0'
    elif mutation == 'score': row['energy_readout_join']['complete_unnormalized_score_contrast_strict_lower'] = '1'
    elif mutation == 'source': changed['source_bindings'].pop(next(iter(changed['source_bindings'])))
    elif mutation == 'thermal_as_ground': changed['scope']['ground_equals_faithful_transfer_reference'] = True
    elif mutation == 'particle': changed['scope']['relativistic_3_plus_1_continuum_or_particles_identified'] = True
    else: changed['scope']['finite_native_ground_preparation'] = 1
    with pytest.raises(ArithmeticError, match='fresh owning replay'):
        packet.check_payload(changed, fresh)


def test_prior_measurement_word_mutation_fails_exact_pin(tmp_path, monkeypatch):
    retained = packet.strict_json.load(packet.measurement_owner.OUTPUT)
    retained['all_native_words'][0]['word_sha256'] = '0'*64
    altered = tmp_path/'measurement.json'
    altered.write_text(packet.measurement_owner.canonical(retained))
    monkeypatch.setattr(packet.measurement_owner, 'OUTPUT', altered)
    with pytest.raises(ArithmeticError, match='receipt SHA differs'):
        packet.owning_measurement()


def test_changed_prior_live_source_is_rejected(monkeypatch):
    pins = packet.measurement_owner.source_bindings()
    pins[next(iter(pins))] = '0'*64
    monkeypatch.setattr(packet.measurement_owner, 'source_bindings', lambda: pins)
    with pytest.raises(ArithmeticError, match='live sources or consumed bounds'):
        packet.owning_measurement()


def test_mid_replay_source_change_fails(fresh, monkeypatch):
    before = fresh['source_bindings']
    after = dict(before, unexpected='0'*64)
    calls = iter((before, after))
    monkeypatch.setattr(packet, 'source_bindings', lambda: next(calls))
    monkeypatch.setattr(packet, 'owning_measurement', lambda: fresh['inherited_complete_energy_measurement'])
    monkeypatch.setattr(packet.preparation, 'run', lambda: fresh['executed_fair_ground_preparation'])
    with pytest.raises(ArithmeticError, match='changed during'):
        packet.produce()


def test_source_inventory_is_local_and_has_no_research_note(fresh):
    paths = set(fresh['source_bindings'])
    assert set(packet.OWNED_FILES) <= paths
    assert 'model/native_energy_measurement/results.json' in paths
    assert 'model/native_fair_ground_readout/results.json' not in paths
    assert not any(path.startswith('docs/') for path in paths)
    assert packet.ROOT == packet.measurement_owner.ROOT


def test_bounded_local_arithmetic_does_not_allocate_large_regions():
    with pytest.raises(ValueError, match='64'):
        reference.dicke_weights(10**122, 10**121, 10**120)
    with pytest.raises(ValueError, match='64'):
        reference.limiting_weights(10**122, F(1, 2))
    with pytest.raises(ValueError):
        reference.history_budget((1,)*65, (F(1),)*65, 30)


def test_duplicate_and_nonfinite_receipts_fail(fresh, tmp_path, monkeypatch):
    monkeypatch.setattr(packet, 'produce', lambda: fresh)
    for raw in ('{"schema":1,"schema":2}', '{"value":NaN}'):
        path = tmp_path/'bad.json'; path.write_text(raw)
        with pytest.raises((ValueError, TypeError)):
            packet.verify(path)
