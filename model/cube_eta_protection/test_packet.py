"""Exact claims, full records and owned cold source custody fail closed."""
from copy import deepcopy

import pytest
from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def replay(): return packet.report()


def test_retained_full_cube_packet_cold_replays(replay):
    retained = strict_json.load(packet.OUTPUT)
    packet.compare(retained, replay); packet.validate_accounting(retained)


@pytest.mark.parametrize('key', ['integrated_gate_accepted', 'physical_particle_identification',
                                'primitive_or_repair_law_changed', 'native_time_dependent_detector_executed',
                                'physical_hardware_execution_or_accuracy_proved'])
def test_scope_promotion_rejects(replay, key):
    changed = deepcopy(replay); changed[key] = True
    with pytest.raises(ValueError): packet.compare(changed, replay)


def test_exact_polynomial_dimension_or_source_bound_cannot_drift(replay):
    for section, key, wrong in [
        ('exact_full_field_eta_protection', 'protected_eta_image_dimension', 1821),
        ('exact_interacting_spectral_census', 'eta3_distinct_energies', 12),
        ('exact_native_source_error_certificate', 'complete_source_state_error_bound', '1/100')]:
        changed = deepcopy(replay); changed[section][key] = wrong
        with pytest.raises(ValueError): packet.compare(changed, replay)


def test_records_padding_and_positive_spectral_lines_cannot_be_dropped(replay):
    changed = deepcopy(replay); changed['native_reference_and_local_phase']['all130_following_repair_records'].pop(0)
    with pytest.raises(ValueError): packet.validate_accounting(changed)
    changed = deepcopy(replay); row = changed['native_reference_and_local_phase']['source_records'][0]
    row['all_first_export_failure_weights'][row['pair_code_event_indices'][0]][3] = 1e-30
    with pytest.raises(ValueError): packet.validate_accounting(changed)
    changed = deepcopy(replay); changed['exact_fixed_number_local_response']['eight_interacting_energy_lines'].pop()
    with pytest.raises(ValueError): packet.validate_accounting(changed)


@pytest.mark.parametrize('bad', [-1e-30, float('nan'), float('inf'), True, 0, '0'])
def test_invalid_numerical_evidence_rejects(bad):
    with pytest.raises(ValueError): packet.compare({'probability': bad}, {'probability': 0.})


def test_source_change_rejected_before_owned_worker(monkeypatch):
    changed = packet.source_bindings(); changed['model/cube_eta_protection/cube.py'] = '0'*64
    monkeypatch.setattr(packet, 'source_bindings', lambda: changed)
    with pytest.raises(ValueError, match='fresh interpreter'): packet.report()


def test_source_manifest_contains_consumed_laws_exact_checks_and_paper():
    paths = packet.source_bindings()
    for name in ('cube.py', 'repair.py', 'spectrum.py', 'source.py', 'packet.py', 'test_cube.py',
                 'test_source.py', 'test_packet.py', 'contract.json', 'README.md'):
        assert 'model/cube_eta_protection/'+name in paths
    for name in ('model/observer/source_scalar/compiler.py', 'model/certified_reference/source.py',
                 'model/field_symmetry/symmetry.py', 'model/particle_transport/algebra.py',
                 'paper/tex_fragments/CUBE_ETA_PROTECTION.tex'):
        assert name in paths
    assert 'model/cube_eta_protection/results.json' not in paths
