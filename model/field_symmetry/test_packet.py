"""Owned receipt replay refuses missing, invented and malformed evidence."""
from copy import deepcopy
import json
import math
from pathlib import Path

import pytest

from model.field_symmetry import packet as p


@pytest.fixture(scope='module')
def fresh():
    return p.report()


def test_actual_producer_round_trip(fresh, tmp_path, monkeypatch):
    retained = json.loads(json.dumps(fresh, allow_nan=False))
    p.compare(retained, fresh)
    path = tmp_path/'receipt.json'
    path.write_text(json.dumps(retained), encoding='utf-8')
    monkeypatch.setattr(p, 'report', lambda: fresh)
    assert p.verify(path) == fresh
    assert all(row['complete_Fock_dimension'] == 4**row['sites']
               for row in fresh['symmetry']['full_carrier_checks'])
    assert fresh['symmetry']['full_carrier_checks'][-1]['complete_Fock_dimension'] == 65536


@pytest.mark.parametrize('mutation', ['empty', 'missing_mode', 'wrong_population', 'wrong_count',
    'bool_count', 'nonfinite', 'negative_probability', 'wrong_charge', 'false_physics',
    'stale_hash', 'failed_primitive', 'wrong_budget', 'extra_key', 'wrong_symmetry'])
def test_mutated_receipts_fail(fresh, mutation):
    bad = deepcopy(fresh)
    row = bad['complete_source_experiments'][0]
    if mutation == 'empty': bad = {}
    elif mutation == 'missing_mode': row['terminal_probabilities'].pop()
    elif mutation == 'wrong_population': row['terminal_probabilities'][0] = .9
    elif mutation == 'wrong_count': row['resources']['native_pulses'] += 1
    elif mutation == 'bool_count': row['steps'] = True
    elif mutation == 'nonfinite': row['complete_output_error'] = float('nan')
    elif mutation == 'negative_probability': row['terminal_probabilities'][255] = -1e-300
    elif mutation == 'wrong_charge': bad['symmetry']['eta_pair_charge'] = 1
    elif mutation == 'false_physics': bad['physical_species_identified'] = True
    elif mutation == 'stale_hash': bad['source_bindings'][next(iter(bad['source_bindings']))] = '0'*64
    elif mutation == 'failed_primitive': row['raw_primitive_residual'] = 1.1e-11
    elif mutation == 'wrong_budget': row['state_specific_ideal_split_error'] = '1/4'
    elif mutation == 'extra_key': row['invented_success'] = True
    elif mutation == 'wrong_symmetry': bad['symmetry']['full_carrier_irreps']['4'][0]['multiplicity'] += 1
    with pytest.raises(ValueError): p.compare(bad, fresh)


def test_duplicate_json_refused_before_replay(tmp_path, monkeypatch):
    path = tmp_path/'bad.json'; path.write_text('{"schema":1,"schema":2}', encoding='utf-8')
    monkeypatch.setattr(p, 'report', lambda: pytest.fail('ambiguous evidence must not execute'))
    with pytest.raises(ValueError): p.verify(path)


def test_warm_source_change_refused_before_production(monkeypatch):
    bad = dict(p._IMPORTED_BINDINGS); bad[next(iter(bad))] = '0'*64
    monkeypatch.setattr(p, 'source_bindings', lambda: bad)
    monkeypatch.setattr(p.symmetry, 'certificate', lambda: pytest.fail('stale code must not execute'))
    with pytest.raises(ValueError, match='since import'): p.report()


def test_public_verifier_owns_replay(tmp_path, monkeypatch):
    path = tmp_path/'bad.json'; path.write_text('{}', encoding='utf-8')
    called = []
    def owned():
        called.append(True)
        return {'schema': 'real-science'}
    monkeypatch.setattr(p, 'report', owned)
    with pytest.raises(ValueError): p.verify(path)
    assert called == [True]


def test_all_owning_sources_and_primitive_inputs_are_bound(fresh):
    bindings = fresh['source_bindings']
    for name in ('model/field_symmetry/readout.py', 'model/field_symmetry/symmetry.py',
                 'model/field_symmetry/test_topology.py', 'model/field_symmetry/test_packet.py',
                 'model/observer/source_scalar/compiler.py', 'model/field_thermal/thermal.py',
                 'model/field_stress/energy.py', 'tools/strict_json.py',
                 'code/source_selection_model/response.json', 'code/m1_source_realization/receipt.json'):
        assert name in bindings
    assert all('\\' not in name for name in bindings)


def test_retained_current_receipt():
    p.verify()


def test_preimport_warm_primitive_cache_cannot_authenticate_changed_data(monkeypatch):
    # A source cache may predate packet import. A newly computed manifest
    # alone cannot authorize numerical products made from the old data.
    p.readout.execute()
    control = p.ROOT/p.readout.compiler.CONTROL_RECEIPT
    read_bytes, read_text = Path.read_bytes, Path.read_text
    def bytes_now(path):
        return b'{}' if path == control else read_bytes(path)
    def text_now(path, *args, **kwargs):
        return '{}' if path == control else read_text(path, *args, **kwargs)
    with monkeypatch.context() as m:
        m.setattr(Path, 'read_bytes', bytes_now)
        m.setattr(Path, 'read_text', text_now)
        # Same state as loading the packet after the source data changed.
        m.setattr(p, '_IMPORTED_BINDINGS', p.source_bindings())
        m.setattr(p.symmetry, 'certificate', lambda: {})
        with pytest.raises((ValueError, KeyError)):
            p.report()


@pytest.mark.parametrize('case', ['centered', 'free', 'H_forbidden', 'failure', 'marginal',
                                 'near_zero_marginal'])
def test_public_verifier_rejects_coordinated_probability_budget_bypass(fresh, tmp_path, monkeypatch, case):
    bad = deepcopy(fresh)
    index = {'centered': 3, 'free': 4, 'near_zero_marginal': 3}.get(case, 0)
    row = bad['complete_source_experiments'][index]
    delta = 1e-11  # Below replay tolerance, above the squared amplitude budget.
    if case == 'near_zero_marginal':
        row['read_one_probabilities'][0] = 1e-14
    elif case == 'marginal':
        # This is within the state-distance probability allowance, but is not
        # the marginal of the retained complete terminal distribution.
        row['read_one_probabilities'][0] += delta
    elif case == 'failure':
        row['terminal_probabilities'][0] -= delta
        row['first_export_failure_probability'] += delta
    else:
        row['terminal_probabilities'][0] -= delta
        row['terminal_probabilities'][1] += delta  # A singly occupied site.
        row['read_one_probabilities'][0] += delta
    path = tmp_path/'coordinated.json'
    path.write_text(json.dumps(bad, allow_nan=False), encoding='utf-8')
    monkeypatch.setattr(p, 'report', lambda: fresh)
    with pytest.raises(ValueError):
        p.verify(path)


def test_retained_probabilities_allow_small_nonzero_population(fresh, tmp_path, monkeypatch):
    valid = deepcopy(fresh)
    row = valid['complete_source_experiments'][3]
    delta = 1e-18  # Nonzero, but its amplitude is below the fixed source budget.
    row['terminal_probabilities'][1] += delta
    row['terminal_probabilities'][0] -= delta
    row['read_one_probabilities'][0] = math.fsum(
        value for address, value in enumerate(row['terminal_probabilities']) if address & 1)
    path = tmp_path/'small-positive.json'
    path.write_text(json.dumps(valid, allow_nan=False), encoding='utf-8')
    monkeypatch.setattr(p, 'report', lambda: fresh)
    assert p.verify(path) == fresh


def test_retained_probabilities_accept_actual_complete_source_leakage(fresh, tmp_path, monkeypatch):
    import numpy as np
    original = p.readout.compiler.native_pulse
    angle = 3e-12
    rotation = np.eye(6, dtype=complex)
    rotation[0, 0] = rotation[5, 5] = math.cos(angle)
    rotation[0, 5], rotation[5, 0] = -math.sin(angle), math.sin(angle)
    monkeypatch.setattr(p.readout.compiler, 'native_pulse',
                        lambda kind, value: rotation@original(kind, value))
    current = deepcopy(fresh)
    current['complete_source_experiments'][0] = p.readout.execute()['report']
    assert current['complete_source_experiments'][0]['first_export_failure_probability'] > 1e-22
    path = tmp_path/'actual-positive-leakage.json'
    path.write_text(json.dumps(current, allow_nan=False), encoding='utf-8')
    monkeypatch.setattr(p, 'report', lambda: current)
    assert p.verify(path) == current
