"""Strict retained-evidence and source-custody tests for the exact certificate."""
from copy import deepcopy

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def replay():
    return packet.report()


def test_retained_receipt_reexecutes_complete_evidence(replay):
    packet.compare(strict_json.load(packet.OUTPUT), replay)
    assert replay['branch_algebra']['dimensions']['total'] == 264
    assert replay['same_candidate_dictionary']['new_absolute_physical_masses_predicted'] == 0
    assert replay['physical_particle_identification'] is False


@pytest.mark.parametrize('key,value', [('integrated_gate_accepted', True),
                                     ('physical_comparison_count_added', 1),
                                     ('bare_native_law_modified', True)])
def test_scope_promotions_reject(replay, key, value):
    changed = deepcopy(replay)
    changed[key] = value
    with pytest.raises(ValueError):
        packet.compare(changed, replay)


def test_exact_branch_count_odd_loop_and_gap_exclusion_are_mandatory(replay):
    edits = [('decomposition', 'integer_input_columns_checked', 15839),
             ('commutant', 'commutant_complex_dimension', 9),
             ('commutant', 'odd_chord', [0, 0, 0, 1])]
    for group, key, value in edits:
        changed = deepcopy(replay)
        changed['branch_algebra']['laws']['C'][group][key] = value
        with pytest.raises(ValueError):
            packet.compare(changed, replay)
    changed = deepcopy(replay)
    changed['same_candidate_dictionary']['common_scale_or_chemical_potential_splits_triplet'] = True
    with pytest.raises(ValueError):
        packet.compare(changed, replay)


@pytest.mark.parametrize('bad', [-1e-20, float('nan'), float('inf'), True, 0, '0'])
def test_negative_nonfinite_or_wrong_typed_residual_reject(bad):
    with pytest.raises(ValueError):
        packet.compare({'residual': bad}, {'residual': 0.0})


def test_changed_source_snapshot_rejects_warm_replay(monkeypatch):
    changed = packet.source_bindings()
    changed['model/charged_defect/branches.py'] = '0'*64
    monkeypatch.setattr(packet, 'source_bindings', lambda: changed)
    with pytest.raises(ValueError, match='fresh interpreter'):
        packet.report()


def test_transitive_source_and_independent_checks_are_bound():
    bindings = packet.source_bindings()
    for name in ('model/charged_defect/branches.py', 'model/charged_defect/dictionary.py',
                 'model/charged_defect/test_branches.py', 'model/charged_defect/test_dictionary.py',
                 'model/charged_defect/test_packet.py', 'model/charged_defect/contract.json',
                 'model/particle_transport/algebra.py', 'model/particle_transport/packet.py',
                 'model/exchange_sector/statistics.py', 'model/observer/native.py',
                 'model/full_fock_field/field.py', 'model/__init__.py',
                 'evidence/native_geometric_source_20260925/spec.json'):
        assert name in bindings
    assert 'model/charged_defect/results.json' not in bindings


def test_duplicate_or_extra_receipt_keys_reject(replay):
    with pytest.raises(ValueError):
        strict_json.loads('{"x": 1, "x": 2}')
    with pytest.raises(ValueError):
        packet.compare({**replay, 'extra': True}, replay)
