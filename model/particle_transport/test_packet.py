"""Hostile receipt and custody checks, separate from the operator witnesses."""
from copy import deepcopy
import json
from pathlib import Path

import pytest

from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def replay():
    return packet.report()


def _leaves(value, path=()):
    if type(value) is dict:
        for key, child in value.items():
            yield from _leaves(child, (*path, key))
    elif type(value) is list:
        for index, child in enumerate(value):
            yield from _leaves(child, (*path, index))
    else:
        yield path, value


def test_committed_receipt_reexecutes_every_component(replay):
    packet.compare(strict_json.load(packet.OUTPUT), replay)
    assert replay['same_candidate_algebra']
    assert replay['source_transport']
    assert replay['number_conserving_detector']
    assert replay['physical_particle_identification'] is False
    assert replay['physical_comparison_count_added'] == 0
    assert replay['failed_dispersion_benchmark_repaired'] is False


def test_every_retained_leaf_matters(replay):
    for path, value in _leaves(replay):
        changed = deepcopy(replay)
        parent = changed
        for key in path[:-1]:
            parent = parent[key]
        if type(value) is bool:
            mutant = not value
        elif type(value) is int:
            mutant = value+1
        elif type(value) is float:
            mutant = value+max(1e-5, abs(value)*1e-4)
        elif type(value) is str:
            mutant = value+' altered'
        else:
            mutant = 'unexpected replacement'
        parent[path[-1]] = mutant
        with pytest.raises(ValueError):
            packet.compare(changed, replay)


@pytest.mark.parametrize('bad', [float('nan'), float('inf'), -float('inf'), True, 0, '0'])
def test_nonfinite_and_wrong_numerical_types_reject(bad):
    with pytest.raises(ValueError):
        packet.compare({'probability': bad}, {'probability': 0.0})


@pytest.mark.parametrize('field', ['preparation_abort_probability', 'continuation_abort_probability'])
@pytest.mark.parametrize('bad', [-1e-10, -5e-324])
def test_replay_tolerance_does_not_accept_negative_retained_outcomes(field, bad):
    # The retained display rounds these tiny abort observations to zero.
    # Closeness to that valid display must not admit a negative probability.
    actual = strict_json.load(packet.OUTPUT)
    changed = deepcopy(actual)
    records = changed['number_conserving_detector']['complete_source_readout']['complete_records']
    assert records[0][field] == 0.0
    records[0][field] = bad
    with pytest.raises(ValueError, match='outcome probability'):
        packet.compare(changed, actual)
    packet.compare(actual, actual)


@pytest.mark.parametrize('actual,bad', [(0., -1e-10), (1., 1.+1e-10)])
def test_pointer_probability_domain_is_stricter_than_replay_distance(actual, bad):
    with pytest.raises(ValueError, match='outcome probability'):
        packet.compare({'pointer_probabilities': [bad]}, {'pointer_probabilities': [actual]})


@pytest.mark.parametrize('bad', [{}, {'a': [1, 2]}, {'a': [1], 'extra': 0}, {'a': [True]}])
def test_missing_extra_and_wrong_typed_payloads_reject(bad):
    with pytest.raises(ValueError):
        packet.compare(bad, {'a': [1]})


def test_source_snapshot_detects_changed_loaded_code(monkeypatch):
    pins = packet.source_bindings()
    pins['model/particle_transport/source.py'] = '0'*64
    monkeypatch.setattr(packet, 'source_bindings', lambda: pins)
    with pytest.raises(ValueError, match='fresh interpreter'):
        packet.report()


def test_warm_compiler_source_is_revalidated(monkeypatch):
    packet.compiler._source()
    original = Path.read_text
    target = (packet.ROOT/packet.compiler.RESPONSE).resolve()

    def changed(path, *args, **kwargs):
        text = original(path, *args, **kwargs)
        if path.resolve() == target:
            value = json.loads(text)
            value['generators'] = value['generators'][:-1]
            return json.dumps(value)
        return text

    with monkeypatch.context() as patch:
        patch.setattr(Path, 'read_text', changed)
        with pytest.raises(ValueError, match='twelve generators'):
            packet.report()
    packet.compiler._source.cache_clear()


def test_warm_native_source_is_revalidated(monkeypatch):
    packet.slow_geometry.source()
    original = Path.read_text
    target = packet.slow_geometry.SPEC.resolve()

    def changed(path, *args, **kwargs):
        text = original(path, *args, **kwargs)
        if path.resolve() == target:
            value = json.loads(text)
            value['source_sha256']['oph_exact/carrier.py'] = '0'*64
            return json.dumps(value)
        return text

    with monkeypatch.context() as patch:
        patch.setattr(Path, 'read_text', changed)
        with pytest.raises(ValueError, match='pinned native geometry source'):
            packet.report()
    packet.slow_geometry.source.cache_clear()


def test_closure_includes_eager_parents_and_actual_source_data():
    pins = packet.source_bindings()
    for name in (
        'model/__init__.py', 'model/particle_transport/__init__.py',
        'model/particle_transport/algebra.py', 'model/particle_transport/source.py',
        'model/particle_transport/detector.py', 'model/particle_transport/test_source.py',
        'model/particle_transport/test_detector.py', 'model/particle_transport/test_algebra.py',
        'model/particle_response/dimer.py', 'model/particle_response/spectrum.py',
        'model/full_fock_field/field.py', 'model/field_thermal/thermal.py',
        'model/observer/source_scalar/compiler.py', packet.compiler.RESPONSE,
        packet.compiler.CONTROL_RECEIPT,
        'evidence/native_geometric_source_20260925/spec.json',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/carrier.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_exact/__init__.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/self_readback_repair_closure.py',
        'evidence/observer_dynamics_20260925/oph-physics-sim/oph_fpe/dynamics/__init__.py',
    ):
        assert name in pins
    assert 'model/particle_transport/results.json' not in pins


def test_duplicate_json_keys_reject():
    with pytest.raises(ValueError):
        strict_json.loads('{"physical_particle_identification": false, "physical_particle_identification": true}')
