"""Owning receipt and source-binding rejection tests."""
from copy import deepcopy

import pytest

from . import packet


@pytest.fixture(scope='module')
def result():
    return packet.produce()


def test_retained_receipt_replays(result):
    assert packet.strict_json.load(packet.OUTPUT) == result


def test_own_proof_and_test_fixtures_are_bound(result):
    expected = {p.relative_to(packet.ROOT).as_posix() for p in packet.HERE.iterdir()
                if p.suffix in ('.py', '.md')}
    assert expected <= result['source_bindings'].keys()
    assert not any(name.startswith('docs/') for name in result['source_bindings'])
    assert 'Lean/Thermodynamics/FiniteConditionalRepair.lean' in result['source_bindings']


def test_mutated_physical_promotion_rejected(tmp_path, result):
    changed = deepcopy(result)
    changed['scope']['physical_clock_or_capacity_selected'] = True
    file = tmp_path/'changed.json'
    file.write_text(packet.canonical(changed))
    with pytest.raises(ValueError, match='differs from exact replay'):
        packet.verify(file)


def test_duplicate_key_receipt_rejected(tmp_path):
    file = tmp_path/'duplicate.json'
    file.write_text('{"schema": "one", "schema": "two"}')
    with pytest.raises(ValueError):
        packet.verify(file)
