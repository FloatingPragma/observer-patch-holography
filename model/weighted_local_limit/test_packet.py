"""Deterministic replay, source custody and scientific-promotion controls."""
from copy import deepcopy
from fractions import Fraction as F

import pytest

from . import packet


@pytest.fixture(scope='module')
def result():
    return packet.produce()


def test_retained_receipt_replays(result):
    assert packet.canonical(packet.strict_json.load(packet.OUTPUT))==packet.canonical(result)


def test_source_closure_binds_own_tests_proof_and_parent(result):
    expected={p.relative_to(packet.ROOT).as_posix() for p in packet.HERE.iterdir()
              if p.suffix in ('.py','.md')}
    assert expected<=result['source_bindings'].keys()
    assert 'model/weighted_source_family/results.json' in result['source_bindings']
    assert 'Lean/Thermodynamics/FiniteConditionalRepair.lean' in result['source_bindings']
    assert not any(name.startswith('docs/') for name in result['source_bindings'])


@pytest.mark.parametrize('change',['physical','discard_sectors','bound','numeric_type'])
def test_changed_science_or_scalar_type_is_rejected(tmp_path,result,monkeypatch,change):
    changed=deepcopy(result)
    if change=='physical':changed['scope']['source_selected_capacity_or_density']=True
    elif change=='discard_sectors':changed['all_domain_theorem']['all_local_number_sectors_retained']=False
    elif change=='bound':changed['all_domain_theorem']['conditioning_probability_lower']='1/B'
    else:changed['family']['declared_mean_occupations_per_cell'][0]=2.0
    output=tmp_path/'changed.json'; output.write_text(packet.canonical(changed))
    monkeypatch.setattr(packet,'produce',lambda:result)
    with pytest.raises(ValueError,match='differs from exact replay'):
        packet.verify(output)


def test_parent_source_custody_is_not_silently_refreshed(monkeypatch):
    original=packet.strict_json.load
    def stale(path):
        value=original(path)
        if path==packet.parent.OUTPUT:
            value['source_bindings']['model/weighted_source_family/family.py']='0'*64
        return value
    monkeypatch.setattr(packet.strict_json,'load',stale)
    with pytest.raises(ValueError,match='source receipt is stale'):
        packet.source_bindings()


def test_duplicate_json_key_rejected(tmp_path):
    output=tmp_path/'duplicate.json'; output.write_text('{"schema":"a","schema":"b"}')
    with pytest.raises(ValueError):packet.verify(output)


def test_receipt_rounding_is_outward():
    lo,hi=packet.rounded_interval((F(1,7),F(2,7)))
    assert F(lo)<=F(1,7)<F(2,7)<=F(hi)
