from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json

import pytest

from . import records as r,packet


def test_complete_original_connected_record_scores():
    rows=r.exact_scores()
    assert rows[0]['raw_score']==F(1,7170)
    assert rows[1]['normalized_response']==F(499,504)
    assert r.program(0,1,2)['declared_resources']['attempts_upper']==7


def test_wrong_reflection_cut_breaks_source_statistic(monkeypatch):
    original=r.local_records.scored_step
    def wrong(value,command,density,side=None):
        return original(value,command,density,'incoming' if side=='outgoing' else side)
    monkeypatch.setattr(r.local_records,'scored_step',wrong)
    with pytest.raises(ArithmeticError,match='reflected past'):
        r.exact_scores(0)


def test_global_attempt_discount_and_finite_resource_budget():
    data=r.abel(0,F(24),512)
    assert data['discount']==F(105,107)
    assert data['discarded_tail_upper']<F(1,10000)
    assert data['working_attempts_per_lag_program']==516
    assert r.abel(1,24,512)['discount']>data['discount']
    with pytest.raises(ValueError):r.program(0,r.source.MAX_ATTEMPTS-3)
    with pytest.raises(ValueError):r.abel(0,24,r.source.MAX_ATTEMPTS-3)
    with pytest.raises(ValueError):r.abel(21,24,1)
    with pytest.raises(ValueError):r.abel(0,24,100000)
    with pytest.raises(ValueError):r.abel(20,24,4096)


@pytest.mark.parametrize('bad',[True,-1,0.0,'2'])
def test_record_count_and_shift_types(bad):
    with pytest.raises(ValueError):r.program(0,bad)
    with pytest.raises(ValueError):r.exact_scores(bad)


@pytest.fixture(scope='module')
def receipt():return packet.produce()


def test_cold_retained_receipt(receipt):
    assert packet.canonical(json.loads(packet.OUTPUT.read_text()))==packet.canonical(receipt)
    assert receipt['controlled_exterior_dressing']['positive_energy_pole_continuation_certified'] is False


def test_transitive_source_pins_are_fresh(receipt):
    assert receipt['source_bindings']==packet.source_bindings()
    for name,digest in receipt['source_bindings'].items():
        assert not name.startswith('docs/')
        assert digest==hashlib.sha256((packet.ROOT/name).read_bytes()).hexdigest()
    assert 'model/native_program/interpreter.py' in receipt['source_bindings']
    assert 'model/native_contact_response/PROOF.md' in receipt['source_bindings']


@pytest.mark.parametrize('mutation',['bool','float','extra','pin','bound'])
def test_strict_typed_receipt_mutation(tmp_path,monkeypatch,receipt,mutation):
    value=deepcopy(receipt)
    if mutation=='bool':value['source']['occupation']=True
    elif mutation=='float':value['source']['occupation']=2.0
    elif mutation=='extra':value['unproved_particle']=True
    elif mutation=='pin':value['source_bindings'][next(iter(value['source_bindings']))]='0'*64
    else:value['connected_comparison_threshold']['all_level_upper']='0'
    path=tmp_path/'bad.json';path.write_text(packet.canonical(value))
    monkeypatch.setattr(packet,'produce',lambda:receipt)
    with pytest.raises(ValueError):packet.verify(path)


def test_duplicate_json_key_is_rejected(tmp_path):
    path=tmp_path/'bad.json';path.write_text('{"schema":1,"schema":1}')
    with pytest.raises(ValueError):packet.verify(path)
