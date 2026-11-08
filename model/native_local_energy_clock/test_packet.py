"""Source and complete-observation mutations fail closed."""
import copy

import pytest

from tools import strict_json
from . import energy, packet


@pytest.fixture(scope='module')
def replayed():
    return packet.produce()


def test_sparse_graph_counts_do_not_allocate_the_port_carrier():
    assert energy.graph_data(10**122,((0,1),(1,10**122-1)))==(2,2,2)
    with pytest.raises(ValueError,match='at most 24'):
        energy.matrices(10**122,((0,1),),1)


def test_current_complete_receipt_replays(replayed,monkeypatch):
    monkeypatch.setattr(packet,'produce',lambda:copy.deepcopy(replayed))
    assert packet.verify()


@pytest.mark.parametrize('mutation',('boolean_integer','missing_zero','heat','source','physical_claim'))
def test_changed_retained_receipt_rejected(tmp_path,mutation,replayed,monkeypatch):
    monkeypatch.setattr(packet,'produce',lambda:copy.deepcopy(replayed))
    value=strict_json.load(packet.OUTPUT)
    if mutation=='boolean_integer':value['exact_H66']['dimension']=True
    elif mutation=='missing_zero':
        value['complete_fixed_finite_controlled_history_bounds'][0]['complete_outcome_labels'][0].pop()
    elif mutation=='heat':
        value['exact_H66']['repair_exchange_energy_change_on_P_reference']='0'
    elif mutation=='source':value['source_bindings'].pop(next(iter(value['source_bindings'])))
    else:value['scope']['symmetric_physical_stress_supplied']=True
    path=tmp_path/'receipt.json';path.write_text(packet.canonical(value))
    with pytest.raises(ValueError,match='differs from replay'):packet.verify(path)


def test_changed_measurement_owner_bound_rejected(monkeypatch):
    original=packet.strict_json.load
    def changed(path):
        value=copy.deepcopy(original(path))
        if path==packet.measurement_owner.OUTPUT:
            value['complete_instrument_budget']['native_words']=True
        return value
    monkeypatch.setattr(packet.strict_json,'load',changed)
    with pytest.raises(ValueError,match='owning complete measurement'):packet.owning_measurement()


@pytest.mark.parametrize('field',('word_sha256','distinct_primitive_checks','resources'))
def test_changed_measurement_execution_witness_rejected(tmp_path,monkeypatch,field):
    value=strict_json.load(packet.measurement_owner.OUTPUT)
    value['all_native_words'][0][field]=False
    path=tmp_path/'measurement.json';path.write_text(packet.canonical(value))
    monkeypatch.setattr(packet.measurement_owner,'OUTPUT',path)
    with pytest.raises(ValueError,match='audited complete measurement receipt bytes'):
        packet.owning_measurement()


def test_source_changes_during_replay_rejected(monkeypatch,replayed):
    source=packet.source_bindings();calls=0
    def changing():
        nonlocal calls
        calls+=1
        return source if calls==1 else dict(source,changed='0'*64)
    monkeypatch.setattr(packet,'source_bindings',changing)
    monkeypatch.setattr(packet,'compute',lambda:{'schema':replayed['schema']})
    monkeypatch.setattr(packet,'source_controls',lambda:replayed['source_controls'])
    monkeypatch.setattr(packet,'owning_measurement',lambda:replayed['owning_measurement_attachment'])
    with pytest.raises(ArithmeticError,match='source bytes changed'):packet.produce()


def test_duplicate_and_nonfinite_receipt_values_rejected(tmp_path):
    for contents in ('{"schema":0,"schema":1}','{"value":NaN}'):
        path=tmp_path/'receipt.json';path.write_text(contents)
        with pytest.raises(ValueError):packet.verify(path)
