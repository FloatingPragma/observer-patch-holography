"""Replay the exact finite interfaces and retain their physical scope."""
import pytest

from tools import strict_json
from .completion_certificate import OUTPUT, _same_json, report


def test_complete_interface_receipt_replays():
    assert _same_json(strict_json.load(OUTPUT), report())


@pytest.mark.parametrize("changed", [False, 0.0])
def test_receipt_equality_preserves_scalar_types(changed):
    assert not _same_json({"physical_predictions": 0},
                          {"physical_predictions": changed})


def test_receipt_equality_rejects_nonfinite_values():
    with pytest.raises(ValueError):
        _same_json({"error": float("nan")}, {"error": float("nan")})


def test_fair_interface_receipts_are_freshly_computed(tmp_path,monkeypatch):
    from types import SimpleNamespace
    from . import completion_certificate as owner
    calls=[];producers=[]
    for name,_ in owner.FAIR_LOCAL_INTERFACES:
        value={'schema':name,'all_sectors':list(range(13)),'failure_retained':True}
        path=tmp_path/(name+'.json');path.write_text(__import__('json').dumps(value))
        def produce(value=value,name=name):
            calls.append(name)
            return value
        producers.append((name,SimpleNamespace(OUTPUT=path,produce=produce)))
    monkeypatch.setattr(owner,'FAIR_LOCAL_INTERFACES',tuple(producers))
    result=owner._replay_fair_interfaces()
    assert list(result)==calls==[name for name,_ in producers]


@pytest.mark.parametrize('change',('float_integer','bool_integer','missing_number','missing_failure','stale_value'))
def test_fair_interface_rejects_changed_owner_receipt(tmp_path,monkeypatch,change):
    import json
    from types import SimpleNamespace
    from . import completion_certificate as owner
    expected={'attempts':1931,'all_sectors':list(range(13)),'failure_retained':True}
    changed=json.loads(json.dumps(expected))
    if change=='float_integer':changed['attempts']=1931.0
    elif change=='bool_integer':changed['all_sectors'][1]=True
    elif change=='missing_number':changed['all_sectors'].pop()
    elif change=='missing_failure':del changed['failure_retained']
    else:changed['attempts']=1930
    path=tmp_path/'results.json';path.write_text(json.dumps(changed))
    producer=SimpleNamespace(OUTPUT=path,produce=lambda:expected)
    monkeypatch.setattr(owner,'FAIR_LOCAL_INTERFACES',(('fair',producer),))
    with pytest.raises(ValueError,match='receipt differs'):owner._replay_fair_interfaces()


def test_fair_interface_rejects_duplicate_receipt_keys(tmp_path,monkeypatch):
    from types import SimpleNamespace
    from . import completion_certificate as owner
    path=tmp_path/'results.json';path.write_text('{"all_sectors":true,"all_sectors":false}')
    monkeypatch.setattr(owner,'FAIR_LOCAL_INTERFACES',(('fair',SimpleNamespace(
        OUTPUT=path,produce=lambda:{'all_sectors':False})),))
    with pytest.raises(ValueError):owner._replay_fair_interfaces()


def test_nested_fair_replays_and_sources_are_registered_for_model_verify():
    from model import verify
    from . import completion_certificate as owner
    bindings=owner.source_bindings()
    assert ('model.particle_identification.completion_certificate',None,
            'model/particle_identification/completion_results.json') in verify.REPLAYS
    for name,producer in owner.FAIR_LOCAL_INTERFACES:
        assert name.startswith('fair_source_local_')
        assert producer.OUTPUT.relative_to(owner.ROOT).as_posix() in bindings
        assert set(producer.source_bindings())<=set(bindings)
