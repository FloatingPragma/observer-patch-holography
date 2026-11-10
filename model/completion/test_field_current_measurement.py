"""The successor changes custody service only, never consumer science."""
import copy
import pytest
from tools import strict_json
from . import field_current_measurement as current


def test_exact_complete_parent_and_private_consumer_namespace():
    from model.native_energy_measurement.packet import source_bindings
    assert current.authenticate(current.ROOT/current.retained.PARENT,current.retained.RETAINED_SHA256,source_bindings())==current.CURRENT_SHA256
    for name in current.SOURCES:
        module=current.owner(name)
        assert module.produce.__code__.co_filename==str(current.ROOT/current.SOURCES[name][0])
        assert module.parent_custody.authenticate is current.authenticate
        import importlib
        original=importlib.import_module('model.'+name+'.packet')
        assert original.parent_custody is current.retained
        assert original.produce.__code__.co_code==module.produce.__code__.co_code


@pytest.mark.parametrize('wrong', ['parent','retained','live_source','science'])
def test_current_custody_rejects_changes(wrong):
    from model.native_energy_measurement.packet import source_bindings
    parent=current.ROOT/current.retained.PARENT;old=current.retained.RETAINED_SHA256;pins=source_bindings()
    if wrong=='parent':parent=current.ROOT/'model/native_program/results.json'
    elif wrong=='retained':old='0'*64
    elif wrong=='live_source':pins[next(iter(pins))]='0'*64
    else:
        retained=strict_json.load(current.ROOT/current.retained.ARCHIVE)
        value=strict_json.load(parent);value['scope']['rounded_instrument_is_exactly_trace_preserving']=True
        with pytest.raises(ValueError):current.retained.compare_payloads(retained,value,pins)
        return
    with pytest.raises(ValueError):current.authenticate(parent,old,pins)
