from copy import deepcopy
import pytest
from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def fresh():return packet.produce()


def test_retained_receipt_matches_fresh_replay(fresh):
    assert packet.check_payload(strict_json.load(packet.OUTPUT),fresh)


@pytest.mark.parametrize('field',('source_bindings','source_ground_dilution','literal_native_endpoint',
    'exact_density_ground_closure','necessary_visibility','finite_density_control','growing_norm_control','scope'))
def test_all_receipt_science_and_sources_are_bound(fresh,field):
    changed=deepcopy(fresh);changed[field]={}
    with pytest.raises(ArithmeticError):packet.check_payload(changed,fresh)


def test_scalar_types_are_distinct_and_nonfinite_values_refused():
    for left,right in ((1,True),(0,False),('1',1)):
        with pytest.raises(ArithmeticError):packet.check_payload({'x':left},{'x':right})
    for value in (1.0,float('nan'),float('inf'),float('-inf')):
        with pytest.raises(TypeError):packet.canonical({'x':value})


def test_duplicate_fields_rejected(tmp_path):
    path=tmp_path/'duplicate.json';path.write_text('{"source":1,"source":2}')
    with pytest.raises(ValueError):strict_json.load(path)


def test_native_closure_and_spatial_dictionary_sources_are_bound(fresh):
    paths=fresh['source_bindings']
    for path in ('model/primitive_observer/transfer.py','model/native_program/interpreter.py',
                 'model/federated_observer/SCALING.md','model/native_spatial_scaling/PROOF.md',
                 'model/native_dilute_ground_limit/test_packet.py'):
        assert path in paths
