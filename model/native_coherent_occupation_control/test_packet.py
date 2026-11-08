from copy import deepcopy
from fractions import Fraction as F
import pytest
from tools import strict_json
from . import packet


@pytest.fixture(scope='module')
def fresh():return packet.produce()


def test_retained_receipt_matches_fresh_owner(fresh):
    assert packet.check_payload(strict_json.load(packet.OUTPUT),fresh)


@pytest.mark.parametrize('field',('source_bindings','exact_six_level_exchange_generator',
    'exact_six_level_phase_generators','complete_native_even_parity','executed_addressed_Bell_services',
    'executed_diagonal_services','literal_word_angle_obstruction','literal_full_source_control_obstruction',
    'complete_22_pulse_service_bounds','ideal_balanced_retry_budget','scope'))
def test_each_scientific_or_custody_section_is_bound(fresh,field):
    changed=deepcopy(fresh);changed[field]={}
    with pytest.raises(ArithmeticError):packet.check_payload(changed,fresh)


def test_scalar_types_are_not_coerced_and_nonfinite_values_are_refused():
    for left,right in ((1,True),(1,1.0),(0,False),('1',1)):
        if type(right)is float:
            with pytest.raises(TypeError):packet.canonical({'x':right})
        else:
            with pytest.raises(ArithmeticError):packet.check_payload({'x':left},{'x':right})
    for value in (float('nan'),float('inf'),float('-inf')):
        with pytest.raises(TypeError):packet.canonical({'x':value})
    assert packet.canonical({'x':(F(1,2),2)})==packet.canonical({'x':['1/2',2]})


def test_duplicate_json_fields_are_rejected(tmp_path):
    value=tmp_path/'receipt.json';value.write_text('{"schema":"a","schema":"b"}')
    with pytest.raises(ValueError):strict_json.load(value)


def test_owner_pins_all_local_sources_and_exact_processor_law(fresh):
    paths=fresh['source_bindings']
    for name in ('derived_control.py','diagonal_control.py','literal_controls.py','test_packet.py','PROOF.md'):
        assert 'model/native_coherent_occupation_control/'+name in paths
    assert 'model/certified_reference/source.py'in paths
    assert 'model/certified_reference/interval.py'in paths
    assert packet.beam.compiler.RESPONSE in paths
    assert packet.beam.compiler.CONTROL_RECEIPT in paths
