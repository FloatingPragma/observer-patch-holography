import copy
import numpy as np
import pytest

from . import packet


def test_cold_artifact_roundtrip_and_archive_tamper_rejection(native_packet,tmp_path):
    value,arrays=native_packet;target=tmp_path/'receipt.json'
    packet.write_pair(target,value,arrays)
    retained,stored=packet.read_pair(target)
    packet.compare(retained,value)
    assert stored.keys()==arrays.keys()
    for name in arrays:assert np.array_equal(stored[name],arrays[name])
    with target.with_suffix('.npz').open('ab') as stream:stream.write(b'corrupt')
    with pytest.raises(ValueError,match='archive hash'):packet.read_pair(target)


def test_missing_outcome_and_complete_state_are_rejected(native_packet):
    value,arrays=native_packet
    changed=copy.deepcopy(value);changed['executed_native_detector']['all_joint_selector_charge_clock_outcomes'].pop()
    with pytest.raises(ValueError,match='outcome coverage'):packet.validate_accounting(changed,arrays)
    changed_arrays=dict(arrays);changed_arrays['phase_clock_complete_success_state']=arrays['phase_clock_complete_success_state'][:-1]
    with pytest.raises(ValueError,match='quantum array'):packet.validate_accounting(value,changed_arrays)


def test_scope_promotions_are_rejected(native_packet):
    value,arrays=native_packet
    changed=copy.deepcopy(value);changed['executed_native_detector']['all78_lines_separately_resolved_by_native_bins']=True
    with pytest.raises(ValueError,match='scope promotion'):packet.validate_accounting(changed,arrays)
    changed=copy.deepcopy(value);changed['physical_identification']=True
    with pytest.raises(ValueError,match='integrated claim'):packet.validate_accounting(changed,arrays)
    changed=copy.deepcopy(value);changed['executed_native_detector']['source_records'][0]['first_failure_probability']=float('nan')
    with pytest.raises(ValueError,match='failure probability'):packet.validate_accounting(changed,arrays)


def test_current_source_bindings_and_changed_source_rejection(monkeypatch):
    bindings=packet.source_bindings()
    assert 'model/native_spectral_detector/contract.cpp' in bindings
    assert 'paper/tex_fragments/NATIVE_SPECTRAL_DETECTOR.tex' in bindings
    packet.authenticate_parents()
    monkeypatch.setattr(packet,'source_bindings',lambda:{**bindings,'changed':'0'*64})
    with pytest.raises(ValueError,match='source changed since import'):packet._produce()
