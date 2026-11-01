from copy import deepcopy

import numpy as np
import pytest

from . import packet


def test_lossless_archive_and_hash_authentication(ground_readout,tmp_path):
    value,arrays=ground_readout
    target=tmp_path/'readout.json'
    packet.write_pair(target,value,arrays)
    restored,stored=packet.read_pair(target)
    packet.compare(value,restored)
    for key in arrays:np.testing.assert_array_equal(stored[key],arrays[key])
    with target.with_suffix('.npz').open('ab') as stream:stream.write(b'changed')
    with pytest.raises(ValueError,match='archive hash'):packet.read_pair(target)


def test_every_repair_and_export_slot_is_required(ground_readout):
    value,arrays=ground_readout
    changed=deepcopy(value)
    changed['executed_native_readout']['all130_repair_outcome_slots'].pop(0)
    with pytest.raises(ValueError,match='repair records'):packet.validate_accounting(changed,arrays)
    reduced=dict(arrays);reduced.pop('ground_tree_failure_weights')
    with pytest.raises(ValueError,match='archive key'):packet.validate_accounting(value,reduced)


def test_unproved_scope_changed_source_and_nonfinite_evidence_are_rejected(ground_readout):
    value,arrays=ground_readout
    changed=deepcopy(value)
    changed['executed_native_readout']['all78_lines_native_resolved']=True
    with pytest.raises(ValueError,match='scope promoted'):packet.validate_accounting(changed,arrays)
    changed=deepcopy(value)
    changed['source_bindings']['model/ground_spectral_source/experiment.py']='0'*64
    with pytest.raises(ValueError):packet.compare(changed,value)
    changed=deepcopy(value)
    changed['executed_native_readout']['settings'][0]['all_six_pointer_probabilities'][0]=float('nan')
    with pytest.raises(ValueError,match='pointer probabilities'):packet.validate_accounting(changed,arrays)


def test_in_run_source_mutation_fails_before_execution(monkeypatch):
    monkeypatch.setattr(packet,'source_bindings',lambda:{})
    with pytest.raises(ValueError,match='source changed since import'):packet._produce()
