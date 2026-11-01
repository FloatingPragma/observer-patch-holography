from copy import deepcopy
import json

import numpy as np
import pytest

from . import packet


def test_lossless_archive_replays_full_states_and_every_failure_table(readout, tmp_path):
    value, arrays = readout
    target = tmp_path/'retained.json'
    packet.write_pair(target, value, arrays)
    retained, restored = packet.read_pair(target)
    packet.compare(value, retained)
    for name in arrays: np.testing.assert_array_equal(arrays[name], restored[name])
    with target.with_suffix('.npz').open('ab') as stream: stream.write(b'changed')
    with pytest.raises(ValueError, match='hash changed'): packet.read_pair(target)


def test_omitted_zero_repair_slot_and_missing_abort_array_are_rejected(readout):
    value, arrays = readout
    altered = deepcopy(value)
    rows = altered['executed_native_readout']['all130_repair_records']
    rows.pop(next(j for j, row in enumerate(rows) if row['probability'] == 0))
    with pytest.raises(ValueError, match='repair records'): packet.validate_accounting(altered, arrays)
    reduced = dict(arrays); reduced.pop('pointer_failure_weights')
    with pytest.raises(ValueError, match='archive key'): packet.validate_accounting(value, reduced)


def test_incorrect_scope_tolerance_nonfinite_or_source_are_rejected(readout):
    value, _ = readout
    for update in ('scope', 'bound', 'nan', 'source'):
        altered = deepcopy(value)
        if update == 'scope': altered['executed_native_readout']['resolved_all_eight_frequencies_from_these_two_settings'] = True
        elif update == 'bound': altered['replay_float_absolute_tolerance'] *= 2
        elif update == 'nan': altered['executed_native_readout']['terminal_native_repair_probability'] = float('nan')
        else: altered['source_bindings']['model/protected_particle_readout/experiment.py'] = '0'*64
        with pytest.raises(ValueError): packet.compare(altered, value)


def test_in_run_source_changes_fail_closed_before_expensive_execution(monkeypatch):
    monkeypatch.setattr(packet, 'source_bindings', lambda: {})
    with pytest.raises(ValueError, match='source changed since import'): packet._produce()
