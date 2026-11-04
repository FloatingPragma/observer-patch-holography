import copy
import hashlib
from pathlib import Path

import pytest

from tools import strict_json
from . import packet


def test_cold_complete_observation_receipt_and_source_custody():
    actual = packet.build()
    stored = strict_json.load(packet.OUTPUT)
    packet.validate(stored, actual)
    assert 'model/moving_contact_scattering/PROOF.md' in stored['source_bindings']
    assert 'model/native_spin_observer/precision.py' in stored['source_bindings']
    for path, digest in stored['source_bindings'].items():
        assert hashlib.sha256((packet.ROOT/path).read_bytes()).hexdigest() == digest
    for key in ('analytic_scattering_phase_used_as_native_gate',
                'fine_spacing_native_program_executed', 'physical_particle_identification',
                'many_body_ground_vacuum_constructed', 'relativistic_interacting_QFT_constructed',
                'Einstein_dynamics_constructed', 'completion_gate_accepted'):
        wrong = copy.deepcopy(stored)
        wrong['scope'][key] = True
        with pytest.raises(ValueError):
            packet.validate(wrong, actual)
    wrong = copy.deepcopy(stored)
    wrong['executed_small_native_controls'][0]['complete_pointer_slots'] -= 1
    with pytest.raises(ValueError):
        packet.validate(wrong, actual)
