from copy import deepcopy
import subprocess
import sys

import pytest

from . import packet


@pytest.fixture(scope='module')
def retained():
    return packet.build()


def test_complete_retained_receipt_and_native_execution(retained):
    assert packet.verify()
    assert retained['complete_native_execution']['first_failure_slots'] == 172368
    assert retained['complete_native_execution']['decoded_core_executed']


@pytest.mark.parametrize('field', ['actual_protected_topological_label_attached',
                                  'reflection_coordinate_is_physical_particle_position',
                                  'continuum_detector_convergence_proved_here',
                                  'full_encoded_twenty_bit_native_execution_claimed'])
def test_physical_or_execution_promotion_rejected(retained, field, tmp_path):
    value = deepcopy(retained)
    value['scope'][field] = True
    path = tmp_path/'forged.json'
    path.write_text(packet.canonical(value))
    with pytest.raises(ValueError):
        packet.verify(path)


def test_wrong_target_path_or_duplicate_key_cannot_reuse_default_receipt(tmp_path):
    path = tmp_path/'duplicate.json'
    path.write_text('{"schema": 1, "schema": 2}')
    result = subprocess.run([sys.executable, '-m', 'model.spin_operational_position',
                             '--check', str(path)], cwd=packet.ROOT,
                            capture_output=True, text=True, check=False)
    assert result.returncode != 0
    assert 'duplicate' in (result.stdout+result.stderr).lower()


def test_source_bindings_keep_code_proofs_and_witnesses(retained):
    names = retained['source_bindings']
    assert not any(x.startswith(('docs/', 'tracking/', 'claims/')) for x in names)
    assert {'model/spin_operational_position/reflectors.json',
            'model/native_ordered_ground/witness.npz',
            'model/native_spin_observer/precision.py',
            'model/observer/source_scalar/compiler.py'} <= names.keys()
