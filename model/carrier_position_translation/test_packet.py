from copy import deepcopy
import subprocess
import sys

import pytest

from . import packet


def test_exact_retained_receipt():
    assert packet.verify()


@pytest.mark.parametrize('field', ['field_translation_generator_joined_to_carrier_position',
                                  'carrier_total_load_is_protected_spin_position',
                                  'universal_particle_Heisenberg_bound_from_native_records'])
def test_false_physical_promotion_fails(field, tmp_path):
    result = deepcopy(packet.build())
    result['scope'][field] = True
    target = tmp_path/'forged.json'
    target.write_text(packet.canonical(result))
    with pytest.raises(ValueError):
        packet.verify(target)


def test_cli_checks_requested_receipt_and_rejects_duplicate_fields(tmp_path):
    target = tmp_path/'duplicate.json'
    target.write_text('{"schema": 1, "schema": 2}')
    run = subprocess.run([sys.executable, '-m', 'model.carrier_position_translation',
                          '--check', str(target)], cwd=packet.ROOT,
                         capture_output=True, text=True, check=False)
    assert run.returncode != 0
    assert 'duplicate' in (run.stdout+run.stderr).lower()


def test_source_bindings_exclude_editorial_docs_and_include_owners():
    bindings = packet.bindings()
    assert not any(name.startswith(('docs/', 'claims/', 'tracking/')) for name in bindings)
    assert {'model/common_source/slow_geometry.py', 'model/observer/native.py',
            'model/protected_spin_limit/algebra.py',
            'model/defect_field_bridge/full_domain.py'} <= bindings.keys()
    assert all(len(value) == 64 and (packet.ROOT/name).is_file()
               for name, value in bindings.items())
