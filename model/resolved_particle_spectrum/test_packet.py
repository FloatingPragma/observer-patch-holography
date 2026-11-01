from copy import deepcopy
import json

import pytest

from . import packet


@pytest.fixture(scope='module')
def receipt(): return packet._produce()


def test_json_roundtrip_retains_all_integer_powers_and_inferred_weights(receipt):
    restored = json.loads(json.dumps(receipt, allow_nan=False))
    packet.compare(receipt, restored); packet.validate(restored)


@pytest.mark.parametrize('change', ['matrix', 'count', 'scope', 'weight'])
def test_tampered_matrix_count_scope_or_inferred_weight_is_rejected(receipt, change):
    value = deepcopy(receipt); result = value['executed_accelerated_resolution']
    if change == 'matrix': result['retained_dyadic_power_matrices'][1]['real_integer_numerators'][0][0] += 1
    elif change == 'count': result['settings'][1]['pointer_settings'][0]['pointer_zero_count'] += 1
    elif change == 'scope': result['actual_raw_native_long_word_executed'] = True
    else: result['weight_reconstruction']['weights'][0]['real_estimate'] = '1/10'
    with pytest.raises(ValueError): packet.validate(value)


def test_in_run_source_edit_fails_before_cached_computation(monkeypatch):
    monkeypatch.setattr(packet, 'source_bindings', lambda: {})
    with pytest.raises(ValueError, match='source changed since import'): packet._produce()
