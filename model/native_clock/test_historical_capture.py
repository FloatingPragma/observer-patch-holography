"""Historical timing provenance is preserved; compatibility is tightly scoped."""
from copy import deepcopy
import hashlib
from pathlib import Path

import pytest

from . import historical_capture as history
from . import host_clock as host


def fixture():
    raw = host.strict_json.load(host.RAW_PATH)
    return raw, host._hashes()


def test_exact_original_capture_and_environment_are_preserved():
    raw, current = fixture()
    assert hashlib.sha256(host.RAW_PATH.read_bytes()).hexdigest() == history.CAPTURE_SHA256
    result = history.source_compatibility(raw, current, host.ROOT, host.canonical)
    assert result['original_capture_data_match'] is True
    assert result['capture_source_manifest_preserved'] == raw['capture_sources']
    assert result['capture_declared_requirements_sha256'] == history.ORIGINAL_REQUIREMENTS_SHA256
    assert result['current_analysis_requirements_sha256'] == history.REVIEWED_REQUIREMENTS_SHA256
    assert result['current_timing_measurement_claimed'] is False
    assert result['unrecorded_runtime_or_hardware_reconstructed'] is False
    assert result['timing_portability_claimed'] is False


@pytest.mark.parametrize('mutation', ['manifest', 'requirements', 'native'])
def test_compatibility_rejects_unreviewed_source_changes(mutation):
    raw, current = fixture()
    if mutation == 'manifest':
        raw['capture_sources'][history.REQUIREMENTS] = '0'*64
    elif mutation == 'requirements':
        current[history.REQUIREMENTS] = '0'*64
    else:
        current['model/observer/native.py'] = '0'*64
    with pytest.raises(ValueError, match='capture source bindings'):
        history.source_compatibility(raw, current, host.ROOT, host.canonical)


def test_changed_counterfactual_is_not_labeled_original_measurement():
    raw, current = fixture()
    raw['runs'][0]['block_endpoints'][-1]['elapsed_ns'] += 1
    result = history.source_compatibility(raw, current, host.ROOT, host.canonical)
    assert result['original_capture_data_match'] is False
    assert result['modified_historical_input_is_a_measurement'] is False


@pytest.mark.parametrize('old,new', [('moves += abs(current)', 'moves += 2*abs(current)'),
                                    ('random.Random(SEED)', 'random.Random(SEED+1)'),
                                    ('clock()-started', 'clock()-started+1')])
def test_capture_scheduler_and_timer_ast_changes_are_rejected(tmp_path, old, new):
    raw, current = fixture()
    for name in (history.SOURCE_SNAPSHOT, history.REQUIREMENTS_SNAPSHOT):
        target = tmp_path/name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes((host.ROOT/name).read_bytes())
    source = (host.ROOT/history.SOURCE).read_text()
    assert old in source
    (tmp_path/history.SOURCE).write_text(source.replace(old, new, 1))
    with pytest.raises(ValueError, match='syntax changed'):
        history.source_compatibility(raw, current, tmp_path, host.canonical)
