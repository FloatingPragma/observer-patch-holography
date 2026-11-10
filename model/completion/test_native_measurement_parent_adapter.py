"""The parent interface cannot admit scientific or unreviewed source changes."""
from copy import deepcopy
import hashlib
import pytest
from . import native_measurement_parent_adapter as adapter
from model.native_energy_measurement import packet


def test_actual_parent_authenticates_retained_bytes_and_exact_live_science():
    assert adapter.authenticate(packet.OUTPUT, adapter.RETAINED_SHA256,
                                packet.source_bindings()) == adapter.REPLAYED_SHA256


@pytest.mark.parametrize('damage', ['science', 'type', 'scope', 'extra', 'delete',
                                   'wrong_append', 'stale_source'])
def test_every_unreviewed_payload_or_source_transition_fails(damage):
    old = adapter.strict_json.load(adapter.ROOT / adapter.ARCHIVE)
    current = adapter.strict_json.load(packet.OUTPUT)
    sources = deepcopy(current['source_bindings'])
    if damage == 'science':
        current['all_native_words'][0]['stage'] += 1
    elif damage == 'type':
        current['all_native_words'][0]['stage'] = False
    elif damage == 'scope':
        current['scope']['rounded_instrument_is_exactly_trace_preserving'] = True
    elif damage == 'extra':
        current['source_bindings']['tools/strict_json.py/extra'] = 'f' * 64
    elif damage == 'delete':
        del current['source_bindings'][next(iter(old['source_bindings']))]
    elif damage == 'wrong_append':
        current['source_bindings'][next(iter(adapter.REVIEWED_APPEND))] = 'f' * 64
    else:
        sources[next(iter(sources))] = 'f' * 64
    with pytest.raises(ValueError):
        adapter.compare_payloads(old, current, sources)


@pytest.mark.parametrize('target', ['archive', 'current', 'source'])
def test_archive_current_and_source_byte_guards(tmp_path, target):
    sources = packet.source_bindings()
    for relative in [adapter.ARCHIVE, adapter.PARENT, *sources]:
        destination = tmp_path / relative
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes((adapter.ROOT / relative).read_bytes())
    relative = {'archive': adapter.ARCHIVE, 'current': adapter.PARENT,
                'source': next(iter(sources))}[target]
    with (tmp_path / relative).open('ab') as stream:
        stream.write(b'\n')
    with pytest.raises(ValueError):
        adapter.authenticate(tmp_path / adapter.PARENT, adapter.RETAINED_SHA256,
                             sources, root=tmp_path)
