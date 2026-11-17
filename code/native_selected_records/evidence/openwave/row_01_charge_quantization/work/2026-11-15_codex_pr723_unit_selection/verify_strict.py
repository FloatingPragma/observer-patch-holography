"""Independent protocol/receipt hardening; imports only the independent verifier."""
import gzip
import hashlib
import json
from pathlib import Path
import random
import sys
from verify import verify

HERE = Path(__file__).resolve().parent


def check(directory):
    edges = json.loads((HERE / 'seams.json').read_text())
    receipt = json.loads((directory / 'RECEIPT.json').read_text())
    raw = (directory / 'EVIDENCE.json.gz').read_bytes()
    data = json.loads(gzip.decompress(raw))
    assert hashlib.sha256(raw).hexdigest() == receipt['evidence_sha256']
    assert len(edges) == 30 and len({tuple(e) for e in edges}) == 30
    assert all(sum(v in e for e in edges) == 5 for v in range(12))
    assert [s['total'] for s in data['sectors']] == list(range(5))
    for sector in data['sectors']:
        size = len(sector['states'])
        for row in sector['branches']:
            assert len(row) == 5 and all(type(v) is int for v in row)
            i, e, coin, current, target = row
            assert 0 <= i < size and 0 <= target < size
            assert 0 <= e < 30 and coin in (0, 1)
    families = ('double', 'quadruple', 'signed', 'shifted')
    assert len(data['histories']) == 96
    for position, record in enumerate(data['histories']):
        trial, control = divmod(position, 2)
        port, family_index = divmod(trial, 4)
        family = families[family_index]
        initial = [2 if family == 'shifted' else 0] * 12
        if family == 'signed':
            initial[port], initial[(port + 1) % 12] = 8, -6
        else:
            initial[port] += 4 if family == 'quadruple' else 2
        assert record['initial_audit_only'] == initial
        assert record['family'] == family
        assert record['seed'] == 72301000 + trial
        assert record['law'] == ('swap_hold' if control else 'native')
        assert record['horizon'] == 384 and record['probe_cap'] == 32
    detector_inputs = [[2] + [0] * 11, [1, 1] + [0] * 10,
                       [0] * 12, [2] * 12]
    assert len(data['detector_trials']) == 4
    for i, record in enumerate(data['detector_trials']):
        assert record['initial_audit_only'] == detector_inputs[i]
        assert record['seed'] == 72301900 + i
        assert record['horizon'] == 0 and record['probe_cap'] == 32
        assert record['law'] == 'native'
    for record in data['histories'] + data['detector_trials']:
        rng = random.Random(record['seed'])
        for a, b, coin, current in record['body']:
            assert [a, b] == edges[rng.randrange(30)]
            assert coin == rng.randrange(2)
        for row in record['probe']:
            assert len(row) == 4 and all(type(v) is int for v in row)
            a, b, coin, current = row
            assert 0 <= a < 12 and 0 <= b < 12 and a != b
            assert [a, b] in edges or [b, a] in edges
            assert coin in (0, 1) and coin == rng.randrange(2)
        assert record['probe_outcome'] in ('conclusive', 'inconclusive')
        assert (record['decoded_preprobe_relative'] is None) == (
            record['probe_outcome'] == 'inconclusive')
    forced = data['forced_timeout']
    assert forced['initial_audit_only'] == detector_inputs[0]
    assert forced['seed'] == 72301999 and forced['horizon'] == 0
    assert forced['probe_cap'] == 2 and len(forced['probe']) == 2
    assert all(row[2] == 0 for row in forced['probe'])
    result = verify(directory)
    census = {
        'states': sum(len(s['states']) for s in data['sectors']),
        'native_branches': sum(len(s['branches']) for s in data['sectors']),
        'binary_states': sum(max(x) <= 1 for s in data['sectors'] for x in s['states']),
        'binary_exit_branches': sum(max(s['states'][i]) <= 1 and
            max(s['states'][target]) > 1 for s in data['sectors']
            for i, _, _, _, target in s['branches']),
        'unreachable_states': sum(d < 0 for s in data['sectors']
            for d in s['distance_to_binary']),
        'max_shortest_entry_events': max(d for s in data['sectors']
            for d in s['distance_to_binary'])}
    native = [r for r in data['histories'] if r['law'] == 'native']
    control = [r for r in data['histories'] if r['law'] == 'swap_hold']
    histories = {
        'native': len(native), 'changed_law': len(control),
        'native_balanced_within_horizon': sum(r['first_balanced_body_event']
            is not None for r in native),
        'changed_law_balanced_within_horizon': sum(r['first_balanced_body_event']
            is not None for r in control),
        'inconclusive_terminal_readings': sum(r['probe_outcome'] != 'conclusive'
            for r in data['histories']),
        'postentry_binary_exits': sum(any(span > 1 for span in
            r['body_spans_from_records'][r['first_balanced_body_event']:])
            for r in native if r['first_balanced_body_event'] is not None),
        'body_events': sum(len(r['body']) for r in data['histories']),
        'native_terminal_probe_events': sum(len(r['probe']) for r in data['histories'])}
    assert receipt['exact_census'] == census
    assert receipt['histories'] == histories
    assert receipt['double_initial_entry_probability'] == '5/30'
    assert receipt['double_initial_mean_entry_events'] == 6
    assert receipt['candidate'] == 'A-native-multiload-unit-selection-v1'
    assert receipt['official_row_status'] == 'untested'
    assert receipt['full_row_solved'] is False
    assert receipt['coherent_memory_used'] is False
    assert receipt['nature_values_used'] == []
    detector = receipt['detector_controls']
    assert [detector[k] for k in ('double_span', 'binary_span', 'zero_span',
                                 'uniform_two_span')] == [[2], [1], [0], [0]]
    assert [r['body_spans_from_records'] for r in data['detector_trials']] == [
        [2], [1], [0], [0]]
    assert detector['forced_timeout_outcome'] == 'inconclusive'
    assert detector['forced_timeout_decoded'] is None
    assert detector['changed_current_rejected_by_independent_event_audit'] is True
    result['protocol_and_summary_verified'] = True
    return result


if __name__ == '__main__':
    print(json.dumps(check(Path(sys.argv[1]) if len(sys.argv) > 1 else HERE), indent=2))
