import json
from pathlib import Path

from . import consensus


def test_native_fork_and_arbitration_match_local_records():
    data = consensus.run()
    assert data == json.loads(Path(consensus.__file__).with_suffix('.json').read_text())
    cases = {(r['double_at'], r['priority']): r for r in data['rows']}
    fork = cases[0, None]
    assert [p['terminal']['pair_flags'] for p in fork['schedules']] == [[1, 0], [0, 1]]
    assert all(p['terminal_retries'] == [0, 0] for p in fork['schedules'])
    for side in (0, 1):
        row = cases[0, side]
        assert row['terminal_state_fidelity'] == 1
        for p in row['schedules']:
            assert p['terminal']['pair_flags'][side] == 1
    assert len(cases) == 9
