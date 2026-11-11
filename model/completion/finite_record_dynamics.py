"""Check the distinct finite inputs to the current completion work.

This checks the declared interfaces and source-bound receipts. Mathematical
replays belong to the three owning packages; no joint physical source or
completion-gate acceptance is certified here.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PATH = Path('model/completion/finite_record_dynamics.json')
CAMPAIGN = 'OPH-FINITE-RECORD-DYNAMICS-v1'
COMPONENTS = {
    'lepton_complete_cycle': ('model/lepton_cycle/results.json',
                             'model/lepton_cycle/test_cycle.py', ['CP-02', 'CP-07', 'CP-11', 'CP-13']),
    'bounded_current_record_observer': ('model/repair_sector_bridge/record_observer_results.json',
                                       'model/repair_sector_bridge/test_record_observer.py', ['CP-02', 'CP-13']),
    'native_closed_cell_epochs': ('model/native_abundance_epochs/results.json',
                                 'model/native_abundance_epochs/test_epochs.py', ['CP-05', 'CP-13']),
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def validate(packet, root=ROOT):
    root = Path(root).resolve()
    require(packet.get('schema') == 'oph.finite_record_dynamics.v1'
            and packet.get('campaign_id') == CAMPAIGN
            and packet.get('scope') == 'FINITE_COMPONENT_FAMILY', 'wrong finite component contract')
    require(all(packet.get(key) is False for key in (
        'whole_program_frozen', 'integrated_gate_accepted', 'common_source_join_established'))
        and packet.get('canonical_premises_discharged') == []
        and packet.get('physical_comparison_count_added') == 0, 'finite input promoted')
    rows = packet.get('components', [])
    require([row.get('id') for row in rows] == list(COMPONENTS), 'finite input coverage changed')
    for row in rows:
        path, check, gates = COMPONENTS[row['id']]
        require(row.get('path') == path and row.get('independent_check') == check
                and row.get('supports') == gates, 'finite input owner changed')
        require(set(row.get('interface', {})) == {'state', 'algebra', 'operation', 'clock', 'refinement'}
                and all(row['interface'].values()) and row.get('remaining'), 'finite input interface missing')
        require(row.get('common_source_transfer_proved') is False, 'separate sources identified')
        require(hashlib.sha256((root / path).read_bytes()).hexdigest() == row.get('sha256'),
                f'finite receipt changed: {path}')
        require((root / check).is_file(), f'missing owning check: {check}')
        receipt = json.loads((root / path).read_text())
        bindings = receipt.get('source_bindings', receipt.get('input_sha256'))
        require(isinstance(bindings, dict) and bindings, f'missing owning source bindings: {path}')
        for name, digest in bindings.items():
            source = (root / name).resolve()
            require(source.is_relative_to(root) and source.is_file()
                    and hashlib.sha256(source.read_bytes()).hexdigest() == digest,
                    f'owning source changed: {name}')
    return {'finite_components': 3, 'common_source_join_established': False,
            'integrated_gate_accepted': False, 'physical_comparison_count_added': 0}


if __name__ == '__main__':
    print(json.dumps(validate(json.loads((ROOT / PATH).read_text())), indent=2))
