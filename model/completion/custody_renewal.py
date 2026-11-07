"""Check a dated custody renewal without accepting new mathematical results.

Historical review snapshots and current live source bindings have different
roles. The renewal conserves the adopted class clauses, finite classifications
and completion decisions. It does not re-execute or re-prove their science.
"""
from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
RECORD = 'model/completion/custody_renewal_20261107.json'
REVIEW = 'model/completion/program_class_v2_review.json'
CLASS = 'model/completion/program_class_v2.json'
MAP = 'model/completion/finite_attachment_map.json'
PROGRAM = 'claims/completion_program.json'
INTERFACES = 'model/completion/component_interfaces.json'
BASE = '98391cb29e965cc7338acf5852f6b95b15c27f3c'
PREVIOUS = '3cf3c14810671a16e61b7915f69079a1773a12a2'
CLASS_DIGEST = 'cdbf861edd24dd8139509f1dde3f9769ab67f4e0636cd0393a090ba044c16480'
CURRENT_SCOPE = (
    'Complete 84-premise ownership. The historical review covered 145 domains, '
    'including 140 realized scoped attachments. The current map contains 156 '
    'domains, including 151 realized classifications. The 2026-11-07 renewal '
    'checks custody and preserves these classifications; it does not add '
    'mathematical acceptance or discharge canonical physical premises.'
)
HISTORICAL_FIELDS = (
    (10, 'model/particle_identification/results.json'),
    (11, 'model/particle_identification/results.json'),
    (11, 'model/protected_particle_readout/results.json'),
    (12, 'model/resolved_particle_spectrum/results.json'),
)
REQUIRED_BINDINGS = {
    MAP, INTERFACES, 'model/native_program/results.json',
    'model/completion/spectrum_einstein_campaign.json',
    'model/matched_einstein_limit/results.json',
    'code/a5_closure/manifests/source_current_order_sensitive_inventory.json',
    'code/a5_closure/source_current_order_sensitive_inventory.py',
    'code/a5_closure/verify_source_current_order_sensitive_inventory.py',
    'model/completion/custody_renewal.py',
    'model/completion/test_custody_renewal.py',
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(value):
    return hashlib.sha256(json.dumps(
        value, sort_keys=True, separators=(',', ':')
    ).encode()).hexdigest()


def file_digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def without_pointer_hashes(value):
    if isinstance(value, dict):
        return {key: without_pointer_hashes(child) for key, child in value.items()
                if key != 'sha256'}
    if isinstance(value, list):
        return [without_pointer_hashes(child) for child in value]
    return value


def program_semantics(program, historical_scope):
    result = copy.deepcopy(program)
    stage = next(row for row in result['stages'] if row['id'] == 'CP-13')
    rows = [row for row in stage['accepted_evidence'] if row['path'] == MAP]
    require(len(rows) == 1, 'unique finite-map scope required')
    require(rows[0]['scope'] in (historical_scope, CURRENT_SCOPE),
            'unreviewed finite-map scope')
    rows[0]['scope'] = historical_scope
    return without_pointer_hashes(result)


def review_history(review):
    return {key: value for key, value in review.items()
            if key not in ('source_bindings', 'current_source_custody_renewal')}


def validate(packet=None, root=ROOT):
    root = Path(root).resolve()
    if packet is None:
        packet = strict_json.load(root / RECORD)
    require(packet.get('schema') == 'oph.current_custody_renewal.v1'
            and packet.get('date') == '2026-11-07', 'wrong dated custody renewal')
    require(packet.get('canonical_base') == BASE
            and packet.get('historical_snapshot_base') == PREVIOUS,
            'wrong immutable predecessor')
    require(packet.get('new_mathematical_acceptance') is False
            and packet.get('canonical_physical_premises_discharged') == []
            and packet.get('gate_statuses_changed') is False,
            'custody renewal cannot promote scientific acceptance')
    require(packet.get('historical_scope') == {'domains': 145, 'realized': 140}
            and packet.get('current_scope_text') == CURRENT_SCOPE,
            'historical and current review scopes must remain distinct')

    finite = strict_json.load(root / MAP)
    counts = {'domains': len(finite['attachments']), 'realized': sum(
        row['classification'] == 'realized_on_claimed_domain'
        for row in finite['attachments'])}
    require(packet.get('current_classification_counts') == counts
            and counts == {'domains': 156, 'realized': 151},
            'current attachment classification counts differ')
    conserved = packet['conserved']
    for key in ('canonical_ownership', 'attachments', 'issue_closure_review',
                'reviewed_sources'):
        require(digest(finite[key]) == conserved['finite_map'][key],
                f'finite-map scientific body changed: {key}')
    interfaces = strict_json.load(root / INTERFACES)
    require(digest({key: value for key, value in interfaces.items()
                    if key != 'source_bindings'}) == conserved['component_interfaces'],
            'component integration scope changed')

    klass = strict_json.load(root / CLASS)
    contract = digest({key: value for key, value in klass.items() if key != 'adoption'})
    require(contract == packet['reviewed_contract_digest'] == CLASS_DIGEST,
            'adopted scientific class clauses changed')
    program = strict_json.load(root / PROGRAM)
    require(digest(program_semantics(program, packet['historical_scope_text']))
            == conserved['completion_program'],
            'completion decisions or scientific contracts changed')
    stage = next(row for row in program['stages'] if row['id'] == 'CP-13')
    require(next(row for row in stage['accepted_evidence'] if row['path'] == MAP)
            ['scope'] == CURRENT_SCOPE, 'current finite-map scope is unexplained')

    review = strict_json.load(root / REVIEW)
    repairs = packet['historical_snapshot_restorations']
    require([(row['refresh_index'], row['receipt']) for row in repairs]
            == list(HISTORICAL_FIELDS), 'historical restoration coverage changed')
    for row in repairs:
        refresh = review['review_refreshes'][row['refresh_index']]
        require(refresh['date'] == '2026-11-01'
                and refresh['source_receipt_snapshots'][row['receipt']]
                == row['restored_sha256'], 'historical review snapshot changed')
    require(digest(review_history(review)) == conserved['review_history'],
            'historical review content changed')

    require(set(packet['current_bindings']) == REQUIRED_BINDINGS,
            'current custody binding coverage changed')
    for relative, expected in packet['current_bindings'].items():
        path = (root / relative).resolve()
        require(path.is_relative_to(root) and path.is_file()
                and file_digest(path) == expected,
                f'current custody binding changed: {relative}')
    pointer = review.get('current_source_custody_renewal', {})
    require(pointer.get('date') == packet['date'] and pointer.get('path') == RECORD
            and pointer.get('role') == 'CUSTODY_ONLY_NO_NEW_MATHEMATICAL_ACCEPTANCE'
            and pointer.get('sha256') == file_digest(root / RECORD),
            'dated review does not bind current renewal')
    require(review['source_bindings'].get(RECORD) == file_digest(root / RECORD)
            and klass['adoption']['independent_review']['sha256']
            == file_digest(root / REVIEW), 'adoption custody does not reach renewal')
    return {'date': packet['date'], 'classification_counts': counts,
            'reviewed_contract_digest': contract,
            'scientific_status_changed': False,
            'new_mathematical_acceptance': False}


def main():
    print(json.dumps(validate(), indent=2))


if __name__ == '__main__':
    main()
