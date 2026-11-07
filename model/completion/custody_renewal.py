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
BASE = '3274bf2f995ff2631a7ca56e448e222b95fc4fe4'
PREVIOUS = '3cf3c14810671a16e61b7915f69079a1773a12a2'
CLASS_DIGEST = '74ceb649290a345795adb92e6e796d8f81429d5620739a396c992f88074697ea'
HISTORICAL_SCOPE = ('Complete 84-premise ownership and 145 reviewed domains, '
                    'including 140 realized scoped attachments. Canonical physical '
                    'discharge and gate statuses are unchanged.')
# Derived from the immutable BASE/PREVIOUS Git objects, never from the renewal.
EXPECTED_CONSERVED = {'completion_program': 'be6707e89b13cd67164e56794e5399b5e2be183b0b6478b9d7d0d6bd4627514f',
 'component_interfaces': '814aa8ca294343f3d9167cec0e1f82587517d363c15a846c2cbe8ef55b594bb0',
 'finite_map': '1a955504363ec1690dcad5c0840dd318500a6b8af28e5117a6c7d542188b1594',
 'review_history': '5c84d0a9cf8e93d4c1baa9ab4c555e1582b3055ff918040e122b0e3c38fedbd6'}
ADOPTION_DIGEST = '5e935518c7e99b3950cd7b5310b4aaf5a193fa7082d555a5c8b184a3dee61493'
BINDING_PATH_DIGESTS = {'model/completion/component_interfaces.json': '1f75d695f9419b42e21d83f586bcc5ef1397fe8b1c11ebc33764790b794d79f1',
 'model/completion/finite_attachment_map.json': 'c8ddcc7952b113f5739c715d065a81ff71670b384dab5b5699b61d3bc9fdd29a',
 'model/completion/program_class_v2_review.json': 'ce6d0c427a05e42f5d45d98888b9b646635b81ef56b3d534407ce15f0c5c3060'}
PREDECESSOR_REVIEW_SHA256 = 'de8088697bb32703fdbb3b1013c86a42d406d26e8f8e58916949b5a585a70291'
HISTORICAL_REVIEW_SHA256 = 'ad7a397e42666ecfd6a7be984f0b9443b4314a6aaf01fdf351a64f03c11471bf'
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
EXPECTED_RESTORATIONS = (('71dd82b993200f13b5bb3b9924021231074dc97da58f1d76881dd67be8c20c64',
  '81f4a4482812bab17f96c2c9e93146d12c02a7177a02447832b78b2b22575864'),
 ('71dd82b993200f13b5bb3b9924021231074dc97da58f1d76881dd67be8c20c64',
  '812b5cd7c70dcca7c3a585fa86858d645a8bd6f734ee41acdccd5a8eaad51ac5'),
 ('61d0b99e000dfd336f0b3264c69c8552e3096ced154cf50cad22dc37603ecff1',
  'd89f9bef2a6f059bf2e4d738bc390e7ab050e29eb29873764c61c6d48eab06de'),
 ('d46c07cae1ddb684d0de31fefe334f31cdbe2bb8d99b037deaa2eb00e6fae259',
  'dadd4fd2b5958478a7adb8b70015222848f6dad730588dacd71184aca22f5fa6'))
PROGRAM_POINTERS = ((('stages', 0, 'accepted_evidence', 0), 'model/completion/program_class_v2.json'),
 (('stages', 0, 'accepted_evidence', 1), 'model/completion/program_class_v2_review.json'),
 (('stages', 1, 'accepted_evidence', 1), 'model/completion/program_class_v2_review.json'),
 (('stages', 2, 'accepted_evidence', 0), 'model/completion/finite_attachment_map.json'),
 (('stages', 2, 'accepted_evidence', 1), 'model/completion/program_class_v2_review.json'),
 (('stages', 2, 'progress_evidence', 22), 'model/completion/finite_attachment_map.json'),
 (('component_campaigns', 21, 'contract'), 'model/completion/component_interfaces.json'),
 (('next_wave', 'evidence', 49), 'model/completion/finite_attachment_map.json'),
 (('adopted_class',), 'model/completion/program_class_v2.json'),
 (('acceptance_review',), 'model/completion/program_class_v2_review.json'))
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


def at_path(value, path):
    for key in path:
        value = value[key]
    return value


def without_source_bindings(packet):
    return {key: value for key, value in packet.items() if key != 'source_bindings'}


def program_semantics(program, historical_scope):
    require(historical_scope == HISTORICAL_SCOPE, 'historical scope authority changed')
    result = copy.deepcopy(program)
    require([result['stages'][index]['id'] for index in range(3)]
            == ['CP-00', 'CP-01', 'CP-13'], 'program stage order changed')
    row = at_path(result, ('stages', 2, 'accepted_evidence', 0))
    require(row['path'] == MAP and row['scope'] in (HISTORICAL_SCOPE, CURRENT_SCOPE),
            'unreviewed finite-map scope')
    row['scope'] = HISTORICAL_SCOPE
    for path, target in PROGRAM_POINTERS:
        row = at_path(result, path)
        require(row.get('path') == target and 'sha256' in row,
                'program custody pointer location changed')
        del row['sha256']
    return result


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
    require(packet.get('predecessor_review_sha256') == PREDECESSOR_REVIEW_SHA256
            and packet.get('historical_review_sha256') == HISTORICAL_REVIEW_SHA256,
            'immutable review authority changed')
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
    require(conserved == EXPECTED_CONSERVED, 'immutable conservation authority changed')
    require(digest(without_source_bindings(finite)) == EXPECTED_CONSERVED['finite_map'],
            'finite-map scientific body changed')
    interfaces = strict_json.load(root / INTERFACES)
    require(digest(without_source_bindings(interfaces)) == EXPECTED_CONSERVED['component_interfaces'],
            'component integration scope changed')

    klass = strict_json.load(root / CLASS)
    contract = digest({key: value for key, value in klass.items() if key != 'adoption'})
    require(contract == packet['reviewed_contract_digest'] == CLASS_DIGEST,
            'adopted scientific class clauses changed')
    adoption = copy.deepcopy(klass['adoption'])
    adoption['independent_review'].pop('sha256')
    require(digest(adoption) == ADOPTION_DIGEST, 'adoption state or meaning changed')
    program = strict_json.load(root / PROGRAM)
    require(digest(program_semantics(program, packet['historical_scope_text']))
            == conserved['completion_program'],
            'completion decisions or scientific contracts changed')
    for path, target in PROGRAM_POINTERS:
        require(at_path(program, path)['sha256'] == file_digest(root / target),
                'program custody target changed: ' + target)
    stage = next(row for row in program['stages'] if row['id'] == 'CP-13')
    require(next(row for row in stage['accepted_evidence'] if row['path'] == MAP)
            ['scope'] == CURRENT_SCOPE, 'current finite-map scope is unexplained')

    review = strict_json.load(root / REVIEW)
    repairs = packet['historical_snapshot_restorations']
    require([(row['refresh_index'], row['receipt']) for row in repairs]
            == list(HISTORICAL_FIELDS), 'historical restoration coverage changed')
    for row, (before, original) in zip(repairs, EXPECTED_RESTORATIONS, strict=True):
        require(row.get('pre_renewal_sha256') == before
                and row.get('restored_sha256') == original,
                'immutable historical restoration changed')
        refresh = review['review_refreshes'][row['refresh_index']]
        require(refresh['date'] == '2026-11-01'
                and refresh['source_receipt_snapshots'][row['receipt']]
                == row['restored_sha256'], 'historical review snapshot changed')
    require(digest(review_history(review)) == conserved['review_history'],
            'historical review content changed')
    for name, data in ((MAP, finite), (INTERFACES, interfaces), (REVIEW, review)):
        require(digest(sorted(data['source_bindings'])) == BINDING_PATH_DIGESTS[name],
                'source binding path coverage changed: ' + name)

    require(set(packet['current_bindings']) == REQUIRED_BINDINGS,
            'current custody binding coverage changed')
    for relative, expected in packet['current_bindings'].items():
        path = (root / relative).resolve()
        require(path.is_relative_to(root) and path.is_file()
                and file_digest(path) == expected,
                f'current custody binding changed: {relative}')
    pointer = review.get('current_source_custody_renewal', {})
    require(set(pointer) == {'date', 'path', 'role', 'sha256'}
            and pointer.get('date') == packet['date'] and pointer.get('path') == RECORD
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
