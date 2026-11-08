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
BASE = '6b20cfa05ad82582bb1c8a786b32b7556b4008d6'
PREVIOUS = '3cf3c14810671a16e61b7915f69079a1773a12a2'
CLASS_DIGEST = '3fcac70fab6e158362a2e00c4dc768251bf763184c9a517adbfc4bd9593101c5'
INTEGRATION_CLASS_DIGEST = '54fe2116bef3a7cb2eab1678ff20281d80f01032d21d5024009424fe5ff8d877'
INTEGRATION_BASELINE = {'commit': '85a0cc208f52b385cfa674b288bcb079156abfec',
 'source_sha256': {'claims/completion_program.json': '1847a49457e31b84894fda42dd59c1ec04fc7a5fe800bc52320b6222ea0d7338',
                   'model/completion/component_interfaces.json': '1ddf4f6ff53ea5541828e05f5132594bca4c587bd5b42ca2f61f71ee57924a5d',
                   'model/completion/finite_attachment_map.json': 'e95ae10bf0b4a486f19f1b0f6af0c38720d8e45af834f21bd871d18000fecfac',
                   'model/completion/program_class_v2.json': '02fadb6bceabe3876776d2fce9da76707df5b95598595fcdda9d7c69b47c0522',
                   'model/completion/program_class_v2_review.json': 'df9401a150f759927da5054d16f300df07f664c1ea0de2a90f0810415952cfd7'}}
ORIGINAL_CLASS_DIGEST = 'ea7a6107552b3abbc141a9f774fd1301762853c96869a7a09db137bd296c5954'
ORIGINAL_CONSERVED = {'completion_program': 'f0945c71f28e476e946506cc0cd5c03f0aead8c4624828bb239980132da81e22',
 'component_interfaces': '814aa8ca294343f3d9167cec0e1f82587517d363c15a846c2cbe8ef55b594bb0',
 'finite_map': '36f79f223fe732f367869626c562a861241086d875763dbaef1bcae37f1f4d4a',
 'review_history': 'feca1d0a11e2ee50be6765a63b75de5f4a8c95a2f889547541f9f6949ea4f384'}
HISTORICAL_SCOPE = ('Complete 84-premise ownership and 145 reviewed domains, '
                    'including 140 realized scoped attachments. Canonical physical '
                    'discharge and gate statuses are unchanged.')
# Derived from the immutable integration tree, with only the four predecessor
# restorations and the named live-pointer/scope normalization; never from the renewal.
EXPECTED_CONSERVED = {'completion_program': '8a41a23ae63f97ac68c65d26fb4d8dc876dc144563e49ea909fa36140a73bacc',
 'component_interfaces': '814aa8ca294343f3d9167cec0e1f82587517d363c15a846c2cbe8ef55b594bb0',
 'finite_map': '3dd36dc4ffb4dc490729740f6e5b260026ae86359b30b534edcf85dfaa8938df',
 'review_history': '0d46ba4927a7811f2da209acdc54c9aa99578d084ce3fed9c656cf064994825d'}
ADOPTION_DIGEST = '03f52dbb33b86b8588f3172b23d2bf4acd01a59af73bf2363ce127bd1e8f44ed'
BINDING_PATH_DIGESTS = {'model/completion/component_interfaces.json': '1f75d695f9419b42e21d83f586bcc5ef1397fe8b1c11ebc33764790b794d79f1',
 'model/completion/finite_attachment_map.json': 'c8ddcc7952b113f5739c715d065a81ff71670b384dab5b5699b61d3bc9fdd29a',
 'model/completion/program_class_v2_review.json': 'ce6d0c427a05e42f5d45d98888b9b646635b81ef56b3d534407ce15f0c5c3060'}
PREDECESSOR_REVIEW_SHA256 = 'c3afbb4442df07154d894707f43c171536117557c17e9a8329f0c52f1556e9d0'
HISTORICAL_REVIEW_SHA256 = 'ad7a397e42666ecfd6a7be984f0b9443b4314a6aaf01fdf351a64f03c11471bf'
CURRENT_SCOPE = ('Complete 84-premise ownership and 156 classified domains, including 151 realized scoped attachments. '
 'These classifications do not discharge canonical physical premises or change gate statuses.')
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
PROGRAM_POINTER_SCHEMA_SHA256 = '978680d08bfdf6c35fedacd5c26043436d27a6865d071a8252c422676a80e433'
PROGRAM_POINTER_COUNT = 736
PROGRAM_TARGET_COUNT = 267
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

# Immutable accepted main supplies the authority for this live-source refresh.
SOURCE_REFRESH_BASELINE = {'commit': '12047f556a3f4d6bae5b48aa5b85980a0c5eda40',
 'source_sha256': {'claims/completion_program.json': '6ded0a1096ec34270311e3069fec700e412cc2c0e31d4ab3fad4896b7a1794cb',
                   'model/completion/component_interfaces.json': '220b389f0b6f956969d0226c2ac121bb6728ec7b979bdf29319c59771279ac72',
                   'model/completion/custody_renewal.py': '2ef536f614e7ef038447a997903f3048bb936313e0fcfad9dc41f6f3d28e3194',
                   'model/completion/custody_renewal_20261107.json': 'ac4665183969ad1f414db6fc797f1669b99aba7a0f8f947d454bd28504c6c71c',
                   'model/completion/finite_attachment_map.json': 'e95ae10bf0b4a486f19f1b0f6af0c38720d8e45af834f21bd871d18000fecfac',
                   'model/completion/program_class_v2.json': '464e85412867bb7a4cdccf42d7a291f8da72daa875b11c02d051d3ffa557af63',
                   'model/completion/program_class_v2_review.json': '7c3dc4b3493a3ea430d4805576051e8b115f0ca83290ebac8133c439711767dc',
                   'model/completion/test_custody_renewal.py': '9554aa1e5a633beb5a600b92822d5da07a0652a7567fe451e77dca53c353aed7'}}
CLASS_LIVE_RECEIPT = 'model/native_program/results.json'
CLASS_SEMANTICS_SHA256 = '3a01c7d2b459de4a43807e6e7f1d0da589bfc81e0e0f9d5684940887cb6ba4b5'
BASELINE_NAVIGATION = {'claim_count': 389,
 'contract_coverage_sha256': 'ddcbfdb1a9dbfcbbf708eab51c0225d18b8135153e054b65b5a0efda31d72dfe',
 'lean_module_count': 676,
 'note': 'Claim and Lean navigation includes native record probability and thermodynamics, '
         'count-derived modular histories, complete transfer localization and spectrum, '
         'full-sector nonquadratic interaction, and regional entropy response. Canonical premise '
         'contracts, physical attachment edges, discharge statuses and finite attachment '
         'classifications are unchanged.'}
NAVIGATION_SOURCE = 'tracking/premise_register.json'
CURRENT_CLAIM_COUNT = 392
CURRENT_COVERAGE_SHA256 = '99bfa21a39222d236dc5d03226f44bd25a74f2d1f7627fe0f082acd088ef349a'


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


def program_pointers(program):
    """Authenticate every live evidence location against the immutable baseline."""
    rows = []
    def walk(value, location=()):
        if isinstance(value, dict):
            if 'path' in value and 'sha256' in value:
                require(type(value['path']) is str, 'program pointer path type changed')
                rows.append((location, value['path']))
            for key in sorted(value):
                walk(value[key], location + (key,))
        elif isinstance(value, list):
            for index, child in enumerate(value):
                walk(child, location + (index,))
    walk(program)
    require(len(rows) == PROGRAM_POINTER_COUNT
            and len({target for _, target in rows}) == PROGRAM_TARGET_COUNT
            and digest(rows) == PROGRAM_POINTER_SCHEMA_SHA256,
            'immutable program pointer schema changed')
    return tuple(rows)


def program_semantics(program, historical_scope):
    require(historical_scope == HISTORICAL_SCOPE, 'historical scope authority changed')
    result = copy.deepcopy(program)
    require([result['stages'][index]['id'] for index in range(3)]
            == ['CP-00', 'CP-01', 'CP-13'], 'program stage order changed')
    row = at_path(result, ('stages', 2, 'accepted_evidence', 0))
    require(row['path'] == MAP and row['scope'] in (HISTORICAL_SCOPE, CURRENT_SCOPE),
            'unreviewed finite-map scope')
    row['scope'] = HISTORICAL_SCOPE
    for path, target in program_pointers(program):
        row = at_path(result, path)
        require(row.get('path') == target and 'sha256' in row,
                'program custody pointer location changed')
        del row['sha256']
    return result


def review_history(review):
    result = {key: copy.deepcopy(value) for key, value in review.items()
              if key not in ('source_bindings', 'current_source_custody_renewal')}
    require(result.get('reviewed_contract_digest') == CLASS_DIGEST,
            'current review class digest changed')
    result['reviewed_contract_digest'] = INTEGRATION_CLASS_DIGEST
    return result


def class_semantics(klass, root):
    """Only the named live member receipt is excluded from class clauses."""
    result = {key: copy.deepcopy(value) for key, value in klass.items() if key != 'adoption'}
    bindings = result['interpreter_binding']['source_hashes']
    require(bindings.get(CLASS_LIVE_RECEIPT) == file_digest(Path(root) / CLASS_LIVE_RECEIPT),
            'class live member receipt differs from its owning bytes')
    bindings[CLASS_LIVE_RECEIPT] = '<live source SHA256>'
    return result


def finite_semantics(finite, root):
    """Keep all classifications; renew only two authenticated navigation fields."""
    result = copy.deepcopy(without_source_bindings(finite))
    coverage = strict_json.load(Path(root) / NAVIGATION_SOURCE)['model_supply']['contract_coverage']
    require(len(coverage['claim_ids']) == CURRENT_CLAIM_COUNT
            and len(coverage['lean_files']) == BASELINE_NAVIGATION['lean_module_count']
            and digest(coverage) == CURRENT_COVERAGE_SHA256,
            'immutable current navigation authority changed')
    navigation = result['existing_navigation_coverage']
    require(type(navigation.get('claim_count')) is int
            and navigation['claim_count'] == CURRENT_CLAIM_COUNT
            and navigation.get('contract_coverage_sha256') == CURRENT_COVERAGE_SHA256,
            'current finite-map navigation changed')
    navigation['claim_count'] = BASELINE_NAVIGATION['claim_count']
    navigation['contract_coverage_sha256'] = BASELINE_NAVIGATION['contract_coverage_sha256']
    return result


def validate(packet=None, root=ROOT):
    root = Path(root).resolve()
    if packet is None:
        packet = strict_json.load(root / RECORD)
    require(packet.get('schema') == 'oph.current_custody_renewal.v1'
            and packet.get('date') == '2026-11-07', 'wrong dated custody renewal')
    require(packet.get('integration_baseline') == INTEGRATION_BASELINE,
            'immutable integration baseline changed')
    require(packet.get('source_refresh_baseline') == SOURCE_REFRESH_BASELINE,
            'immutable accepted-main refresh authority changed')
    require(packet.get('original_conserved') == ORIGINAL_CONSERVED
            and packet.get('original_reviewed_contract_digest') == ORIGINAL_CLASS_DIGEST,
            'immutable original conservation authority changed')
    require(packet.get('program_pointer_schema') == {
        'locations': PROGRAM_POINTER_COUNT, 'targets': PROGRAM_TARGET_COUNT,
        'sha256': PROGRAM_POINTER_SCHEMA_SHA256},
        'immutable program pointer authority changed')
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
    require(digest(finite_semantics(finite, root)) == EXPECTED_CONSERVED['finite_map'],
            'finite-map scientific body changed')
    interfaces = strict_json.load(root / INTERFACES)
    require(digest(without_source_bindings(interfaces)) == EXPECTED_CONSERVED['component_interfaces'],
            'component integration scope changed')

    klass = strict_json.load(root / CLASS)
    contract = digest({key: value for key, value in klass.items() if key != 'adoption'})
    require(contract == packet['reviewed_contract_digest'] == CLASS_DIGEST,
            'adopted scientific class clauses changed')
    require(digest(class_semantics(klass, root)) == CLASS_SEMANTICS_SHA256,
            'immutable accepted class clauses changed outside the named live source pin')
    adoption = copy.deepcopy(klass['adoption'])
    adoption['independent_review'].pop('sha256')
    require(adoption['independent_review'].get('reviewed_contract_digest') == CLASS_DIGEST,
            'adoption does not bind the current exact class digest')
    adoption['independent_review']['reviewed_contract_digest'] = INTEGRATION_CLASS_DIGEST
    require(digest(adoption) == ADOPTION_DIGEST, 'adoption state or meaning changed')
    program = strict_json.load(root / PROGRAM)
    require(digest(program_semantics(program, packet['historical_scope_text']))
            == conserved['completion_program'],
            'completion decisions or scientific contracts changed')
    for path, target in program_pointers(program):
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
