"""Freeze the bounded same-W instrument without awarding platform status."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PROVIDER = HERE.parent / 'native_walk_bridge'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run():
    names = ('NATIVE_PHASE_PROBE.json', 'SAME_WALK_RECORDS.json', 'INSTRUMENT_CONTRACT.json')
    receipts = {name: json.loads((HERE / name).read_text()) for name in names}
    checked = {}
    for name, receipt in receipts.items():
        assert receipt['checks'] and all(v is True for v in receipt['checks'].values())
        assert receipt['official_rows_changed'] == []
        for file, expected in receipt['source_sha256'].items():
            assert sha(ROOT / file) == expected, (name, file)
            checked[file] = expected
    expected_provider = {
        'common_walk.py': '4c8b6087aab6212ae12af126640ba19bbfec63215be5234be64f016cc604d75d',
        'compile_walk.py': '26371a86d04d7853ab03dbccb46f89a65c809fb6d75480312b0ed61b0259ada2',
        'COMMON_WALK_RECEIPT.json': '39e2d7a6dace9ff549d8e68974239c11f03316550a36f1c708302dd562048052',
        'COMMON_WALK_PROGRAM.json': '0439ac27c173d074871122cf893afd2617243003906f51865490cd6e4f73d7e7',
    }
    for file, expected in expected_provider.items():
        assert sha(PROVIDER / file) == expected, file
        checked[str((PROVIDER / file).relative_to(ROOT))] = expected
    owner = json.loads((PROVIDER / 'COMMON_WALK_RECEIPT.json').read_text())
    assert owner['checks'] and all(v is True for v in owner['checks'].values())
    for file, expected in owner.get('source_sha256', {}).items():
        assert sha(ROOT / file) == expected, file
        checked[file] = expected
    paths = [HERE / name for name in names]
    paths += [HERE / name for name in ('native_phase_probe.py', 'same_walk_records.py',
                                      'instrument_contract.py', 'METHOD.md', 'freeze_manifest.py')]
    review = HERE.parent / 'common_walk_review' / 'KINEMATIC_REVIEW.json'
    assert review.exists()
    paths.append(review)
    for path in paths:
        checked[str(path.relative_to(ROOT))] = sha(path)
    return {
        'schema': 'oph.common_walk.bounded_instrument_manifest.v1',
        'model_id': 'native-local12-current-v0-proposed',
        'worker': 'DYN02',
        'state': 'frozen bounded instrumentation; not complete platform validation',
        'checks': {'three_receipts_pass': True, 'all_receipt_source_pins_current': True,
                   'exact_provider_final_freeze_matches': True, 'official_rows_unchanged': True},
        'scope': ['same-W exact complete phase records', 'two-time sharp position instruments with full backaction',
                  'native all-column primitive controls and C/E pointer instrument',
                  'outward controlled-source wrapper bound and added reset failures'],
        'remaining_dependencies': [
            'separate quantitative original common-W and preparation compiler bounds for a numerical end-to-end physical error budget',
            'DYN03 frozen current-fiber intertwiner and coherent preparation before claiming the augmented-current experiment'],
        'no_independent_H_or_calibration': True,
        'official_rows_changed': [],
        'source_sha256': dict(sorted(checked.items())),
    }


if __name__ == '__main__':
    result = run()
    path = HERE / 'MANIFEST.json'
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'checks': result['checks'], 'pinned_files': len(result['source_sha256']),
                      'manifest_sha256': sha(path)}, indent=2))
