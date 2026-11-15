"""Bind the full static native-charge criterion to exact public evidence."""
from pathlib import Path
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / 'code/openwave_charge_quantization'
OUTPUT = ROOT / 'claims/openwave_native_charge.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate():
    public = json.loads((PACKAGE / 'PUBLIC_MANIFEST.json').read_text())
    for name, expected in public['files_sha256'].items():
        path = PACKAGE / name
        if path.is_symlink() or sha(path) != expected:
            raise ValueError('Changed public evidence: ' + name)
    validation = json.loads((PACKAGE / 'VALIDATION.json').read_text())
    replay = json.loads((PACKAGE / 'RESULT.json').read_text())
    criterion = json.loads((PACKAGE / 'workspace/openwave-tuning/row_completion/ROW_01/CRITERION.json').read_text())
    assert validation['row'] == 1 and validation['full_criterion_met'] is True
    assert validation['criterion_sha256'] == hashlib.sha256(criterion['criterion'].encode()).hexdigest()
    assert replay['candidate_id'] == validation['candidate_id']
    assert set(replay['executions']) == {'candidate', 'auditor_root', 'auditor_dyn17', 'codex'}
    assert all(x['all_source_and_regenerated_output_bytes_match'] and x['outside_workspace_read_rejected']
               for x in replay['executions'].values())
    reviews = validation['independent_reviews']
    assert len({r['reviewer'] for r in reviews}) >= 2
    for review in reviews:
        path = PACKAGE / 'workspace' / review['manifest']
        assert sha(path) == review['manifest_sha256']
        data = json.loads(path.read_text())
        assert data['producer_independent'] is True and data['full_criterion_met'] is True
        for key in ('files', 'source_sha256', 'files_sha256'):
            for name, expected in data.get(key, {}).items():
                source = path.parent / name if key == 'files_sha256' else PACKAGE / 'workspace' / name
                assert sha(source) == expected, name
    return {'schema': 'oph.openwave.native_charge.full_criterion.v1',
            'candidate_id': validation['candidate_id'], 'row': 1,
            'criterion': criterion['criterion'], 'full_criterion_met': True,
            'independent_full_reviews': reviews, 'scope': validation['scope'],
            'complete_internal_rows_on_this_candidate': [1],
            'official_OpenWave_awards_added': 0, 'official_status': 'untested',
            'transfer_to_canonical_CAR_or_coherent_Dirac_field_source': False,
            'all31_complete': False,
            'public_replay': 'python3 code/openwave_charge_quantization/replay.py',
            'source_sha256': {str(path.relative_to(ROOT)): sha(path) for path in
                [Path(__file__), PACKAGE / 'PUBLIC_MANIFEST.json', PACKAGE / 'SNAPSHOT_MANIFEST.json',
                 PACKAGE / 'VALIDATION.json', PACKAGE / 'RESULT.json']}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true')
    group.add_argument('--check', action='store_true')
    args = parser.parse_args()
    text = json.dumps(generate(), indent=2, sort_keys=True) + '\n'
    if args.write:
        OUTPUT.write_text(text)
    elif OUTPUT.read_text() != text:
        raise SystemExit('Stale native charge full-criterion ledger')
    print('PASS_FULL_STATIC_NATIVE_CHARGE_CLAIM')


if __name__ == '__main__':
    main()
