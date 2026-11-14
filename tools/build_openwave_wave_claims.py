"""Bind two source-specific wave criteria to complete public review evidence."""
from pathlib import Path
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / 'code/openwave_wave_rows'
OUTPUT = ROOT / 'claims/openwave_wave_rows.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate():
    public = json.loads((PACKAGE / 'PUBLIC_MANIFEST.json').read_text())
    for name, expected in public['files_sha256'].items():
        path = PACKAGE / name
        if path.is_symlink() or sha(path) != expected:
            raise ValueError('Changed public wave evidence: ' + name)
    validation = json.loads((PACKAGE / 'VALIDATION.json').read_text())
    replay = json.loads((PACKAGE / 'RESULT.json').read_text())
    assert validation['same_dynamics_claimed'] is False
    assert set(replay['executions']) == {'maxwell', 'size3', 'kg_root', 'kg_opposing'}
    for result in replay['executions'].values():
        assert result['all_source_and_regenerated_output_bytes_match'] is True
        assert result['runtime']['external_project_read_probe_rejected'] is True
        assert result['runtime']['all_project_modules_inside_snapshot'] is True
    rows = []
    for row in validation['rows']:
        assert row['full_criterion_met'] is True and row['candidate_id']
        assert row['criterion_sha256'] == hashlib.sha256(row['literal_criterion'].encode()).hexdigest()
        for relative, expected in row['evidence_files'].items():
            assert sha(PACKAGE / 'workspace' / relative) == expected, relative
        reviews = row['independent_reviews']
        assert len({r['reviewer'] for r in reviews}) >= 2
        for review in reviews:
            assert review['producer_independent'] is True and review['full_criterion_met'] is True
            path = PACKAGE / 'workspace' / review['review_file']
            assert sha(path) == review['review_sha256']
            decision = json.loads(path.read_text())
            passed = decision.get('full_criterion_met') is True or decision.get('decision') == 'PASS_FULL_LITERAL_CRITERION_ON_DECLARED_NATIVE_ARCHITECTURE'
            assert passed, review['review_file']
        rows.append({'row': row['row'], 'candidate_id': row['candidate_id'],
                     'literal_criterion': row['literal_criterion'],
                     'full_criterion_met': True, 'scope': row['scope'],
                     'independent_full_reviews': reviews, 'official_status': 'untested'})
    assert {r['row'] for r in rows} == {29, 30}
    return {'schema': 'oph.openwave.wave_rows.full_criteria.v1', 'rows': rows,
            'same_source_transfer_between_rows': False,
            'official_OpenWave_awards_added': 0, 'all31_complete': False,
            'one_interacting_continuum_construction_established': False,
            'public_replay': 'python3 code/openwave_wave_rows/replay.py --mode all',
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
        raise SystemExit('Stale source-specific wave claim ledger')
    print('PASS_SOURCE_SPECIFIC_MAXWELL_AND_KLEIN_GORDON_CLAIMS')


if __name__ == '__main__':
    main()
