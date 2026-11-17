"""Bind the common-source charge criterion to its exact public evidence."""
from pathlib import Path
import argparse
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
PACKAGE = ROOT / 'code/openwave_common_charge'
OUTPUT = ROOT / 'claims/openwave_common_charge.json'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def generate():
    public = json.loads((PACKAGE / 'PUBLIC_MANIFEST.json').read_text())
    for name, expected in public['files_sha256'].items():
        path = PACKAGE / name
        if path.is_symlink() or sha(path) != expected:
            raise ValueError('Changed public charge evidence: ' + name)
    validation = json.loads((PACKAGE / 'VALIDATION.json').read_text())
    replay = json.loads((PACKAGE / 'RESULT.json').read_text())
    assert validation['row'] == 1 and validation['full_criterion_met'] is True
    assert validation['criterion_sha256'] == hashlib.sha256(validation['literal_criterion'].encode()).hexdigest()
    assert set(replay['executions']) == {'topology', 'native', 'terminal', 'terminal_opposing', 'terminal_channel'}
    for result in replay['executions'].values():
        assert result['all_source_and_regenerated_output_bytes_match'] is True or (
            result['all_source_bytes_match'] is True and
            result['scientific_replay_matches_with_declared_roundoff'] is True)
        assert result['runtime']['external_project_read_probe_rejected'] is True
        assert result['runtime']['all_project_modules_inside_snapshot'] is True
    for relative, expected in validation['evidence_files'].items():
        assert sha(PACKAGE / 'workspace' / relative) == expected, relative
    reviews = validation['independent_reviews']
    assert len({r['reviewer'] for r in reviews}) >= 2
    for review in reviews:
        assert review['producer_independent'] is True and review['full_criterion_met'] is True
        path = PACKAGE / 'workspace' / review['review_file']
        assert sha(path) == review['review_sha256']
        decision = json.loads(path.read_text())
        assert decision['full_criterion_met'] is True and decision['producer_independent'] is True
    assert validation['original_huge_field_prefix_executed'] is False
    assert validation['new_terminal_program_executed'] is True
    assert validation['global_torus_bundle_claimed'] is False
    assert validation['terminal_annex_independent_reviews'] == 2
    assert len(set(validation['terminal_annex_review_manifests'])) == 2
    central = validation['central_alternative_validation']
    assert sha(PACKAGE / 'workspace' / central['path']) == central['sha256']
    assert json.loads((PACKAGE / 'workspace' / central['path']).read_text())['terminal_annex']['twice_independently_reviewed'] is True
    scope = validation['authoritative_scope_clarification']
    scope_pins = json.loads((PACKAGE / scope['input_index']).read_text())
    for name, digest in scope_pins.items():
        assert sha(PACKAGE / 'scope_inputs' / digest) == digest, name
    assert sha(PACKAGE / scope['scope_file']) == scope['scope_sha256']
    current = json.loads((PACKAGE / 'scope_inputs' / scope['current_central_acceptance_sha256']).read_text())
    assert current['full_criterion_met'] is True
    assert current['additional_full_review']['full_criterion_met'] is True
    assert current['additional_full_review']['review_sha256'] == scope['additional_Fable_full_review_sha256']
    assert scope['scientific_executions_added'] == 0
    return {
        'schema': 'oph.openwave.common_charge.full_criterion.v1',
        'row': 1, 'candidate_id': validation['candidate_id'],
        'literal_criterion': validation['literal_criterion'],
        'full_criterion_met': True, 'scope': validation['scope'],
        'independent_full_reviews': reviews,
        'terminal_annex_independent_reviews': validation['terminal_annex_independent_reviews'],
        'original_huge_field_prefix_executed': False,
        'terminal_observable_program_executed': True,
        'source_occupation_integer_ancestry_declared': True,
        'authoritative_scope_clarification': scope,
        'global_torus_bundle_claimed': False,
        'official_status': 'untested', 'official_OpenWave_awards_added': 0,
        'all31_complete': False, 'one_interacting_continuum_construction_established': False,
        'public_replay': 'python3 code/openwave_common_charge/replay.py --mode all',
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
        raise SystemExit('Stale common-source charge claim ledger')
    print('PASS_COMMON_SOURCE_CHARGE_CLAIM')


if __name__ == '__main__':
    main()
