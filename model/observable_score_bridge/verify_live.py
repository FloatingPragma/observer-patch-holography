"""Verify externally applied dashboard statuses without modifying the website."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json

from .bridge import HERE, IDS, ROOT, authenticate, build, require


def verify_snapshot(live, proposed):
    baseline = json.loads((HERE / 'dashboard_inputs.json').read_text())
    old = {row['id']: row for row in baseline['observables']}
    rows = live['observables']
    current = {row['id']: row for row in rows}
    require(len(rows) == len(current) == 235, 'live catalogue has duplicate or changed membership')
    require(current.keys() == old.keys(), 'live observable membership changed')
    require(all(row['kind'] != 'input' and row['category'] != 'Foundations' for row in rows),
            'live observable filter changed')
    changed = {key for key in old if old[key]['status'] != current[key]['status']}
    require(changed == IDS, 'status changes differ from the eight reviewed recommendations')
    require(all(old[key]['status'] == 'computation' and current[key]['status'] == 'exectheorem'
                for key in IDS), 'an applied status differs from the recommendation')
    counts = Counter(row['status'] for row in rows)
    require(counts['theorem'] == 38 and counts['exectheorem'] == 20,
            'live completed counts changed')
    require(counts['conditional'] == 75 and counts['computation'] == 9,
            'live conditional/computation counts changed')
    selected_rows = live['selected_source_rows']
    selected = {row['id']: row for row in selected_rows}
    require(len(selected_rows) == len(selected) == 8 and selected.keys() == IDS,
            'live selected rows differ from the reviewed targets')
    require(all(selected[key]['status'] == current[key]['status'] for key in IDS),
            'live detail and summary statuses disagree')
    proposals = {row['id']: row for row in proposed['updates']}

    # The captured public text has not applied the two required semantic corrections.
    reference = selected['exec-coherent-screen']
    require('untouched reference alone unchanged' in reference['name'],
            'live reference title changed; review the wording again')
    require('coherence of an untouched reference qubit' in reference['why']
            and 'different reference coherences' in reference['summary'],
            'live reference explanation changed; review the wording again')
    integer = selected['exec-integer-lift']
    require(integer['value'] == 'exact identities'
            and 'Diagnostic only (not promoted)' in integer['summary'],
            'live integer title/value correction changed; review again')
    require('already carry a Standard Model hypercharge' in integer['why']
            and '68,515,840' in integer['summary'],
            'live integer explanation changed; review the wording again')
    corrections = []
    for key, reason in [
        ('exec-coherent-screen',
         'The title is corrected, but summary and explanation attribute a change to the untouched reference alone.'),
        ('exec-integer-lift',
         'The value and diagnostic label are corrected, but the explanation asserts physical hypercharge and the summary does not distinguish imported large-run counts from the executed finite reference.')]:
        row = selected[key]
        corrections.append(dict(
            id=key, status='exectheorem', status_change_required=False,
            current_name=row['name'], current_value=row['value'],
            current_summary=row['summary'], current_why=row['why'],
            reason=reason,
            replacement_summary=proposals[key]['summary'],
            replacement_why=proposals[key]['why'],
            evidence=proposals[key]['evidence']))
    return dict(
        schema='oph.observable_score_bridge.live_verification.v1',
        source_url=live['source_url'], source_javascript_asset=live['source_javascript_asset'],
        source_asset_sha256=live['source_asset_sha256'],
        live_completed=58, live_total=235, live_completed_percent=100 * 58 / 235,
        live_status_counts=dict(counts), all_eight_status_updates_applied=True,
        applied_ids=sorted(IDS), status_changes_needed=0,
        all_required_wording_corrections_applied=False,
        remaining_text_corrections=corrections,
        score_increase_from_remaining_text_corrections=0,
        website_written_by_this_package=False,
        identity_of_external_site_editor='not inferred',
        physical_completion_gates_accepted=[], new_empirical_observations=0)


def certificate():
    live = json.loads((HERE / 'live_status.json').read_text())
    proposed = authenticate(HERE / 'results.json')
    require(proposed == build(), 'original recommendation receipt is stale')
    result = verify_snapshot(live, proposed)
    paths = [HERE / name for name in ['live_status.json', 'verify_live.py',
             'LIVE_VERIFICATION.md', 'test_live.py', 'results.json']]
    result['source_bindings'] = {
        str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in paths}
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    result = certificate()
    data = json.dumps(result, indent=2, ensure_ascii=False, sort_keys=True) + '\n'
    path = HERE / 'live_verification.json'
    if args.write:
        path.write_text(data)
    if args.check:
        require(path.read_text() == data, 'live verification receipt is stale')
    print(json.dumps({key: result[key] for key in ['live_completed', 'live_total',
        'all_eight_status_updates_applied', 'status_changes_needed',
        'all_required_wording_corrections_applied']}))


if __name__ == '__main__':
    main()
