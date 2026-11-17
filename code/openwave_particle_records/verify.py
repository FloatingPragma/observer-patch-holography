"""Verify archived particle-record evidence and its bounded claim projection.

This is a custody and receipt-consistency check, not a scientific rerun.
Use --workspace to check every leaf of the original source manifests too.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text())


def leaves(path):
    if path.suffix == '.json':
        return {name: value['sha256'] if isinstance(value, dict) else value
                for name, value in read(path)['files'].items()}
    return dict((line.split(maxsplit=1)[1].lstrip('*'), line.split()[0])
                for line in path.read_text().splitlines() if line.strip())


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify(workspace=None):
    inventory = read(HERE / 'EVIDENCE.json')
    count = 0
    original_count = 0
    for provider in inventory['providers']:
        base = HERE / provider['snapshot']
        manifest = base / provider['manifest_name']
        require(digest(manifest) == provider['manifest_sha256'], str(manifest))
        entries = leaves(manifest)
        require(len(entries) == provider['sealed_leaf_count'], 'leaf count')
        for name, expected in provider['archived_files'].items():
            require(entries.get(name) == expected, 'archive not covered: ' + name)
            require(digest(base / name) == expected, 'archive mismatch: ' + name)
            count += 1
        if workspace is not None:
            original = Path(workspace) / provider['workspace_path']
            require(digest(original / provider['manifest_name']) == provider['manifest_sha256'],
                    'original manifest: ' + str(original))
            for name, expected in entries.items():
                require(digest(original / name) == expected, 'original leaf: ' + str(original / name))
                original_count += 1

    def receipt(name):
        return read(HERE / 'evidence' / name / 'RECEIPT.json')

    static = receipt('static_candidate_join')
    require(static['all_pass'] and static['native_calls_total'] == 435080, 'static calls')
    require(static['official_rows_solved'] == [] and not static['G6_confirmed'], 'static scope')
    compact = receipt('compact_static_adjudication')
    require(compact['pair_all_columns_checked'] == 82575360, 'compact columns')
    require(compact['native_calls_total'] == 12798, 'columns are not calls')
    clock = receipt('compact_internal_clock')['scientific']
    require(clock['localization']['exact_integer_invariant_span_checked'], 'localization')
    require(clock['localization']['all_integer_times_proved'], 'all integer times')
    require(not any(clock['energy']['mass_clock_gates'].values()), 'clock failures')
    require(not clock['onset']['autonomous_onset_established'], 'onset scope')
    require(clock['source']['native_calls'] == 303660, 'arc calls')
    spin = receipt('compact_field_spin_join__existing_audit')
    require(spin['C_complete_instrument_intertwiner'] and not spin['E_coarse_record_intertwiner'],
            'C/E source separation')
    require(spin['L6_tensor_Spin2_fundamental_multiplicity'] == 0, 'spatial-spin control')
    flavour = receipt('neutral_flavour_join')
    require(not flavour['full_PMNS_criterion_satisfied'] and not flavour['row_complete'], 'PMNS scope')
    require(flavour['primary_flights'] == 840 and flavour['sampled_primary_timeouts'] == 39,
            'retained timeouts')
    boundary = receipt('neutral_source_completion')
    require(boundary['literal_native_calls'] == 6400 and boundary['full_rows_solved'] == [], 'boundary scope')
    require(all(v == 0 for v in boundary['seeded_accepts_by_trigger'].values()), 'short accepts')
    assessment = read(HERE / 'ASSESSMENT.json')
    require([row['row'] for row in assessment['rows']] == list(range(1, 11)), 'row coverage')
    ids = {p['id'] for p in inventory['providers']}
    for row in assessment['rows']:
        require(not row['full_row_solved'] and row['official'] == 'untested', 'row promotion')
        require(set(row['evidence']).issubset(ids), 'missing row provider')
    require(all(not p['completed_evidence'] for p in assessment['unexecuted_proposals']), 'proposal promotion')
    return {'providers': len(inventory['providers']), 'archived_leaves_verified': count,
            'original_leaves_verified': original_count, 'rows': 10,
            'official_awards_added': 0, 'scientific_rerun': False}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workspace', type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.workspace), indent=2))
