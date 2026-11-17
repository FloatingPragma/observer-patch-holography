"""Authenticate separate component theorem domains and scientific sources.

Scheduling and programme adoption are excluded from the mathematical map.
The official CCF-LT-2 dependency graph owns its source attachments.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
from pathlib import Path
from tools import strict_json

ROOT = Path(__file__).resolve().parents[2]
RETAINED = 'model/completion/theorem_map.json'
RETAINED_SHA256 = '4c64cec4dd0252f7b3d8c8e16b26a8498404cd1e9c496bc679555b33a92ed05d'
OUTPUT = 'model/completion/theorem_navigation_current.json'
DOCUMENT = 'model/CONSTRUCTION.md'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build():
    reviewed = strict_json.load(ROOT / RETAINED)
    current = copy.deepcopy(reviewed)
    # The scalar adoption row is programme metadata, not a mathematical
    # premise. The official acceptance graph owns source selection.
    current['mappings'] = [row for row in current['mappings'] if row['id'] != 'TM-11']
    current['sources'].pop('construction', None)
    for row in current['sources'].values():
        path = (ROOT / row['path']).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():
            raise ValueError('missing or escaping reviewed source')
        raw = path.read_bytes()
        actual = hashlib.sha256(raw).hexdigest()
        if any(anchor not in raw.decode() for anchor in row['reviewed_anchors']):
            raise ValueError('reviewed source anchor missing: ' + row['path'])
        row['sha256'] = actual
    # Campaign scheduling and documentary source selection are not mathematical
    # hypotheses of these component results. Authenticate the listed scientific
    # sources and retain their quantified domains directly.
    current.pop('campaign', None)
    current.pop('review_date', None)
    for row in current['mappings']:
        row.pop('owners', None)
    return {'schema': 'oph.current_theorem_navigation.v1',
            'role': 'CURRENT_NAVIGATION_CUSTODY_ONLY',
            'retained_review': {'path': RETAINED, 'sha256': RETAINED_SHA256},
            'new_theorem_review_or_acceptance': False,
            'scientific_mappings_and_reference_domains_exactly_preserved': True,
            'official_model': 'CCF-LT-2',
            'programme_adoption_row_is_excluded': True,
            'retained_documentary_pins_claimed_current': False,
            'reviewed_map_with_current_navigation': current}


def validate(packet=None):
    packet = strict_json.load(ROOT / OUTPUT) if packet is None else packet
    expected = build()
    if json.dumps(packet, sort_keys=True, allow_nan=False) != json.dumps(expected, sort_keys=True, allow_nan=False):
        raise ValueError('current theorem navigation changed outside admitted custody')
    return packet['reviewed_map_with_current_navigation']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--renew', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.renew:
        (ROOT / OUTPUT).write_text(json.dumps(build(), indent=2, sort_keys=True) + '\n')
    validate()
    print('Verified retained theorem mappings and current navigation; no theorem accepted')


if __name__ == '__main__':
    main()
