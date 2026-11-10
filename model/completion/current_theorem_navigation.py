"""Current source navigation for an immutable, dated theorem-domain review.

Only the reviewed construction explanation may renew its documentary digest.
Every executable/Lean source and every theorem mapping retains the old review's
bytes or exact semantic value. This projection adds no theorem acceptance.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
from pathlib import Path
from tools import strict_json
from .family import binding

ROOT = Path(__file__).resolve().parents[2]
RETAINED = 'model/completion/theorem_map.json'
RETAINED_SHA256 = '4c64cec4dd0252f7b3d8c8e16b26a8498404cd1e9c496bc679555b33a92ed05d'
OUTPUT = 'model/completion/theorem_navigation_current.json'
DOCUMENT = 'model/CONSTRUCTION.md'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def build():
    if digest(ROOT / RETAINED) != RETAINED_SHA256:
        raise ValueError('dated theorem review changed')
    reviewed = strict_json.load(ROOT / RETAINED)
    current = copy.deepcopy(reviewed)
    for row in current['sources'].values():
        path = (ROOT / row['path']).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():
            raise ValueError('missing or escaping reviewed source')
        raw = path.read_bytes()
        actual = hashlib.sha256(raw).hexdigest()
        if actual != row['sha256'] and row['path'] != DOCUMENT:
            raise ValueError('reviewed scientific source changed: ' + row['path'])
        if any(anchor not in raw.decode() for anchor in row['reviewed_anchors']):
            raise ValueError('reviewed source anchor missing: ' + row['path'])
        row['sha256'] = actual
    live_campaign = binding()
    if {k: v for k, v in live_campaign.items() if k != 'sha256'} != {
            k: v for k, v in reviewed['campaign'].items() if k != 'sha256'}:
        raise ValueError('theorem review campaign identity changed')
    current['campaign'] = live_campaign
    return {'schema': 'oph.current_theorem_navigation.v1',
            'role': 'CURRENT_NAVIGATION_CUSTODY_ONLY',
            'retained_review': {'path': RETAINED, 'sha256': RETAINED_SHA256},
            'new_theorem_review_or_acceptance': False,
            'scientific_mappings_and_reference_domains_exactly_preserved': True,
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
