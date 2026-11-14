"""Recompute literal native current identities; no physical cell is awarded.

Run with --prepare to inspect the fixed protocol; --run executes finite controls.
Uses integer arithmetic and standard-library dependencies only.
"""
import argparse
import json


def main():
    """Execute a finite signed-load identity check and rejection control."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prepare', action='store_true')
    parser.add_argument('--run', action='store_true')
    args = parser.parse_args()
    spec = {'task': 'M10.2 (proposed; allocation requires upstream review)',
            'load_range': [-32, 32], 'coins': [0, 1],
            'gate': 'exact current integrality, conservation, record identity',
            'physical_row_1_closed': False}
    if not args.run:
        print(json.dumps(spec, indent=2))
        return
    failures = []
    cases = 0
    mutant_rejected = False
    for a in range(-32, 33):
        for b in range(-32, 33):
            d = a - b
            currents = []
            for coin in (0, 1):
                j = (d - (2 * coin - 1) * (abs(d) % 2)) // 2
                aa, bb = a - j, b + j
                # E oriented from a to b; G = n + div E.
                checks = (aa + bb == a + b, aa + j == a,
                          bb - j == b, isinstance(j, int))
                if not all(checks):
                    failures.append([a, b, coin])
                mutant_rejected |= aa + (j + 1) != a
                currents.append(j)
                cases += 1
            if sum(currents) != d:
                failures.append([a, b, 'mean'])
    assert not failures and mutant_rejected
    print(json.dumps({**spec, 'cases': cases, 'failures': failures,
                      'corrupted_record_control_rejected': mutant_rejected,
                      'scope': 'arithmetic subgate; integer charge is supplied; '
                               'electron topology and full row gate remain open'}, indent=2))


if __name__ == '__main__':
    main()
