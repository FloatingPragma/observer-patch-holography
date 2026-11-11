"""Regenerate or replay the quenched hadron binding receipt.

--write runs the fast gates, analyses the recorded tables and writes results.json.
--check does the same and compares with results.json: named checks and the hash of the tables exactly, numbers of
the analysis to a relative 1e-6. The numbers of the fast gates are held to their own tolerances by the named checks.
It does not rerun the production.
"""
import argparse
import json
from pathlib import Path

from .certificate import certificate, close

COMPARED = ('production', 'melting', 'production_gates', 'headline', 'measured', 'convention', 'literature', 'input_sha256')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    args = parser.parse_args()
    value = certificate()
    failed = [name for name, ok in value['checks'].items() if not ok]
    if failed:
        raise SystemExit('failed checks: ' + ', '.join(failed))
    path = Path(__file__).with_name('results.json')
    if args.write:
        path.write_text(json.dumps(value, indent=1, sort_keys=True) + '\n')
    else:
        stored = json.loads(path.read_text())
        value = json.loads(json.dumps(value))
        diff = []
        if stored['checks'] != value['checks']:
            diff.append('named checks differ')
        for key in COMPARED:
            diff += close(stored[key], value[key], path=key)
        if diff:
            raise SystemExit('quenched hadron binding differs from its recomputation: ' + '; '.join(diff[:5]))
    h = value['headline']
    print(f"PASS: {len(value['checks'])} checks: gamma algebra, gauge sector, Wilson-Dirac operator, hadron "
          'contractions, static potential, hadron masses, chiral limit and melting; '
          f"m_N / m_rho = {h['mN_over_mrho_chiral']['value']:.3f} +- {h['mN_over_mrho_chiral']['err']:.3f} in the chiral limit")


if __name__ == '__main__':
    main()
