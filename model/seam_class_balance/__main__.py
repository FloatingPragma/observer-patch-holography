"""Regenerate or replay the class balance receipt."""
import argparse
import json
from pathlib import Path

from .certificate import certificate


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
    text = json.dumps(value, indent=1, sort_keys=True) + '\n'
    if args.write:
        path.write_text(text)
    else:
        old = json.loads(path.read_text())
        new = json.loads(text)
        if old['checks'] != new['checks'] or old['measured_points_sha256'] != new['measured_points_sha256']:
            raise SystemExit('class balance receipt differs from its recomputation')
    print(f"PASS: {len(value['checks'])} checks: one register, Cartan circle on complexes, flux densities, Wilson minimum, "
          'second reading, estimator, gates, fitted laws and sector laws')


if __name__ == '__main__':
    main()
