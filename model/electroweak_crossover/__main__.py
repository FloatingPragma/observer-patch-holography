"""Regenerate or replay the electroweak crossover receipt from the production record."""
import argparse
import json
import math
from pathlib import Path

from .certificate import certificate


def same(a, b, rel=1e-8) -> bool:
    """Equality of two JSON values with a relative tolerance on numbers."""
    if isinstance(a, bool) or isinstance(b, bool) or a is None or b is None or isinstance(a, str) or isinstance(b, str):
        return a == b
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        return math.isclose(a, b, rel_tol=rel, abs_tol=1e-12)
    if isinstance(a, list) and isinstance(b, list):
        return len(a) == len(b) and all(same(x, y, rel) for x, y in zip(a, b))
    if isinstance(a, dict) and isinstance(b, dict):
        return a.keys() == b.keys() and all(same(a[k], b[k], rel) for k in a)
    return False


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
    elif not same(json.loads(path.read_text()), value):
        raise SystemExit('electroweak crossover results differ from their recomputation')
    row = value['rows']['declared']
    print(f"PASS: {len(value['checks'])} checks: lattice gates, lattice-continuum relations, the reduction map, "
          f"the continuum extrapolation and the validation points; T_c = {row['T_c']:.2f} GeV for the declared chain")


if __name__ == '__main__':
    main()
