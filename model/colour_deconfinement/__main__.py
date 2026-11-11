"""Regenerate or replay the colour deconfinement receipt, or run a small demonstration scan."""
import argparse
import json
import math
from pathlib import Path

from .certificate import certificate


def same(a, b, rel=1e-9) -> bool:
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


def demo():
    from . import su3
    print('lattice 6^3 x 2, 30 + 60 compound sweeps per coupling')
    for beta in (4.8, 5.0, 5.1, 5.2, 5.4):
        _, loops, kernel = su3.series((6, 6, 6, 2), beta, 90, n_or=2, seed=int(beta * 100), start='hot')
        print(f'beta = {beta:.2f}   <|L|> = {abs(loops[30:]).mean():.4f}   ({kernel} kernel)')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--check', action='store_true')
    mode.add_argument('--demo', action='store_true')
    args = parser.parse_args()
    if args.demo:
        return demo()
    value = certificate()
    failed = [name for name, ok in value['checks'].items() if not ok]
    if failed:
        raise SystemExit('failed checks: ' + ', '.join(failed))
    path = Path(__file__).with_name('results.json')
    if args.write:
        path.write_text(json.dumps(value, indent=1, sort_keys=True) + '\n')
    else:
        stored = json.loads(path.read_text())
        fresh = json.loads(json.dumps(value))
        stored.pop('replay', None)
        fresh.pop('replay', None)
        if not same(stored, fresh):
            raise SystemExit('colour deconfinement receipt differs from its recomputation')
    print(f"PASS: {len(value['checks'])} checks: plaquette action, centre classes, centre transformation, "
          'strong coupling, kernel gates, recorded tables and scale')


if __name__ == '__main__':
    main()
