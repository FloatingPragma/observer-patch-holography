import argparse
import json
import numpy as np

from . import bridge


def equal(left, right, path=()):
    if type(left) is not type(right):
        return False
    if isinstance(left, dict):
        return set(left) == set(right) and all(equal(left[k], right[k], path+(k,)) for k in left)
    if isinstance(left, list):
        return len(left) == len(right) and all(equal(a, b, path) for a, b in zip(left, right, strict=True))
    if isinstance(left, float):
        tolerance = (2e-30+2e-8*max(abs(left), abs(right)) if any('probability' in key for key in path)
                     else 2e-8*max(1, abs(left), abs(right)))
        return np.isfinite(left) and np.isfinite(right) and abs(left-right) <= tolerance
    return left == right


def check_archive(stored, values):
    names = {'W', 'site_kernel', 'scalar_modes', 'raw_mode_kernel', 'canonical_mode_kernel'}
    if set(stored) != names:
        raise ArithmeticError('the complete source contact array inventory differs')
    for name in sorted(names):
        actual, expected = stored[name], values[name]
        if actual.shape != expected.shape or not np.isfinite(actual).all():
            raise ArithmeticError('a complete source contact array shape or value is invalid')
        if not np.allclose(actual, expected, atol=2e-8, rtol=0):
            raise ArithmeticError('a complete source contact array differs')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--write', action='store_true')
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    if args.write:
        bridge.write()
        print('Wrote the common defect/Higgs source and complete contact packet.')
    elif args.check:
        fresh = bridge.certificate()
        saved = json.loads((bridge.HERE/'results.json').read_text())
        if not equal(fresh, saved):
            raise ArithmeticError('the common source receipt differs on replay')
        with np.load(bridge.HERE/'operators.npz') as stored:
            check_archive(stored, bridge.spatial_contact())
        print('Reproduced the common defect/Higgs source and complete contact packet.')
    else:
        print(json.dumps(bridge.certificate(), indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
