"""Read-only receipt comparison with bounded JSON-path diagnostics."""
from itertools import islice
import json


def difference_paths(expected, actual, path=''):
    """Yield JSON pointers to changed values, missing keys and list lengths."""
    if isinstance(expected, dict) and isinstance(actual, dict):
        for key in sorted(expected.keys() | actual.keys()):
            child = path + '/' + key.replace('~', '~0').replace('/', '~1')
            if key not in expected or key not in actual:
                yield child
            else:
                yield from difference_paths(expected[key], actual[key], child)
    elif isinstance(expected, list) and isinstance(actual, list):
        if len(expected) != len(actual):
            yield path
        for index, (left, right) in enumerate(zip(expected, actual)):
            yield from difference_paths(left, right, path + '/' + str(index))
    elif type(expected) is not type(actual) or expected != actual:
        yield path


def check_receipt(path, rendered, label):
    """Retain byte equality and report the first twenty changed JSON paths."""
    stored = path.read_text()
    if stored == rendered:
        return
    paths = list(islice(difference_paths(json.loads(stored), json.loads(rendered)), 20))
    raise SystemExit(json.dumps({'error': 'stale ' + label + ' receipt',
                                'difference_paths': paths,
                                'serialization_only': not paths}, sort_keys=True))
