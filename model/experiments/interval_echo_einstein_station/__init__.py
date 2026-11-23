"""A finite observer measures a regional entropy response across a physical cut."""

import math
from numbers import Integral, Real


def _mismatches(actual, stored, path, rel, absolute):
    if isinstance(stored, dict):
        if not isinstance(actual, dict) or actual.keys() != stored.keys():
            yield f"{path}: keys differ"
            return
        for key, value in stored.items():
            yield from _mismatches(actual[key], value, f"{path}.{key}", rel, absolute)
    elif isinstance(stored, list):
        if not isinstance(actual, (list, tuple)) or len(actual) != len(stored):
            yield f"{path}: length differs"
            return
        for index, (left, right) in enumerate(zip(actual, stored)):
            yield from _mismatches(left, right, f"{path}[{index}]", rel, absolute)
    else:
        if isinstance(stored, float):
            equal = (isinstance(actual, Real) and not isinstance(actual, Integral)
                     and math.isfinite(stored) and math.isclose(float(actual), stored, rel_tol=rel, abs_tol=absolute))
        elif type(stored) is int:
            equal = isinstance(actual, Integral) and not isinstance(actual, bool) and actual == stored
        else:
            equal = type(actual) is type(stored) and actual == stored
        if not equal:
            yield f"{path}: {actual!r} != {stored!r}"


def assert_same_replay(actual, stored, rel=1e-9, absolute=1e-12):
    """Compare structures and typed counts exactly and allow rounding only on finite readouts."""
    found = list(_mismatches(actual, stored, "$", rel, absolute))
    if found:
        raise AssertionError(f"{len(found)} replay mismatch(es): " + "; ".join(found[:12]))
